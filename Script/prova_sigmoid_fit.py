import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Point this at whichever VI extraction you want to run phenology on -
# fluxnet_all_highlat_s2_indices.csv, fluxnet_all_highlat_hls_indices.csv, etc.
# The script only needs a site_id, a date column, and one or more of the VI
# columns below; missing columns are skipped automatically.
INPUT_CSV = DATA_DIR / "fluxnet_all_highlat_hls_indices.csv"
OUTPUT_CSV = DATA_DIR / "phenology_sos_eos_by_site_year_index.csv"

VI_COLUMNS = ['NDVI', 'EVI', 'NIRv', 'ARI1', 'CRI1', 'PRI']

# Minimum number of clear observations required in a phase (green-up or
# senescence) before attempting a sigmoid fit at all.
MIN_POINTS_PER_PHASE = 5
# Minimum valid_pixel_frac to keep an observation, if that column exists.
MIN_VALID_FRAC = 0.5
# Savitzky-Golay smoothing window (odd, in days on the interpolated daily grid).
SG_WINDOW = 15
SG_POLYORDER = 3
# Fractions defining the phenology thresholds, per the project's SOS/EOS convention.
F_LOW, F_HIGH = 0.10, 0.90

if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run the matching extraction script first.")


# ---------------------------------------------------------------------------
# 1. Logistic model and its inverse
# ---------------------------------------------------------------------------
def logistic(t, vmin, vmax, k, t0):
    """Generic 4-parameter logistic. k>0 = rising (green-up), k<0 = falling
    (senescence) - same functional form fits both phases."""
    return vmin + (vmax - vmin) / (1 + np.exp(-k * (t - t0)))


def invert_logistic(vmin, vmax, k, t0, frac):
    """DOY at which the logistic reaches `frac` of its (vmin -> vmax) amplitude."""
    frac = np.clip(frac, 1e-6, 1 - 1e-6)
    return t0 - (1.0 / k) * np.log((1.0 / frac) - 1.0)


def fit_phase(doy, vi, rising):
    """Fit the logistic to one phase. Returns (params, r2) or (None, None)
    if the fit fails or there isn't enough data."""
    if len(doy) < MIN_POINTS_PER_PHASE:
        return None, None

    vmin_obs, vmax_obs = np.nanmin(vi), np.nanmax(vi)
    amp = max(vmax_obs - vmin_obs, 1e-6)
    t0_guess = doy[np.argmin(np.abs(vi - (vmin_obs + amp / 2)))]
    k_guess = (2.0 / amp) if rising else -(2.0 / amp)

    p0 = [vmin_obs, vmax_obs, k_guess, t0_guess]
    bounds = (
        [vmin_obs - amp, vmax_obs - amp, -5 if not rising else 1e-4, doy.min() - 30],
        [vmin_obs + amp, vmax_obs + amp, -1e-4 if not rising else 5, doy.max() + 30]
    )

    try:
        popt, _ = curve_fit(logistic, doy, vi, p0=p0, bounds=bounds, maxfev=5000)
        pred = logistic(doy, *popt)
        ss_res = np.sum((vi - pred) ** 2)
        ss_tot = np.sum((vi - np.mean(vi)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        return popt, r2
    except (RuntimeError, ValueError):
        return None, None


def threshold_crossing_linear(doy, vi, frac_target, vmin, vmax, rising):
    """Fallback when the sigmoid fit fails: linear interpolation of the DOY
    at which the (normalized) VI crosses frac_target, using the raw points
    directly rather than a fitted curve."""
    norm = (vi - vmin) / max(vmax - vmin, 1e-6)
    order = np.argsort(doy)
    doy_s, norm_s = doy[order], norm[order]
    if not rising:
        # senescence: normalized VI decreases with DOY -> flip for np.interp,
        # which requires the x-array (norm_s here) to be increasing.
        doy_s, norm_s = doy_s[::-1], norm_s[::-1]
    if norm_s[0] > norm_s[-1]:
        return np.nan
    return float(np.interp(frac_target, norm_s, doy_s))


# ---------------------------------------------------------------------------
# 2. Per-phase extraction wrapper (fit, fall back to linear interp on failure)
# ---------------------------------------------------------------------------
def extract_phase_thresholds(doy, vi, rising, f_low, f_high):
    popt, r2 = fit_phase(doy, vi, rising)
    if popt is not None:
        vmin, vmax, k, t0 = popt
        t_low = invert_logistic(vmin, vmax, k, t0, f_low)
        t_high = invert_logistic(vmin, vmax, k, t0, f_high)
        method = 'sigmoid'
    else:
        vmin, vmax = np.nanmin(vi), np.nanmax(vi)
        t_low = threshold_crossing_linear(doy, vi, f_low, vmin, vmax, rising)
        t_high = threshold_crossing_linear(doy, vi, f_high, vmin, vmax, rising)
        method, r2 = 'linear_fallback', np.nan
    return t_low, t_high, method, r2


# ---------------------------------------------------------------------------
# 3. Smoothing + peak detection on an irregular time series
# ---------------------------------------------------------------------------
def smooth_and_find_peak(doy, vi):
    """Interpolate onto a daily grid, Savitzky-Golay smooth, and return the
    DOY of the seasonal maximum - used only to split green-up vs senescence,
    the actual threshold fits still use the original (unsmoothed) points."""
    daily_doy = np.arange(doy.min(), doy.max() + 1)
    daily_vi = np.interp(daily_doy, doy, vi)

    window = min(SG_WINDOW, len(daily_doy) - (1 - len(daily_doy) % 2))
    if window >= SG_POLYORDER + 2 and window % 2 == 1 and window <= len(daily_doy):
        smoothed = savgol_filter(daily_vi, window_length=window, polyorder=SG_POLYORDER)
    else:
        smoothed = daily_vi  # too few points to smooth meaningfully

    return daily_doy[np.argmax(smoothed)]


# ---------------------------------------------------------------------------
# 4. Main loop: site x year x VI index
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['date'])
df['year'] = df['date'].dt.year
df['doy'] = df['date'].dt.dayofyear

if 'valid_pixel_frac' in df.columns:
    df = df[df['valid_pixel_frac'].isna() | (df['valid_pixel_frac'] >= MIN_VALID_FRAC)]

available_vi = [c for c in VI_COLUMNS if c in df.columns]
print(f"Fitting phenology for indices: {available_vi}")

results = []
site_years = df.groupby(['site_id', 'year'])
total = site_years.ngroups

for gi, ((site_id, year), group) in enumerate(site_years, start=1):
    for vi_col in available_vi:
        sub = group[['doy', vi_col]].dropna().drop_duplicates(subset='doy').sort_values('doy')
        if len(sub) < 2 * MIN_POINTS_PER_PHASE:
            continue

        doy_all = sub['doy'].to_numpy(dtype=float)
        vi_all = sub[vi_col].to_numpy(dtype=float)

        peak_doy = smooth_and_find_peak(doy_all, vi_all)

        greenup_mask = doy_all <= peak_doy
        senesc_mask = doy_all >= peak_doy

        sos10, sos90, gu_method, gu_r2 = extract_phase_thresholds(
            doy_all[greenup_mask], vi_all[greenup_mask], rising=True, f_low=F_LOW, f_high=F_HIGH)
        eos90, eos10, se_method, se_r2 = extract_phase_thresholds(
            doy_all[senesc_mask], vi_all[senesc_mask], rising=False, f_low=F_HIGH, f_high=F_LOW)

        gsl = eos10 - sos90 if np.isfinite(eos10) and np.isfinite(sos90) else np.nan

        results.append({
            'site_id': site_id, 'year': year, 'vi_index': vi_col,
            'n_obs': len(sub), 'peak_doy': peak_doy,
            'SOS10': sos10, 'SOS90': sos90, 'EOS90': eos90, 'EOS10': eos10,
            'growing_season_length': gsl,
            'greenup_method': gu_method, 'greenup_r2': gu_r2,
            'senescence_method': se_method, 'senescence_r2': se_r2
        })

    if gi % 25 == 0 or gi == total:
        print(f"  processed {gi}/{total} site-years")

results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_CSV, index=False)

n_sigmoid = (results_df[['greenup_method', 'senescence_method']] == 'sigmoid').sum().sum()
n_fallback = (results_df[['greenup_method', 'senescence_method']] == 'linear_fallback').sum().sum()
print(f"\nDone. {len(results_df)} site-year-index rows written to '{OUTPUT_CSV}'.")
print(f"Phase fits: {n_sigmoid} sigmoid, {n_fallback} linear-interpolation fallback.")