import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

INPUT_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"   # output of the FLUXNET + Landsat merge script
OUTPUT_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"

# Landsat (HLSL30) has no red-edge band, so only these three VIs exist in
# the merged file - no ARI1/CRI1/PRI here (those need Sentinel-2).
VI_COLUMNS = ['NDVI', 'EVI', 'NIRv']
PHYSICAL_BOUNDS = {'NDVI': (-1, 1), 'EVI': (-1, 1), 'NIRv': (-1, 1)}

MIN_POINTS_FOR_FIT = 12   # double-logistic has 6 free params; need real margin over that
MIN_VALID_FRAC = 0.5      # QC floor on Fmask-derived clear-pixel fraction

# Percentiles for each phase, as % of the seasonal amplitude (vmin->vmax):
# leaf_out_X = curve has risen TO X% (SOS-X convention); EOS_X = curve is
# AT X% on the falling limb (standard Zhang/TIMESAT convention - EOS90 is
# early decline/still mostly green, EOS10 is late decline/near dormancy).
PERCENTILES = [10, 50, 90]

if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run the FLUXNET + Landsat merge script first.")


# ---------------------------------------------------------------------------
# 1. Double-logistic model (Zhang et al. 2003 / Beck et al. 2006 / TIMESAT)
# ---------------------------------------------------------------------------
# One continuous curve for the whole season instead of two separately-fit
# sigmoids glued together at a hand-picked peak DOY. S/i_sos control the
# green-up (start-of-season) inflection point and rate; A/i_eos control the
# senescence (end-of-season) inflection and rate. 'i' is the kinetic
# (steepness) parameter requested for each phase - i_sos for green-up,
# i_eos for senescence - both in units of 1/day: a larger i means a
# faster, more abrupt transition, a smaller i a slower, more gradual one.
def double_logistic(t, vmin, vmax, S, i_sos, A, i_eos):
    # Clip the exponent argument to avoid float overflow warnings for
    # points far from the inflection point; doesn't change the result
    # (exp of a huge/very-negative number saturates to 0 either way).
    z_S = np.clip(-i_sos * (t - S), -500, 500)
    z_A = np.clip(-i_eos * (t - A), -500, 500)
    return vmin + (vmax - vmin) * (1.0 / (1.0 + np.exp(z_S)) - 1.0 / (1.0 + np.exp(z_A)))


def _residuals(params, t, y, w):
    return w * (double_logistic(t, *params) - y)


def fit_double_logistic(doy, vi, weights):
    """Robust (soft-L1 loss) weighted fit of the double logistic.
    Returns (popt, r2, rmse, corr) or (None, None, None, None)."""
    if len(doy) < MIN_POINTS_FOR_FIT:
        return None, None, None, None

    vmin_obs, vmax_obs = np.nanmin(vi), np.nanmax(vi)
    amp = max(vmax_obs - vmin_obs, 1e-6)

    peak_guess = doy[np.argmax(vi)]
    trough_guess_lo, trough_guess_hi = doy.min(), doy.max()

    # Rate-parameter guess is in units of 1/day and must NOT scale with the
    # index's amplitude - a fixed "typical transition takes ~2-3 weeks"
    # guess (i = 0.15) works whether the index ranges over 0-1 (NDVI) or a
    # very different scale for another index.
    rate_guess = 0.15

    p0 = [
        vmin_obs, vmax_obs,
        max(peak_guess - 20, trough_guess_lo - 60),   # S: green-up inflection guess
        rate_guess,                                    # i_sos > 0
        min(peak_guess + 20, trough_guess_hi + 60),    # A: senescence inflection guess
        rate_guess,                                     # i_eos > 0 (sign handled inside the model)
    ]
    lower = [vmin_obs - amp, vmax_obs - amp, doy.min() - 60, 1e-4, doy.min() - 60, 1e-4]
    upper = [vmin_obs + amp, vmax_obs + amp, doy.max() + 60, 5.0, doy.max() + 60, 5.0]
    # Defensive clip: guarantees p0 is inside (lower, upper) even in edge
    # cases (e.g. very short doy spans) so least_squares never raises on
    # "initial guess outside bounds".
    p0 = list(np.clip(p0, lower, upper))

    try:
        result = least_squares(
            _residuals, p0, args=(doy, vi, weights),
            bounds=(lower, upper), loss='soft_l1', f_scale=0.1 * amp, max_nfev=5000
        )
        popt = result.x
        pred = double_logistic(doy, *popt)
        ss_res = np.sum((vi - pred) ** 2)
        ss_tot = np.sum((vi - np.mean(vi)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rmse = np.sqrt(np.mean((vi - pred) ** 2))
        corr = np.corrcoef(vi, pred)[0, 1] if np.std(pred) > 0 else np.nan
        return popt, r2, rmse, corr
    except (RuntimeError, ValueError):
        return None, None, None, None


def crossing_doy(t_grid, y_grid, frac_target, rising):
    """DOY at which the fitted curve crosses `frac_target` of its local
    (min->max) amplitude, on the rising (green-up) or falling (senescence)
    side of the curve's peak. Works on the smooth fitted curve rather than
    the noisy raw points, and doesn't assume a fixed calendar direction."""
    peak_idx = int(np.argmax(y_grid))
    if rising:
        t_seg, y_seg = t_grid[:peak_idx + 1], y_grid[:peak_idx + 1]
    else:
        t_seg, y_seg = t_grid[peak_idx:], y_grid[peak_idx:]

    if len(t_seg) < 2:
        return np.nan

    vmin_seg, vmax_seg = y_seg.min(), y_seg.max()
    norm = (y_seg - vmin_seg) / max(vmax_seg - vmin_seg, 1e-9)
    if not rising:
        t_seg, norm = t_seg[::-1], norm[::-1]
    if norm[0] > norm[-1]:
        return np.nan
    return float(np.interp(frac_target, norm, t_seg))


def segment_fit_metrics(doy_obs, vi_obs, popt, peak_doy, rising):
    """r2/rmse/corr of the fitted curve against the raw QC'd points,
    restricted to one phase (green-up: doy <= peak; senescence: doy >=
    peak) - i.e. phase-specific fit quality, separate from the whole-season
    r2/rmse/corr already returned by fit_double_logistic."""
    mask = (doy_obs <= peak_doy) if rising else (doy_obs >= peak_doy)
    if mask.sum() < 3:
        return np.nan, np.nan, np.nan

    t_seg, y_seg = doy_obs[mask], vi_obs[mask]
    pred_seg = double_logistic(t_seg, *popt)

    ss_res = np.sum((y_seg - pred_seg) ** 2)
    ss_tot = np.sum((y_seg - np.mean(y_seg)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = np.sqrt(np.mean((y_seg - pred_seg) ** 2))
    corr = (np.corrcoef(y_seg, pred_seg)[0, 1]
            if len(y_seg) > 1 and np.std(pred_seg) > 0 and np.std(y_seg) > 0 else np.nan)
    return r2, rmse, corr


# ---------------------------------------------------------------------------
# 2. Load + QC
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
df = df.dropna(subset=['date']).copy()
df['year'] = df['date'].dt.year
df['doy'] = df['date'].dt.dayofyear

# Only rows with an actual Landsat observation are usable for the curve fit
# (most FLUXNET days have no matching overpass - that's the sparse/dense
# join from the merge script working as intended).
if 'has_landsat_obs' in df.columns:
    df = df[df['has_landsat_obs']].copy()

# QC floor on Fmask clear-pixel fraction. Rows with no valid_pixel_frac info
# are kept (older extractions may not have this column); everything else
# below the floor is dropped outright rather than fit on cloud/shadow noise.
if 'valid_pixel_frac' in df.columns:
    df = df[df['valid_pixel_frac'].isna() | (df['valid_pixel_frac'] >= MIN_VALID_FRAC)]
    df['fit_weight'] = df['valid_pixel_frac'].fillna(1.0).clip(lower=0.05)
else:
    df['fit_weight'] = 1.0

# Physically-bounded indices: clip anything outside [-1, 1] to NaN.
for vi_col, (lo, hi) in PHYSICAL_BOUNDS.items():
    if vi_col in df.columns:
        df.loc[(df[vi_col] < lo) | (df[vi_col] > hi), vi_col] = np.nan

available_vi = [c for c in VI_COLUMNS if c in df.columns]
print(f"Fitting double-logistic phenology curves for indices: {available_vi}")


# ---------------------------------------------------------------------------
# 3. Main loop: site x year x VI index
# ---------------------------------------------------------------------------
results = []
site_years = df.groupby(['site_id', 'year'])
total = site_years.ngroups

for gi, ((site_id, year), group) in enumerate(site_years, start=1):
    for vi_col in available_vi:
        sub = group[['doy', vi_col, 'fit_weight']].dropna(subset=['doy', vi_col])
        sub = sub.drop_duplicates(subset='doy').sort_values('doy')
        if len(sub) < MIN_POINTS_FOR_FIT:
            continue

        doy_all = sub['doy'].to_numpy(dtype=float)
        vi_all = sub[vi_col].to_numpy(dtype=float)
        w_all = sub['fit_weight'].to_numpy(dtype=float)

        popt, r2, rmse, corr = fit_double_logistic(doy_all, vi_all, w_all)

        row = {'site_id': site_id, 'year': year, 'vi_index': vi_col, 'n_obs': len(sub)}

        if popt is not None:
            t_grid = np.linspace(doy_all.min() - 10, doy_all.max() + 10, 2000)
            y_grid = double_logistic(t_grid, *popt)
            peak_doy = float(t_grid[np.argmax(y_grid)])

            # Green-up (leaf-out) percentiles: DOY at which the curve has
            # risen to X% of its seasonal amplitude (SOS-X convention).
            leafout = {p: crossing_doy(t_grid, y_grid, p / 100.0, rising=True) for p in PERCENTILES}
            # EOS-X percentiles (standard Zhang/TIMESAT convention): DOY on
            # the falling limb where the curve's VALUE is X% of the way
            # from vmin to vmax - i.e. EOS90 is where the curve is STILL AT
            # 90% (early decline, canopy still mostly green), and EOS10 is
            # where the curve has fallen DOWN TO 10% (late decline, near
            # dormancy). Same crossing_doy(frac_target=X/100) call as
            # leaf-out, just on the falling side - no inversion.
            eos = {p: crossing_doy(t_grid, y_grid, p / 100.0, rising=False) for p in PERCENTILES}

            gu_r2, gu_rmse, gu_corr = segment_fit_metrics(doy_all, vi_all, popt, peak_doy, rising=True)
            se_r2, se_rmse, se_corr = segment_fit_metrics(doy_all, vi_all, popt, peak_doy, rising=False)

            row.update({
                'method': 'double_logistic',
                'r2': r2, 'rmse': rmse, 'corr': corr,
                'vmin': popt[0], 'vmax': popt[1],
                'S': popt[2], 'greenup_kinetic_i': popt[3],
                'A': popt[4], 'senescence_kinetic_i': popt[5],
                'greenup_corr': gu_corr, 'greenup_rmse': gu_rmse,
                'leaf_out_10': leafout[10], 'leaf_out_50': leafout[50], 'leaf_out_90': leafout[90],
                'senescence_corr': se_corr, 'senescence_rmse': se_rmse,
                'EOS10': eos[10], 'EOS50': eos[50], 'EOS90': eos[90],
            })
        else:
            # No double-logistic curve could be fit (too few points or the
            # optimizer didn't converge) - report NaN across the board
            # rather than estimating percentile crossings some other way.
            # A crossing DOY derived without an actual fitted curve isn't
            # comparable to the ones that are, so it's not reported at all.
            row.update({
                'method': 'fit_failed',
                'r2': np.nan, 'rmse': np.nan, 'corr': np.nan,
                'vmin': np.nan, 'vmax': np.nan,
                'S': np.nan, 'greenup_kinetic_i': np.nan,
                'A': np.nan, 'senescence_kinetic_i': np.nan,
                'greenup_corr': np.nan, 'greenup_rmse': np.nan,
                'leaf_out_10': np.nan, 'leaf_out_50': np.nan, 'leaf_out_90': np.nan,
                'senescence_corr': np.nan, 'senescence_rmse': np.nan,
                'EOS10': np.nan, 'EOS50': np.nan, 'EOS90': np.nan,
            })

        results.append(row)

    if gi % 25 == 0 or gi == total:
        print(f"  processed {gi}/{total} site-years")

results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_CSV, index=False)

n_fit = (results_df['method'] == 'double_logistic').sum()
n_failed = (results_df['method'] == 'fit_failed').sum()
print(f"\nDone. {len(results_df)} site-year-index rows written to '{OUTPUT_CSV}'.")
print(f"Fits: {n_fit} double-logistic succeeded, {n_failed} failed "
      "(too few points or non-convergent - all metrics NaN for those rows).")