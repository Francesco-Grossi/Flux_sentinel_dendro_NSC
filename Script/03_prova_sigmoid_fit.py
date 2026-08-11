import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Two SEPARATE output folders: one with the fitted curve drawn on top of the
# points, one with raw points only. Same site-years, same file basenames,
# so the two folders are easy to flip between when checking a fit visually.
FIT_PLOT_DIR = Path(__file__).resolve().parent.parent / "figure" / "sigmoid_fit"
RAW_PLOT_DIR = Path(__file__).resolve().parent.parent / "figure" / "sigmoid_raw"
FIT_PLOT_DIR.mkdir(parents=True, exist_ok=True)
RAW_PLOT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_CSV = DATA_DIR / "fluxnet_all_highlat_hls_indices.csv"
OUTPUT_CSV = DATA_DIR / "phenology_sos_eos_by_site_year_index.csv"

VI_COLUMNS = ['NDVI', 'EVI', 'NIRv', 'ARI1', 'CRI1', 'PRI']
# Indices with a well-defined physical range, checked/clipped before fitting.
# ARI1/CRI1 are reciprocal indices (1/reflectance terms) with no fixed
# physical bound, so they get their own outlier handling below instead.
PHYSICAL_BOUNDS = {'NDVI': (-1, 1), 'EVI': (-1, 1), 'NIRv': (-1, 1), 'PRI': (-1, 1)}

MIN_POINTS_FOR_FIT = 12   # double-logistic has 6 free params; need real margin over that
MIN_VALID_FRAC = 0.5      # QC floor on Fmask-derived clear-pixel fraction
SG_WINDOW = 15
SG_POLYORDER = 3
F_LOW, F_HIGH = 0.10, 0.90

if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run the matching extraction script first.")


# ---------------------------------------------------------------------------
# 1. Double-logistic model (Zhang et al. 2003 / Beck et al. 2006 / TIMESAT)
# ---------------------------------------------------------------------------
# One continuous curve for the whole season instead of two separately-fit
# sigmoids glued together at a hand-picked peak DOY. S/mS control the
# green-up inflection point and rate, A/mA control the senescence inflection
# and rate. This also removes the need for a hardcoded EOS_MIN_DOY cutoff:
# the model doesn't assume which calendar direction is "green-up" vs
# "senescence", so it works the same way for a boreal site (peak in summer)
# and a Mediterranean site (peak in winter/spring, e.g. ES-Abr).
def double_logistic(t, vmin, vmax, S, mS, A, mA):
    # Clip the exponent argument to avoid float overflow warnings for
    # points far from the inflection point; doesn't change the result
    # (exp of a huge/very-negative number saturates to 0 either way).
    z_S = np.clip(-mS * (t - S), -500, 500)
    z_A = np.clip(-mA * (t - A), -500, 500)
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
    # guess (mS = mA = 0.15) works whether the index ranges over 0-1 (NDVI)
    # or ~1e-5 (the currently mis-scaled EVI in this extraction). Scaling by
    # 1/amp, as an earlier version of this did, sends the guess to the
    # thousands for tiny-amplitude indices, which lands outside the bounds
    # below and makes least_squares raise before it even starts.
    rate_guess = 0.15

    p0 = [
        vmin_obs, vmax_obs,
        max(peak_guess - 20, trough_guess_lo - 60),   # S: green-up inflection guess
        rate_guess,                                    # mS > 0
        min(peak_guess + 20, trough_guess_hi + 60),    # A: senescence inflection guess
        rate_guess,                                     # mA > 0 (sign handled inside the model)
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


def threshold_crossing_linear(doy, vi, frac_target, vmin, vmax, rising):
    """Fallback: linear interpolation directly on the raw QC'd points,
    used only when the double-logistic fit doesn't have enough points
    or fails to converge."""
    norm = (vi - vmin) / max(vmax - vmin, 1e-6)
    order = np.argsort(doy)
    doy_s, norm_s = doy[order], norm[order]
    if not rising:
        doy_s, norm_s = doy_s[::-1], norm_s[::-1]
    if norm_s[0] > norm_s[-1]:
        return np.nan
    return float(np.interp(frac_target, norm_s, doy_s))


def smooth_and_find_peak(doy, vi):
    """Used only by the linear fallback path to split green-up/senescence."""
    daily_doy = np.arange(doy.min(), doy.max() + 1)
    daily_vi = np.interp(daily_doy, doy, vi)

    window = min(SG_WINDOW, len(daily_doy) - (1 - len(daily_doy) % 2))
    if window >= SG_POLYORDER + 2 and window % 2 == 1 and window <= len(daily_doy):
        smoothed = savgol_filter(daily_vi, window_length=window, polyorder=SG_POLYORDER)
    else:
        smoothed = daily_vi

    return daily_doy[np.argmax(smoothed)]


# ---------------------------------------------------------------------------
# 2. Load + QC
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['date'])
df['year'] = df['date'].dt.year
df['doy'] = df['date'].dt.dayofyear

# QC floor on Fmask clear-pixel fraction. Rows with no valid_pixel_frac info
# are kept (older extractions may not have this column); everything else
# below the floor is dropped outright rather than fit on cloud/shadow noise.
if 'valid_pixel_frac' in df.columns:
    df = df[df['valid_pixel_frac'].isna() | (df['valid_pixel_frac'] >= MIN_VALID_FRAC)]
    df['fit_weight'] = df['valid_pixel_frac'].fillna(1.0).clip(lower=0.05)
else:
    df['fit_weight'] = 1.0

# Physically-bounded indices: clip anything outside [-1, 1] to NaN. About 2%
# of NDVI rows in the raw extraction fall outside this range (up to ~23),
# almost all at very low valid_pixel_frac - residual cloud/shadow/snow
# contamination that Fmask didn't fully catch.
for vi_col, (lo, hi) in PHYSICAL_BOUNDS.items():
    if vi_col in df.columns:
        df.loc[(df[vi_col] < lo) | (df[vi_col] > hi), vi_col] = np.nan

# ARI1/CRI1 are reciprocal indices (1/reflectance terms) that blow up when
# the denominator reflectance is near zero - raw values up to ~1.9e6 are
# expected, not a bug. Drop the extreme raw outliers, then rescale each
# site-year robustly (2nd-98th percentile) onto [0, 1] so magnitudes are
# comparable across site-years and usable by the fit.
for vi_col in ['ARI1', 'CRI1']:
    if vi_col in df.columns:
        df.loc[(df[vi_col] < -100) | (df[vi_col] > 100), vi_col] = np.nan

        def scale_group(g):
            valid = g.dropna()
            if len(valid) < 3:
                return g
            p2, p98 = np.nanpercentile(valid, [2, 98])
            if np.isfinite(p2) and np.isfinite(p98) and p98 > p2:
                return np.clip((g - p2) / (p98 - p2), 0, 1)
            return g

        df[vi_col] = df.groupby(['site_id', 'year'])[vi_col].transform(scale_group)

# NOTE: EVI in this extraction is ~1e4x smaller than physically expected
# (mean ~0.00004 vs a normal ~0.3-0.5). That traces back to the EVI formula
# in the HLS extraction script, where the "+1" denominator constant is
# added after reflectance was already scaled to 0-1, swamping the numerator.
# It's still internally self-consistent (relative shape/seasonality is
# fine), so the fit below still works on it, but absolute EVI values here
# should not be compared to literature EVI values until that's fixed
# upstream in the extraction script.

available_vi = [c for c in VI_COLUMNS if c in df.columns]
print(f"Fitting phenology for indices: {available_vi}")


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

        if popt is not None:
            t_grid = np.linspace(doy_all.min() - 10, doy_all.max() + 10, 2000)
            y_grid = double_logistic(t_grid, *popt)

            sos10 = crossing_doy(t_grid, y_grid, F_LOW, rising=True)
            sos90 = crossing_doy(t_grid, y_grid, F_HIGH, rising=True)
            eos90 = crossing_doy(t_grid, y_grid, F_HIGH, rising=False)
            eos10 = crossing_doy(t_grid, y_grid, F_LOW, rising=False)
            method = 'double_logistic'
        else:
            # Fallback: naive peak-split + linear threshold crossing on raw points.
            peak_doy = smooth_and_find_peak(doy_all, vi_all)
            greenup_mask = doy_all <= peak_doy
            senesc_mask = doy_all >= peak_doy

            vmin_gu = np.nanmin(vi_all[greenup_mask]) if greenup_mask.any() else np.nan
            vmax_gu = np.nanmax(vi_all[greenup_mask]) if greenup_mask.any() else np.nan
            vmin_se = np.nanmin(vi_all[senesc_mask]) if senesc_mask.any() else np.nan
            vmax_se = np.nanmax(vi_all[senesc_mask]) if senesc_mask.any() else np.nan

            sos10 = threshold_crossing_linear(doy_all[greenup_mask], vi_all[greenup_mask], F_LOW, vmin_gu, vmax_gu, True)
            sos90 = threshold_crossing_linear(doy_all[greenup_mask], vi_all[greenup_mask], F_HIGH, vmin_gu, vmax_gu, True)
            eos90 = threshold_crossing_linear(doy_all[senesc_mask], vi_all[senesc_mask], F_HIGH, vmin_se, vmax_se, False)
            eos10 = threshold_crossing_linear(doy_all[senesc_mask], vi_all[senesc_mask], F_LOW, vmin_se, vmax_se, False)
            method, r2, rmse, corr = 'linear_fallback', np.nan, np.nan, np.nan
            popt = (np.nan,) * 6

        gsl = eos10 - sos90 if np.isfinite(eos10) and np.isfinite(sos90) else np.nan

        results.append({
            'site_id': site_id, 'year': year, 'vi_index': vi_col,
            'n_obs': len(sub), 'method': method,
            'r2': r2, 'rmse': rmse, 'corr': corr,
            'SOS10': sos10, 'SOS90': sos90, 'EOS90': eos90, 'EOS10': eos10,
            'growing_season_length': gsl,
            'vmin': popt[0], 'vmax': popt[1], 'S': popt[2], 'mS': popt[3], 'A': popt[4], 'mA': popt[5],
        })

    if gi % 25 == 0 or gi == total:
        print(f"  processed {gi}/{total} site-years")

results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_CSV, index=False)

n_fit = (results_df['method'] == 'double_logistic').sum()
n_fallback = (results_df['method'] == 'linear_fallback').sum()
print(f"\nDone. {len(results_df)} site-year-index rows written to '{OUTPUT_CSV}'.")
print(f"Fits: {n_fit} double-logistic, {n_fallback} linear-interpolation fallback.")


# ---------------------------------------------------------------------------
# 4. Plots - two separate folders: raw only, and raw + fitted curve
# ---------------------------------------------------------------------------
target_indices = ['NDVI', 'EVI', 'ARI1', 'CRI1']
unique_sites = df['site_id'].unique()


def make_plot(site_id, with_fit):
    site_df = df[df['site_id'] == site_id]
    site_results = results_df[results_df['site_id'] == site_id]
    years = sorted(site_df['year'].dropna().unique())
    if not years:
        return

    colors = plt.cm.viridis(np.linspace(0, 1, max(1, len(years))))
    fig, axes = plt.subplots(len(target_indices), 1, figsize=(10, 4 * len(target_indices)), sharex=True)
    if len(target_indices) == 1:
        axes = [axes]

    title = "Vegetation Index Values + Double-Logistic Fit" if with_fit else "Raw Vegetation Index Values"
    fig.suptitle(f"{title} - {site_id}", fontsize=18, fontweight='bold')

    for i, vi_col in enumerate(target_indices):
        ax = axes[i]
        if vi_col not in site_df.columns:
            ax.set_visible(False)
            continue

        y_label = f"{vi_col} (Normalized)" if vi_col in ['ARI1', 'CRI1'] else f"{vi_col} Value"
        ax.set_ylabel(y_label, fontsize=12)
        ax.set_title(vi_col, fontsize=14)

        for j, year in enumerate(years):
            color = colors[j]
            raw_sub = site_df[site_df['year'] == year][['doy', vi_col]].dropna().sort_values('doy')
            if raw_sub.empty:
                continue

            doy_all = raw_sub['doy'].to_numpy(float)
            vi_all = raw_sub[vi_col].to_numpy(float)

            if not with_fit:
                ax.scatter(doy_all, vi_all, color=color, s=20, alpha=0.6, label=str(int(year)))
                continue

            row = site_results[(site_results['year'] == year) & (site_results['vi_index'] == vi_col)]
            if row.empty:
                ax.scatter(doy_all, vi_all, color=color, s=20, alpha=0.6, label=str(int(year)))
                continue
            row = row.iloc[0]

            ax.scatter(doy_all, vi_all, color=color, s=20, alpha=0.4)
            if row['method'] == 'double_logistic':
                t_fine = np.linspace(doy_all.min(), doy_all.max(), 300)
                y_fine = double_logistic(t_fine, row['vmin'], row['vmax'], row['S'], row['mS'], row['A'], row['mA'])
                label = f"{int(year)} (RMSE: {row['rmse']:.3f}, r: {row['corr']:.2f})"
                ax.plot(t_fine, y_fine, color=color, linewidth=2, label=label)
            else:
                ax.scatter([], [], color=color, label=f"{int(year)} (fallback, no smooth fit)")

        ax.legend(fontsize=8, loc="best", ncol=2)
        ax.grid(True, linestyle='--', alpha=0.6)

    axes[-1].set_xlabel("Day of Year", fontsize=12)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    out_dir = FIT_PLOT_DIR if with_fit else RAW_PLOT_DIR
    suffix = "fit" if with_fit else "raw"
    plt.savefig(out_dir / f"{site_id}_{suffix}_values.png", dpi=150)
    plt.close(fig)


print("\nGenerating raw-only plots...")
for site_id in unique_sites:
    make_plot(site_id, with_fit=False)
print(f"Saved raw-only plots to '{RAW_PLOT_DIR}'.")

print("\nGenerating raw + fit plots...")
for site_id in unique_sites:
    make_plot(site_id, with_fit=True)
print(f"Saved raw + fit plots to '{FIT_PLOT_DIR}'.")