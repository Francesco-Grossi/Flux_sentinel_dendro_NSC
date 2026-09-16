"""
PIPELINE STEP 5 - Fit a double-logistic phenology curve (Zhang/Beck/TIMESAT)
to each site-year x VI-index (NDVI, EVI, NIRv) series, and extract leaf_out
(green-up) and EOS (senescence) percentile-crossing DOYs plus the kinetic
(steepness) parameters for both phases.

Output: data/phenology_double_logistic_by_site_year_index.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
INPUT_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"
OUTPUT_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"

VI_COLUMNS = ['NDVI', 'EVI', 'NIRv']
PHYSICAL_BOUNDS = {'NDVI': (-1, 1), 'EVI': (-1, 1), 'NIRv': (-1, 1)}
MIN_POINTS_FOR_FIT = 12
MIN_VALID_FRAC = 0.5
PERCENTILES = [10, 50, 90]

if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run 04_merge_fluxnet_landsat.py first.")


def double_logistic(t, vmin, vmax, S, i_sos, A, i_eos):
    z_S = np.clip(-i_sos * (t - S), -500, 500)
    z_A = np.clip(-i_eos * (t - A), -500, 500)
    return vmin + (vmax - vmin) * (1.0 / (1.0 + np.exp(z_S)) - 1.0 / (1.0 + np.exp(z_A)))


def _residuals(params, t, y, w):
    return w * (double_logistic(t, *params) - y)


def fit_double_logistic(doy, vi, weights):
    if len(doy) < MIN_POINTS_FOR_FIT:
        return None, None, None, None
    vmin_obs, vmax_obs = np.nanmin(vi), np.nanmax(vi)
    amp = max(vmax_obs - vmin_obs, 1e-6)
    peak_guess = doy[np.argmax(vi)]
    trough_lo, trough_hi = doy.min(), doy.max()
    rate_guess = 0.15
    p0 = [vmin_obs, vmax_obs, max(peak_guess - 20, trough_lo - 60), rate_guess,
          min(peak_guess + 20, trough_hi + 60), rate_guess]
    lower = [vmin_obs - amp, vmax_obs - amp, doy.min() - 60, 1e-4, doy.min() - 60, 1e-4]
    upper = [vmin_obs + amp, vmax_obs + amp, doy.max() + 60, 5.0, doy.max() + 60, 5.0]
    p0 = list(np.clip(p0, lower, upper))
    try:
        result = least_squares(_residuals, p0, args=(doy, vi, weights), bounds=(lower, upper),
                                loss='soft_l1', f_scale=0.1 * amp, max_nfev=5000)
        popt = result.x
        pred = double_logistic(doy, *popt)
        ss_res, ss_tot = np.sum((vi - pred) ** 2), np.sum((vi - np.mean(vi)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rmse = np.sqrt(np.mean((vi - pred) ** 2))
        corr = np.corrcoef(vi, pred)[0, 1] if np.std(pred) > 0 else np.nan
        return popt, r2, rmse, corr
    except (RuntimeError, ValueError):
        return None, None, None, None


def crossing_doy(t_grid, y_grid, frac_target, rising):
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


df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
df = df.dropna(subset=['date']).copy()
df['year'] = df['date'].dt.year
df['doy'] = df['date'].dt.dayofyear

if 'has_landsat_obs' in df.columns:
    df = df[df['has_landsat_obs']].copy()
if 'valid_pixel_frac' in df.columns:
    df = df[df['valid_pixel_frac'].isna() | (df['valid_pixel_frac'] >= MIN_VALID_FRAC)]
    df['fit_weight'] = df['valid_pixel_frac'].fillna(1.0).clip(lower=0.05)
else:
    df['fit_weight'] = 1.0

for vi_col, (lo, hi) in PHYSICAL_BOUNDS.items():
    if vi_col in df.columns:
        df.loc[(df[vi_col] < lo) | (df[vi_col] > hi), vi_col] = np.nan

available_vi = [c for c in VI_COLUMNS if c in df.columns]
print(f"Fitting double-logistic phenology curves for: {available_vi}")

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
            leafout = {p: crossing_doy(t_grid, y_grid, p / 100.0, rising=True) for p in PERCENTILES}
            eos = {p: crossing_doy(t_grid, y_grid, p / 100.0, rising=False) for p in PERCENTILES}
            row.update({
                'method': 'double_logistic', 'r2': r2, 'rmse': rmse, 'corr': corr,
                'vmin': popt[0], 'vmax': popt[1], 'S': popt[2], 'greenup_kinetic_i': popt[3],
                'A': popt[4], 'senescence_kinetic_i': popt[5],
                'leaf_out_10': leafout[10], 'leaf_out_50': leafout[50], 'leaf_out_90': leafout[90],
                'EOS10': eos[10], 'EOS50': eos[50], 'EOS90': eos[90],
            })
        else:
            row.update({'method': 'fit_failed', 'r2': np.nan, 'rmse': np.nan, 'corr': np.nan,
                        'vmin': np.nan, 'vmax': np.nan, 'S': np.nan, 'greenup_kinetic_i': np.nan,
                        'A': np.nan, 'senescence_kinetic_i': np.nan,
                        'leaf_out_10': np.nan, 'leaf_out_50': np.nan, 'leaf_out_90': np.nan,
                        'EOS10': np.nan, 'EOS50': np.nan, 'EOS90': np.nan})
        results.append(row)

    if gi % 25 == 0 or gi == total:
        print(f"  processed {gi}/{total} site-years")

results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_CSV, index=False)
n_fit = (results_df['method'] == 'double_logistic').sum()
print(f"\nDone. {n_fit}/{len(results_df)} rows fit successfully -> '{OUTPUT_CSV}'.")
