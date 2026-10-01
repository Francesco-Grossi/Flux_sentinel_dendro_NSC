"""
Shared double-logistic phenology fitting, used by every step that turns a
greenness time series into leaf-out / EOS dates:
    step 05  satellite VIs (HLS NDVI / NIRv)
    step 19  tower broadband NDVI
    step 20  PhenoCam GCC
NOT a pipeline step itself. Using one fitting routine for all sources keeps
EOS90 / EOS50 / EOS10 defined identically, so they can be compared directly
(step 21) and swapped as the outcome in steps 12-17.

Sparse or noisy series give unstable autumn fits. Four safeguards, in order:

1. Dormant-season background (Beck et al. 2006). Per site x index, the winter
   background is the median value of clear observations on frozen days
   (15-day running mean TA_F < 0 C). Values below it are raised to it
   (residual snow/cloud dips), and the background is inserted as
   pseudo-observations on (a) snow-covered dates ('snowy') and (b) frozen
   days with no clear observation within +/-WINTER_FILL_STEP days.
   Sites without frozen periods get no winter fill.
2. Climatology gap-fill (as in MSLSP30NA, Bolton et al. 2020). A multi-year
   curve is fit per site x index on all years pooled by DOY. In each year,
   gaps > GAP_DAYS between observations are filled with that curve (shifted
   by the year's median offset from it) at low weight. Only gaps are filled;
   where a year has data, its own observations decide the dates.
3. Transition dates are read off the fitted curve over the full year
   (DOY 1-366), relative to the fitted amplitude (peak - vmin).
4. QC flags (qc_pass / qc_reason). Failing site-years keep their fit
   parameters but get NaN phenology dates, so downstream steps drop them:
     - >= MIN_OBS_AFTER_DOY200 real observations after DOY 200
     - fit R2 >= MIN_R2 on the real (not pseudo) observations
     - amplitude >= MIN_AMP_TO_RMSE x fit RMSE (seasonal signal above noise)
     - leaf_out_10 < 50 < 90 < peak < EOS90 < EOS50 < EOS10
     - EOS90 >= MIN_EOS90_DOY
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

PERCENTILES = [10, 50, 90]
MIN_REAL_OBS = 12           # real (non-pseudo) clear observations per site-year; 8 let noisy years in

# 1. dormant-season background
FROZEN_TA_C = 0.0
TA_WINDOW_DAYS = 15
MIN_BACKGROUND_OBS = 5
WINTER_FILL_STEP = 8
W_SNOW, W_WINTER = 0.5, 0.3

# 2. climatology gap-fill
GAP_DAYS = 30
CLIM_FILL_STEP = 10
W_CLIM = 0.2

# 4. QC
MIN_OBS_AFTER_DOY200 = 3
MIN_R2 = 0.8
MIN_AMP_TO_RMSE = 2.0
MIN_EOS90_DOY = 180

T_GRID = np.arange(1.0, 366.01, 0.25)

OUTPUT_COLUMNS = ['site_id', 'year', 'vi_index', 'n_obs', 'n_obs_after_doy200', 'n_winter_fill', 'n_clim_fill',
                  'background', 'method', 'r2', 'rmse', 'corr', 'vmin', 'vmax', 'S', 'greenup_kinetic_i', 'A',
                  'senescence_kinetic_i', 'peak_doy', 'amplitude',
                  'leaf_out_10', 'leaf_out_50', 'leaf_out_90', 'EOS10', 'EOS50', 'EOS90', 'qc_pass', 'qc_reason']


def double_logistic(t, vmin, vmax, S, i_sos, A, i_eos):
    z_S = np.clip(-i_sos * (t - S), -500, 500)
    z_A = np.clip(-i_eos * (t - A), -500, 500)
    return vmin + (vmax - vmin) * (1.0 / (1.0 + np.exp(z_S)) - 1.0 / (1.0 + np.exp(z_A)))


def _residuals(params, t, y, w):
    return w * (double_logistic(t, *params) - y)


def fit_double_logistic(doy, vi, weights):
    """Robust weighted fit; returns params or None. S (green-up) is kept in
    the first part of the year and A (senescence) after it."""
    vmin_obs, vmax_obs = np.nanpercentile(vi, 2), np.nanpercentile(vi, 98)
    amp = max(vmax_obs - vmin_obs, 1e-6)
    peak_guess = float(np.clip(doy[np.argmax(vi)], 120, 240))
    p0 = [vmin_obs, vmax_obs, peak_guess - 40, 0.1, peak_guess + 60, 0.1]
    lower = [vmin_obs - amp, vmin_obs, 1, 0.01, 120, 0.01]
    upper = [vmax_obs, vmax_obs + amp, 240, 1.0, 366, 1.0]
    p0 = list(np.clip(p0, np.array(lower) + 1e-6, np.array(upper) - 1e-6))
    try:
        result = least_squares(_residuals, p0, args=(doy, vi, weights), bounds=(lower, upper),
                               loss='soft_l1', f_scale=0.1 * amp, max_nfev=5000)
    except (RuntimeError, ValueError):
        return None
    popt = result.x
    return popt if popt[2] < popt[4] else None


def transition_dates(popt):
    """Crossing DOYs relative to the fitted amplitude: frac = (y - vmin)/(peak - vmin)."""
    y = double_logistic(T_GRID, *popt)
    ipk = int(np.argmax(y))
    amp = y[ipk] - popt[0]
    out = {'peak_doy': T_GRID[ipk], 'amplitude': amp}
    if amp <= 0:
        return out
    frac = (y - popt[0]) / amp
    for p in PERCENTILES:
        rise = np.where(frac[:ipk + 1] >= p / 100)[0]
        fall = np.where(frac[ipk:] <= p / 100)[0]
        out[f'leaf_out_{p}'] = T_GRID[rise[0]] if len(rise) and rise[0] > 0 else np.nan
        out[f'EOS{p}'] = T_GRID[ipk + fall[0]] if len(fall) else np.nan
    return out


def fill_points(doy_obs, lo=1, hi=366, gap=GAP_DAYS, step=CLIM_FILL_STEP):
    """DOYs inside gaps > `gap` days between observations (incl. year edges)."""
    edges = np.concatenate([[lo - 1], np.sort(doy_obs), [hi + 1]])
    pts = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a > gap:
            pts.extend(np.arange(a + step / 2, b - step / 2 + 1e-9, step))
    return np.array(pts, dtype=float)


def qc_check(row, n_after_200):
    reasons = []
    if n_after_200 < MIN_OBS_AFTER_DOY200:
        reasons.append('few_autumn_obs')
    if not (row.get('r2', np.nan) >= MIN_R2):
        reasons.append('low_r2')
    if not (row.get('amplitude', 0) >= MIN_AMP_TO_RMSE * row.get('rmse', np.inf)):
        reasons.append('amp_below_noise')
    seq = [row.get(k, np.nan) for k in ['leaf_out_10', 'leaf_out_50', 'leaf_out_90', 'peak_doy',
                                         'EOS90', 'EOS50', 'EOS10']]
    if np.any(np.isnan(seq)) or np.any(np.diff(seq) <= 0):
        reasons.append('dates_missing_or_out_of_order')
    if not (row.get('EOS90', np.nan) >= MIN_EOS90_DOY):
        reasons.append('eos90_too_early')
    return len(reasons) == 0, ';'.join(reasons)


def add_frozen_flag(df, ta_col='TA_F'):
    """Adds 'frozen' (15-day running mean air temperature < 0 C) per site.
    df must be daily, sorted by (site_id, date)."""
    ta_roll = df.groupby('site_id')[ta_col].transform(
        lambda s: s.rolling(TA_WINDOW_DAYS, center=True, min_periods=5).mean())
    df['frozen'] = ta_roll < FROZEN_TA_C
    return df


def fit_site_index(site_id, sdf, value_col, label=None):
    """All site-years of one site x index. sdf: one row per day (or per
    observation date) with columns year, doy, frozen, snowy, fit_weight and
    `value_col` (NaN where there is no clear observation). Returns a list of
    row dicts with OUTPUT_COLUMNS keys; vi_index is `label` or `value_col`."""
    label = label or value_col
    rows = []
    obs = sdf.loc[sdf[value_col].notna() & ~sdf['snowy'], ['year', 'doy', 'frozen', value_col, 'fit_weight']] \
        .rename(columns={value_col: 'vi', 'fit_weight': 'w'})
    if obs.empty:
        return rows

    # 1. dormant background
    winter_obs = obs.loc[obs['frozen'], 'vi']
    background = winter_obs.median() if len(winter_obs) >= MIN_BACKGROUND_OBS else np.nan
    if np.isnan(background) and (sdf['frozen'].any() or sdf['snowy'].any()):
        background = obs['vi'].quantile(0.10)   # frozen/snowy site but few clear winter observations
    if not np.isnan(background):
        obs['vi'] = obs['vi'].clip(lower=background)

    def winter_pseudo(ydf, yobs_doy):
        if np.isnan(background):
            return pd.DataFrame(columns=['doy', 'vi', 'w'])
        snow_doy = ydf.loc[ydf['snowy'], 'doy'].to_numpy(float)
        frozen_doy = ydf.loc[ydf['frozen'], 'doy'].to_numpy(float)
        frozen_doy = frozen_doy[::WINTER_FILL_STEP]
        if len(yobs_doy) and len(frozen_doy):
            near = np.abs(frozen_doy[:, None] - yobs_doy[None, :]).min(axis=1) <= WINTER_FILL_STEP
            frozen_doy = frozen_doy[~near]
        return pd.concat([pd.DataFrame({'doy': snow_doy, 'vi': background, 'w': W_SNOW}),
                          pd.DataFrame({'doy': frozen_doy, 'vi': background, 'w': W_WINTER})])

    # 2. site climatology: all years pooled (real obs + winter pseudo-obs)
    pooled = [obs[['doy', 'vi', 'w']]]
    for year, ydf in sdf.groupby('year'):
        pooled.append(winter_pseudo(ydf, obs.loc[obs['year'] == year, 'doy'].to_numpy(float)))
    pooled = pd.concat([x for x in pooled if len(x)], ignore_index=True).astype(float)
    clim = fit_double_logistic(pooled['doy'].to_numpy(), pooled['vi'].to_numpy(), pooled['w'].to_numpy()) \
        if len(pooled) >= 3 * MIN_REAL_OBS else None

    for year, ydf in sdf.groupby('year'):
        yobs = obs[obs['year'] == year]
        row = {'site_id': site_id, 'year': year, 'vi_index': label, 'n_obs': len(yobs),
               'n_obs_after_doy200': int((yobs['doy'] > 200).sum()), 'background': background}
        if len(yobs) < MIN_REAL_OBS:
            continue
        real = yobs[['doy', 'vi', 'w']].astype(float)
        winter = winter_pseudo(ydf, real['doy'].to_numpy())
        parts = [real, winter]
        if clim is not None:
            offset = np.median(real['vi'] - double_logistic(real['doy'].to_numpy(), *clim))
            have = np.concatenate([real['doy'].to_numpy(), winter['doy'].to_numpy(float)])
            cdoy = fill_points(have)
            parts.append(pd.DataFrame({'doy': cdoy, 'vi': double_logistic(cdoy, *clim) + offset,
                                       'w': W_CLIM}))
        fit_df = pd.concat([x for x in parts if len(x)], ignore_index=True).astype(float).sort_values('doy')
        row.update({'n_winter_fill': len(winter), 'n_clim_fill': len(parts[2]) if clim is not None else 0})

        popt = fit_double_logistic(fit_df['doy'].to_numpy(), fit_df['vi'].to_numpy(), fit_df['w'].to_numpy())
        if popt is None:
            row.update({'method': 'fit_failed', 'qc_pass': False, 'qc_reason': 'fit_failed'})
            rows.append(row)
            continue

        pred = double_logistic(real['doy'].to_numpy(), *popt)
        ss_res = np.sum((real['vi'] - pred) ** 2)
        ss_tot = np.sum((real['vi'] - real['vi'].mean()) ** 2)
        row.update({
            'method': 'double_logistic',
            'r2': 1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
            'rmse': float(np.sqrt(np.mean((real['vi'] - pred) ** 2))),
            'corr': np.corrcoef(real['vi'], pred)[0, 1] if np.std(pred) > 0 else np.nan,
            'vmin': popt[0], 'vmax': popt[1], 'S': popt[2], 'greenup_kinetic_i': popt[3],
            'A': popt[4], 'senescence_kinetic_i': popt[5],
        })
        row.update(transition_dates(popt))
        ok, reason = qc_check(row, row['n_obs_after_doy200'])
        row.update({'qc_pass': ok, 'qc_reason': reason})
        if not ok:
            for p in PERCENTILES:
                row[f'leaf_out_{p}'] = np.nan
                row[f'EOS{p}'] = np.nan
        rows.append(row)
    return rows


def write_results(rows, output_csv):
    results_df = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)
    results_df.to_csv(output_csv, index=False)
    print(f"\nDone -> '{output_csv}'.")
    for vi, d in results_df.groupby('vi_index'):
        reasons = d.loc[~d['qc_pass'].astype(bool), 'qc_reason'].str.split(';').explode().value_counts()
        print(f"  {vi}: {int(d['qc_pass'].sum())}/{len(d)} site-years pass QC "
              f"({d.loc[d['qc_pass'].astype(bool), 'site_id'].nunique()} sites); fail reasons: {reasons.to_dict()}")
    return results_df
