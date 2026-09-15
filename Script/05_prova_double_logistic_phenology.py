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

MIN_POINTS_FOR_FIT = 15   # double-logistic has 6 free params; 12 was too thin a
                           # margin and let genuinely underdetermined fits through
MIN_VALID_FRAC = 0.5      # QC floor on Fmask-derived clear-pixel fraction

# Upper bound on the kinetic (steepness) rate parameters during fitting.
# The 10-90% transition of a logistic phase takes roughly 4.4/i days, so
# i=5.0 (the old bound) allows a full green-up or senescence transition to
# complete in under a day - a near-vertical step, not a real phenological
# transition. That's exactly the failure mode seen in fits that report a
# deceptively high r/low RMSE (a handful of points can look "well fit" by
# a step) while being visually nonsensical. Capping at 1.0 still allows a
# fast ~4-5 day transition, comfortably covering real abrupt events (e.g.
# frost-triggered senescence), while ruling out step-function degenerate
# solutions.
MAX_KINETIC_RATE = 1.0

# How far outside the actual observed DOY range a fitted inflection point
# (S or A) is allowed to fall before the fit is rejected as implausible -
# an inflection point extrapolated far past the data isn't meaningfully
# constrained by any observation.
MAX_INFLECTION_EXTRAPOLATION_DAYS = 30

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
        max(peak_guess - 20, trough_guess_lo - MAX_INFLECTION_EXTRAPOLATION_DAYS),   # S: green-up inflection guess
        rate_guess,                                    # i_sos > 0
        min(peak_guess + 20, trough_guess_hi + MAX_INFLECTION_EXTRAPOLATION_DAYS),   # A: senescence inflection guess
        rate_guess,                                     # i_eos > 0 (sign handled inside the model)
    ]
    lower = [vmin_obs - amp, vmax_obs - amp, doy.min() - MAX_INFLECTION_EXTRAPOLATION_DAYS, 1e-4,
             doy.min() - MAX_INFLECTION_EXTRAPOLATION_DAYS, 1e-4]
    upper = [vmin_obs + amp, vmax_obs + amp, doy.max() + MAX_INFLECTION_EXTRAPOLATION_DAYS, MAX_KINETIC_RATE,
             doy.max() + MAX_INFLECTION_EXTRAPOLATION_DAYS, MAX_KINETIC_RATE]
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


def fit_is_plausible(popt, doy_min, doy_max):
    """Post-fit sanity check, independent of (and complementary to) the
    numeric r2/rmse/corr metrics. A fit can score deceptively well on those
    metrics - a handful of scattered points can make a near-vertical step
    look "well correlated" - while still being a biologically nonsensical
    curve (e.g. senescence starting before green-up, or a transition so
    abrupt it's pinned at the fitting bound). Returns (True, None) if the
    fit passes, or (False, reason) if it should be rejected outright."""
    vmin, vmax, S, i_sos, A, i_eos = popt

    if not (S < A):
        return False, "S >= A (senescence inflection at or before green-up inflection)"

    lo = doy_min - MAX_INFLECTION_EXTRAPOLATION_DAYS
    hi = doy_max + MAX_INFLECTION_EXTRAPOLATION_DAYS
    if not (lo <= S <= hi):
        return False, "S extrapolated far outside the observed DOY range"
    if not (lo <= A <= hi):
        return False, "A extrapolated far outside the observed DOY range"

    # A kinetic rate pinned at (or essentially at) the upper bound means the
    # optimizer wanted an even steeper/more step-like transition than
    # allowed - a sign the data doesn't actually constrain a smooth curve,
    # not a confirmed "fast but real" transition.
    pin_tolerance = 1e-3
    if i_sos >= MAX_KINETIC_RATE - pin_tolerance:
        return False, "green-up kinetic rate pinned at the upper bound (near-step transition)"
    if i_eos >= MAX_KINETIC_RATE - pin_tolerance:
        return False, "senescence kinetic rate pinned at the upper bound (near-step transition)"

    return True, None


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

        # Even a numerically successful fit can be an implausible curve
        # (see fit_is_plausible) - check that BEFORE treating popt as
        # usable, so a degenerate fit is reported the same way as a failed
        # one (NaN metrics, excluded downstream) rather than contributing
        # spurious percentile crossings.
        rejection_reason = None
        if popt is not None:
            plausible, rejection_reason = fit_is_plausible(popt, doy_all.min(), doy_all.max())
            if not plausible:
                popt = None

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

            # RMSE alone is scale-dependent (0.05 means something different
            # for an index ranging 0.6-0.9 vs. 0.1-0.9) - normalizing by the
            # OBSERVED amplitude gives a comparable "fraction of the range"
            # error, useful as an additional filter alongside corr.
            observed_amp = max(float(np.nanmax(vi_all) - np.nanmin(vi_all)), 1e-6)
            relative_rmse = rmse / observed_amp

            row.update({
                'method': 'double_logistic',
                'r2': r2, 'rmse': rmse, 'relative_rmse': relative_rmse, 'corr': corr,
                'vmin': popt[0], 'vmax': popt[1],
                'S': popt[2], 'greenup_kinetic_i': popt[3],
                'A': popt[4], 'senescence_kinetic_i': popt[5],
                'greenup_corr': gu_corr, 'greenup_rmse': gu_rmse,
                'leaf_out_10': leafout[10], 'leaf_out_50': leafout[50], 'leaf_out_90': leafout[90],
                'senescence_corr': se_corr, 'senescence_rmse': se_rmse,
                'EOS10': eos[10], 'EOS50': eos[50], 'EOS90': eos[90],
            })
        else:
            # No usable double-logistic curve - either the optimizer never
            # converged/had too few points (fit_double_logistic returned
            # None), or it converged to a numerically-fine but implausible
            # curve (fit_is_plausible rejected it, rejection_reason set).
            # Either way: report NaN across the board rather than
            # estimating percentile crossings some other way. A crossing
            # DOY derived without a genuinely usable fitted curve isn't
            # comparable to the ones that are, so it's not reported at all.
            row.update({
                'method': 'fit_rejected_implausible' if rejection_reason else 'fit_failed',
                'rejection_reason': rejection_reason,
                'r2': np.nan, 'rmse': np.nan, 'relative_rmse': np.nan, 'corr': np.nan,
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
n_implausible = (results_df['method'] == 'fit_rejected_implausible').sum()
print(f"\nDone. {len(results_df)} site-year-index rows written to '{OUTPUT_CSV}'.")
print(f"Fits: {n_fit} double-logistic succeeded, {n_failed} failed "
      f"(too few points or non-convergent), {n_implausible} rejected as implausible "
      "(converged, but to a biologically nonsensical curve - see 'rejection_reason'). "
      "All metrics NaN for failed/rejected rows.")
if n_implausible:
    print("\nImplausible-fit rejection reasons:")
    print(results_df.loc[results_df['method'] == 'fit_rejected_implausible', 'rejection_reason']
          .value_counts().to_string())