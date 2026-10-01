"""
PIPELINE STEP 12 - Same hypothesis test as step 08 (does GPP have an
opposite effect on EOS90 before vs. after the summer solstice, and does
splitting it improve prediction?), but run SEPARATELY for deciduous vs.
evergreen sites, to check whether the effect is a general phenomenon or
specific to one leaf habit.

Deciduous vs. evergreen is read off each site's IGBP land-cover class
(carried through from step 1's raw download into the step 4 merged table):
    deciduous: DBF (Deciduous Broadleaf Forest), DNF (Deciduous Needleleaf Forest)
    evergreen: ENF (Evergreen Needleleaf Forest), EBF (Evergreen Broadleaf Forest)
Everything else (grassland, savanna, shrubland, mixed forest, ...) doesn't
map cleanly onto a binary deciduous/evergreen split and is EXCLUDED from
this script rather than force-labeled - if you need those included too,
they should probably be a third group, not folded into either side.

This reuses the exact same Wald-contrast / reparameterization / LRT / LOSO-CV
logic as step 08 - see that script's docstring for the full explanation of
each test. Here it's just looped over (leaf_type x vi_index) instead of
just (vi_index).

Inputs: data/fluxnet_landsat_merged.csv                        (step 4, for IGBP per site)
        data/phenology_flux_predictors_by_site_year_index.csv  (step 6)
        data/rate_of_change_predictors_by_site_year_index.csv  (step 7)

Outputs: data/eos90_deciduous_vs_evergreen_contrast_tests.csv
         data/eos90_deciduous_vs_evergreen_model_comparison.csv
         data/eos90_deciduous_vs_evergreen_cv.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                              # step 4
PREDICTORS_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # step 6
RATE_CSV = DATA_DIR / "rate_of_change_predictors_by_site_year_index.csv"        # step 7

OUTPUT_CONTRAST_CSV = DATA_DIR / "eos90_deciduous_vs_evergreen_contrast_tests.csv"
OUTPUT_COMPARISON_CSV = DATA_DIR / "eos90_deciduous_vs_evergreen_model_comparison.csv"
OUTPUT_CV_CSV = DATA_DIR / "eos90_deciduous_vs_evergreen_cv.csv"

VI_INDICES = ['NDVI', 'NIRv']
MIN_OBS, MIN_SITES = 15, 4  # relaxed vs. step 08's 20/5 - splitting by leaf type roughly
                             # halves (or worse) the sample, so the pooled-analysis
                             # thresholds are usually too strict to fit anything at all
BASE_PREDICTORS = ['leaf_out_90', 'temperature_senescence_rate', 'photoperiod_senescence_rate']
TARGET = 'EOS90'

DECIDUOUS_IGBP = {'DBF', 'DNF'}
EVERGREEN_IGBP = {'ENF', 'EBF'}

for p in (FLUX_CSV, PREDICTORS_CSV, RATE_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run steps 1-7 first.")


def classify_land_cover(igbp_raw):
    """Same alias mapping used in step 1, kept local here so this script
    doesn't depend on importing step 1 as a module."""
    if igbp_raw is None or pd.isna(igbp_raw):
        return None
    s = str(igbp_raw).strip().upper()
    aliases = {
        'DECIDUOUS BROADLEAF FOREST': 'DBF', 'DECIDUOUS BROADLEAF FORESTS': 'DBF',
        'DECIDUOUS NEEDLELEAF FOREST': 'DNF', 'DECIDUOUS NEEDLELEAF FORESTS': 'DNF',
        'EVERGREEN NEEDLELEAF FOREST': 'ENF', 'EVERGREEN NEEDLELEAF FORESTS': 'ENF',
        'EVERGREEN BROADLEAF FOREST': 'EBF', 'EVERGREEN BROADLEAF FORESTS': 'EBF',
        'MIXED FOREST': 'MF', 'MIXED FORESTS': 'MF',
        'GRASSLAND': 'GRA', 'GRASSLANDS': 'GRA',
    }
    return aliases.get(s, s[:3])


def leaf_habit(igbp_code):
    if igbp_code in DECIDUOUS_IGBP:
        return 'deciduous'
    if igbp_code in EVERGREEN_IGBP:
        return 'evergreen'
    return None


def zscore(s):
    return (s - s.mean()) / s.std()


def _scalar(x):
    return float(np.asarray(x).reshape(-1)[0])


# ---------------------------------------------------------------------------
# 1. Classify each site's leaf habit from its IGBP class
# ---------------------------------------------------------------------------
flux_igbp = pd.read_csv(FLUX_CSV, usecols=['site_id', 'igbp'])
site_igbp = flux_igbp.groupby('site_id')['igbp'].first()
site_habit = site_igbp.apply(classify_land_cover).apply(leaf_habit)
site_habit = site_habit.dropna()
print(f"Leaf-habit classification: {int((site_habit == 'deciduous').sum())} deciduous sites, "
      f"{int((site_habit == 'evergreen').sum())} evergreen sites "
      f"({site_igbp.shape[0] - len(site_habit)} sites excluded - land cover is neither).")

# ---------------------------------------------------------------------------
# 2. Merge predictors + rates + leaf habit
# ---------------------------------------------------------------------------
pred = pd.read_csv(PREDICTORS_CSV)
rate = pd.read_csv(RATE_CSV)[['site_id', 'year', 'vi_index',
                               'temperature_senescence_rate', 'photoperiod_senescence_rate']]
merged = pred.merge(rate, on=['site_id', 'year', 'vi_index'], how='inner')
merged['leaf_habit'] = merged['site_id'].map(site_habit)
merged = merged.dropna(subset=['leaf_habit'])
print(f"{len(merged)} site-year-index rows with a deciduous/evergreen label.\n")

contrast_rows, comparison_rows, cv_rows = [], [], []

for leaf_type in ['deciduous', 'evergreen']:
    for vi in VI_INDICES:
        print(f"\n{'=' * 70}\n{leaf_type.upper()} | {vi}\n{'=' * 70}")
        needed = [TARGET] + BASE_PREDICTORS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90', 'site_id']
        d = merged[(merged['leaf_habit'] == leaf_type) & (merged['vi_index'] == vi)][needed] \
            .dropna().reset_index(drop=True)
        n_obs, n_sites = len(d), d['site_id'].nunique()
        if n_obs < MIN_OBS or n_sites < MIN_SITES:
            print(f"Skipping {leaf_type}/{vi}: only {n_obs} obs / {n_sites} sites "
                  f"(need >= {MIN_OBS} obs, >= {MIN_SITES} sites).")
            continue
        print(f"n_obs={n_obs}, n_sites={n_sites}")

        for p in BASE_PREDICTORS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90']:
            d[p + '_z'] = zscore(d[p])

        # --- Claim 1a: Wald contrasts on the split-window model ---
        formula_split = (f"{TARGET} ~ " + " + ".join(p + '_z' for p in BASE_PREDICTORS)
                          + " + gpp_sos10_to_solstice_z + gpp_solstice_to_eos90_z")
        try:
            fit_split = smf.mixedlm(formula_split, d, groups=d['site_id']).fit(reml=True)
        except Exception as e:
            print(f"  Model failed to converge ({e}); skipping {leaf_type}/{vi}.")
            continue

        fe_names = [n for n in fit_split.params.index if n != 'Group Var']
        idx_pre = fe_names.index('gpp_sos10_to_solstice_z')
        idx_post = fe_names.index('gpp_solstice_to_eos90_z')
        contrast_diff = np.zeros((1, len(fe_names)))
        contrast_diff[0, idx_post], contrast_diff[0, idx_pre] = 1.0, -1.0
        contrast_sum = np.zeros((1, len(fe_names)))
        contrast_sum[0, idx_post], contrast_sum[0, idx_pre] = 1.0, 1.0

        t_diff = fit_split.t_test(contrast_diff)
        t_sum = fit_split.t_test(contrast_sum)

        contrast_rows.append({'leaf_habit': leaf_type, 'vi_index': vi,
                              'test': 'beta_post - beta_pre = 0 (effects differ)',
                              'estimate': _scalar(t_diff.effect), 'std_err': _scalar(t_diff.sd),
                              't_value': _scalar(t_diff.tvalue), 'p_value': _scalar(t_diff.pvalue),
                              'n_obs': n_obs, 'n_sites': n_sites})
        contrast_rows.append({'leaf_habit': leaf_type, 'vi_index': vi,
                              'test': 'beta_post + beta_pre = 0 (exact mirror image)',
                              'estimate': _scalar(t_sum.effect), 'std_err': _scalar(t_sum.sd),
                              't_value': _scalar(t_sum.tvalue), 'p_value': _scalar(t_sum.pvalue),
                              'n_obs': n_obs, 'n_sites': n_sites})
        print(f"[Contrast] effects differ: diff={_scalar(t_diff.effect):.2f}, p={_scalar(t_diff.pvalue):.4g}")

        # --- Claim 1b: total vs. asymmetry reparameterization ---
        d['gpp_total_raw'] = d['gpp_sos10_to_solstice'] + d['gpp_solstice_to_eos90']
        d['gpp_asymmetry_raw'] = d['gpp_solstice_to_eos90'] - d['gpp_sos10_to_solstice']
        d['gpp_total_z'] = zscore(d['gpp_total_raw'])
        d['gpp_asymmetry_z'] = zscore(d['gpp_asymmetry_raw'])

        formula_asym = (f"{TARGET} ~ " + " + ".join(p + '_z' for p in BASE_PREDICTORS)
                         + " + gpp_total_z + gpp_asymmetry_z")
        fit_asym = smf.mixedlm(formula_asym, d, groups=d['site_id']).fit(reml=True)
        contrast_rows.append({'leaf_habit': leaf_type, 'vi_index': vi,
                              'test': 'gpp_asymmetry_z (timing, net of total)',
                              'estimate': float(fit_asym.params['gpp_asymmetry_z']),
                              'std_err': float(fit_asym.bse['gpp_asymmetry_z']),
                              't_value': float(fit_asym.tvalues['gpp_asymmetry_z']),
                              'p_value': float(fit_asym.pvalues['gpp_asymmetry_z']),
                              'n_obs': n_obs, 'n_sites': n_sites})
        print(f"[Reparam] asymmetry (timing): coef={fit_asym.params['gpp_asymmetry_z']:.2f}, "
              f"p={fit_asym.pvalues['gpp_asymmetry_z']:.4g}")

        # --- Claim 2a: nested model comparison (ML) ---
        d['total_gpp_growing_season_z'] = zscore(d['gpp_total_raw'])
        f_base = f"{TARGET} ~ " + " + ".join(p + '_z' for p in BASE_PREDICTORS)
        f_split = f_base + " + gpp_sos10_to_solstice_z + gpp_solstice_to_eos90_z"

        fit_base_ml = smf.mixedlm(f_base, d, groups=d['site_id']).fit(reml=False)
        fit_split_ml = smf.mixedlm(f_split, d, groups=d['site_id']).fit(reml=False)
        stat = max(2 * (fit_split_ml.llf - fit_base_ml.llf), 0.0)
        p_lrt = chi2.sf(stat, 2)

        for name, fit in [('baseline_no_gpp', fit_base_ml), ('split_gpp_pre_post', fit_split_ml)]:
            comparison_rows.append({'leaf_habit': leaf_type, 'vi_index': vi, 'model': name,
                                    'n_params': int(fit.df_modelwc), 'loglik': fit.llf,
                                    'aic': fit.aic, 'bic': fit.bic, 'n_obs': n_obs, 'n_sites': n_sites})
        print(f"[LRT] baseline -> split_gpp: chi2={stat:.2f}, df=2, p={p_lrt:.4g}")

        # --- Claim 2b: leave-one-site-out CV ---
        models_for_cv = {
            'baseline_no_gpp': BASE_PREDICTORS + [],
            'split_gpp_pre_post': BASE_PREDICTORS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90'],
        }
        sites = d['site_id'].unique()
        preds = {name: np.full(len(d), np.nan) for name in models_for_cv}
        for held_out_site in sites:
            train_mask = d['site_id'] != held_out_site
            test_mask = ~train_mask
            if test_mask.sum() == 0 or train_mask.sum() < MIN_OBS:
                continue
            for name, raw_preds in models_for_cv.items():
                train, test = d.loc[train_mask].copy(), d.loc[test_mask].copy()
                means = {p: train[p].mean() for p in raw_preds}
                stds = {p: train[p].std() for p in raw_preds}
                for p in raw_preds:
                    train[p + '_zc'] = (train[p] - means[p]) / stds[p]
                    test[p + '_zc'] = (test[p] - means[p]) / stds[p]
                formula_cv = f"{TARGET} ~ " + " + ".join(p + '_zc' for p in raw_preds)
                try:
                    fit_cv = smf.mixedlm(formula_cv, train, groups=train['site_id']).fit(reml=True)
                except Exception:
                    continue
                coef_names = ['Intercept'] + [p + '_zc' for p in raw_preds]
                X_test = np.column_stack([np.ones(len(test))] + [test[p + '_zc'].to_numpy() for p in raw_preds])
                coefs = np.array([fit_cv.params.get(c, 0.0) for c in coef_names])
                preds[name][test.index.to_numpy()] = X_test @ coefs

        y_true = d[TARGET].to_numpy()
        print("[Leave-one-site-out CV]")
        for name, yhat in preds.items():
            mask = ~np.isnan(yhat)
            if mask.sum() == 0:
                continue
            resid = y_true[mask] - yhat[mask]
            rmse = float(np.sqrt(np.mean(resid ** 2)))
            ss_res, ss_tot = float(np.sum(resid ** 2)), float(np.sum((y_true[mask] - y_true[mask].mean()) ** 2))
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            cv_rows.append({'leaf_habit': leaf_type, 'vi_index': vi, 'model': name,
                            'n_predicted': int(mask.sum()), 'loso_rmse_days': rmse, 'loso_r2': r2})
            print(f"  {name:22s} n={int(mask.sum()):4d}  RMSE={rmse:6.2f}d  R2={r2:.3f}")

pd.DataFrame(contrast_rows).to_csv(OUTPUT_CONTRAST_CSV, index=False)
pd.DataFrame(comparison_rows).to_csv(OUTPUT_COMPARISON_CSV, index=False)
pd.DataFrame(cv_rows).to_csv(OUTPUT_CV_CSV, index=False)
print(f"\n\nWritten: '{OUTPUT_CONTRAST_CSV.name}', '{OUTPUT_COMPARISON_CSV.name}', '{OUTPUT_CV_CSV.name}'")

if not contrast_rows:
    print("\nWARNING: no leaf_type/vi_index combination had enough data to fit. "
          "This usually means too few sites survive both the deciduous/evergreen "
          "classification AND the stringent phenology-fit quality gate together - "
          "consider relaxing MIN_OBS/MIN_SITES above, or checking how many sites "
          "actually have a DBF/DNF/ENF/EBF land cover in your data.")
