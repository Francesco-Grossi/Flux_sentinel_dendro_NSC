"""
PIPELINE STEP 8 - Formal test of the core hypothesis:

  CLAIM 1 - GPP has an OPPOSITE effect on EOS90 before vs. after the summer
            solstice (gpp_sos10_to_solstice vs. gpp_solstice_to_eos90).
  CLAIM 2 - Splitting GPP into these two windows IMPROVES prediction of
            EOS90 beyond environmental variables alone (temperature/
            photoperiod senescence rate, radiation, mean temperature,
            respiration) or a single total_gpp_growing_season term.

Claim 1 tested two ways on a mixed-effects model with a site random
intercept (site_id):
  (a) Wald contrasts: is beta_post significantly different from beta_pre?
      Is beta_pre + beta_post consistent with 0 (exact mirror image)?
  (b) Reparameterization: gpp_total (=pre+post) vs. gpp_asymmetry (=post-pre).
      If asymmetry is significant net of total, the TIMING of GPP around
      the solstice matters, not just its sum.

Claim 2 tested two ways:
  (a) Nested-model comparison (ML fit): baseline (no GPP) vs. +total_gpp
      vs. +split GPP, via AIC/BIC and a likelihood-ratio test.
  (b) Leave-one-site-out cross-validation: refit on all-but-one site,
      predict the held-out site with FIXED EFFECTS ONLY (no access to that
      site's own random intercept - mimics predicting a genuinely new
      site), compare RMSE/R^2 across models.

Inputs: data/phenology_flux_predictors_by_site_year_index.csv  (step 6)
        data/rate_of_change_predictors_by_site_year_index.csv  (step 7)

Outputs: data/eos90_opposite_effect_contrast_tests.csv
         data/eos90_gpp_model_comparison.csv
         data/eos90_gpp_cv_results.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
PREDICTORS_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # step 6
RATE_CSV = DATA_DIR / "rate_of_change_predictors_by_site_year_index.csv"        # step 7

OUTPUT_CONTRAST_CSV = DATA_DIR / "eos90_opposite_effect_contrast_tests.csv"
OUTPUT_COMPARISON_CSV = DATA_DIR / "eos90_gpp_model_comparison.csv"
OUTPUT_CV_CSV = DATA_DIR / "eos90_gpp_cv_results.csv"

VI_INDICES = ['NDVI', 'NIRv']
MIN_OBS, MIN_SITES = 20, 5
BASE_PREDICTORS = ['leaf_out_90', 'temperature_senescence_rate', 'photoperiod_senescence_rate']
TARGET = 'EOS90'

for p in (PREDICTORS_CSV, RATE_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run steps 1-7 first.")


def zscore(s):
    return (s - s.mean()) / s.std()


def _scalar(x):
    return float(np.asarray(x).reshape(-1)[0])


pred = pd.read_csv(PREDICTORS_CSV)
rate = pd.read_csv(RATE_CSV)[['site_id', 'year', 'vi_index',
                               'temperature_senescence_rate', 'photoperiod_senescence_rate']]
merged = pred.merge(rate, on=['site_id', 'year', 'vi_index'], how='inner')
print(f"{len(merged)} site-year-index rows available.\n")

contrast_rows, comparison_rows, cv_rows = [], [], []

for vi in VI_INDICES:
    print(f"\n{'=' * 70}\n{vi}\n{'=' * 70}")
    needed = [TARGET] + BASE_PREDICTORS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90', 'site_id']
    d = merged[merged['vi_index'] == vi][needed].dropna().reset_index(drop=True)
    n_obs, n_sites = len(d), d['site_id'].nunique()
    if n_obs < MIN_OBS or n_sites < MIN_SITES:
        print(f"Skipping {vi}: only {n_obs} obs / {n_sites} sites.")
        continue
    print(f"n_obs={n_obs}, n_sites={n_sites}")

    for p in BASE_PREDICTORS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90']:
        d[p + '_z'] = zscore(d[p])

    # --- Claim 1a: Wald contrasts on the split-window model ---
    formula_split = (f"{TARGET} ~ " + " + ".join(p + '_z' for p in BASE_PREDICTORS)
                      + " + gpp_sos10_to_solstice_z + gpp_solstice_to_eos90_z")
    fit_split = smf.mixedlm(formula_split, d, groups=d['site_id']).fit(reml=True)

    fe_names = [n for n in fit_split.params.index if n != 'Group Var']
    idx_pre = fe_names.index('gpp_sos10_to_solstice_z')
    idx_post = fe_names.index('gpp_solstice_to_eos90_z')
    contrast_diff = np.zeros((1, len(fe_names)))
    contrast_diff[0, idx_post], contrast_diff[0, idx_pre] = 1.0, -1.0
    contrast_sum = np.zeros((1, len(fe_names)))
    contrast_sum[0, idx_post], contrast_sum[0, idx_pre] = 1.0, 1.0

    t_diff = fit_split.t_test(contrast_diff)
    t_sum = fit_split.t_test(contrast_sum)

    contrast_rows.append({'vi_index': vi, 'test': 'beta_post - beta_pre = 0 (effects differ)',
                          'estimate': _scalar(t_diff.effect), 'std_err': _scalar(t_diff.sd),
                          't_value': _scalar(t_diff.tvalue), 'p_value': _scalar(t_diff.pvalue)})
    contrast_rows.append({'vi_index': vi, 'test': 'beta_post + beta_pre = 0 (exact mirror image)',
                          'estimate': _scalar(t_sum.effect), 'std_err': _scalar(t_sum.sd),
                          't_value': _scalar(t_sum.tvalue), 'p_value': _scalar(t_sum.pvalue)})
    print(f"\n[Contrast] effects differ: diff={_scalar(t_diff.effect):.2f}, p={_scalar(t_diff.pvalue):.4g}")
    print(f"[Contrast] mirror image:   sum={_scalar(t_sum.effect):.2f}, p={_scalar(t_sum.pvalue):.4g}")

    # --- Claim 1b: total vs. asymmetry reparameterization ---
    d['gpp_total_raw'] = d['gpp_sos10_to_solstice'] + d['gpp_solstice_to_eos90']
    d['gpp_asymmetry_raw'] = d['gpp_solstice_to_eos90'] - d['gpp_sos10_to_solstice']
    d['gpp_total_z'] = zscore(d['gpp_total_raw'])
    d['gpp_asymmetry_z'] = zscore(d['gpp_asymmetry_raw'])

    formula_asym = (f"{TARGET} ~ " + " + ".join(p + '_z' for p in BASE_PREDICTORS)
                     + " + gpp_total_z + gpp_asymmetry_z")
    fit_asym = smf.mixedlm(formula_asym, d, groups=d['site_id']).fit(reml=True)
    contrast_rows.append({'vi_index': vi, 'test': 'gpp_asymmetry_z (timing, net of total)',
                          'estimate': float(fit_asym.params['gpp_asymmetry_z']),
                          'std_err': float(fit_asym.bse['gpp_asymmetry_z']),
                          't_value': float(fit_asym.tvalues['gpp_asymmetry_z']),
                          'p_value': float(fit_asym.pvalues['gpp_asymmetry_z'])})
    contrast_rows.append({'vi_index': vi, 'test': 'gpp_total_z (sum, net of asymmetry)',
                          'estimate': float(fit_asym.params['gpp_total_z']),
                          'std_err': float(fit_asym.bse['gpp_total_z']),
                          't_value': float(fit_asym.tvalues['gpp_total_z']),
                          'p_value': float(fit_asym.pvalues['gpp_total_z'])})
    print(f"[Reparam] asymmetry (timing): coef={fit_asym.params['gpp_asymmetry_z']:.2f}, "
          f"p={fit_asym.pvalues['gpp_asymmetry_z']:.4g}")
    print(f"[Reparam] total (sum):        coef={fit_asym.params['gpp_total_z']:.2f}, "
          f"p={fit_asym.pvalues['gpp_total_z']:.4g}")

    # --- Claim 2a: nested model comparison (ML) ---
    d['total_gpp_growing_season_z'] = zscore(d['gpp_total_raw'])
    f_base = f"{TARGET} ~ " + " + ".join(p + '_z' for p in BASE_PREDICTORS)
    f_total = f_base + " + total_gpp_growing_season_z"
    f_split = f_base + " + gpp_sos10_to_solstice_z + gpp_solstice_to_eos90_z"

    fit_base_ml = smf.mixedlm(f_base, d, groups=d['site_id']).fit(reml=False)
    fit_total_ml = smf.mixedlm(f_total, d, groups=d['site_id']).fit(reml=False)
    fit_split_ml = smf.mixedlm(f_split, d, groups=d['site_id']).fit(reml=False)

    def lrt(fit_restricted, fit_full, df):
        stat = max(2 * (fit_full.llf - fit_restricted.llf), 0.0)
        return stat, chi2.sf(stat, df)

    stat_bt, p_bt = lrt(fit_base_ml, fit_total_ml, 1)
    stat_bs, p_bs = lrt(fit_base_ml, fit_split_ml, 2)
    stat_ts, p_ts = lrt(fit_total_ml, fit_split_ml, 1)

    for name, fit in [('baseline_no_gpp', fit_base_ml), ('total_gpp', fit_total_ml),
                       ('split_gpp_pre_post', fit_split_ml)]:
        comparison_rows.append({'vi_index': vi, 'model': name, 'n_params': int(fit.df_modelwc),
                                'loglik': fit.llf, 'aic': fit.aic, 'bic': fit.bic})
    print(f"\n[AIC] baseline={fit_base_ml.aic:.1f}  total_gpp={fit_total_ml.aic:.1f}  split_gpp={fit_split_ml.aic:.1f}")
    print(f"[LRT] baseline->total_gpp:  chi2={stat_bt:.2f}, p={p_bt:.4g}")
    print(f"[LRT] baseline->split_gpp:  chi2={stat_bs:.2f}, p={p_bs:.4g}")
    print(f"[LRT] total_gpp->split_gpp: chi2={stat_ts:.2f}, p={p_ts:.4g}")

    # --- Claim 2b: leave-one-site-out CV ---
    models_for_cv = {
        'baseline_no_gpp': BASE_PREDICTORS + [],
        'total_gpp': BASE_PREDICTORS + ['gpp_total_raw'],
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
    print("\n[Leave-one-site-out CV]")
    for name, yhat in preds.items():
        mask = ~np.isnan(yhat)
        resid = y_true[mask] - yhat[mask]
        rmse = float(np.sqrt(np.mean(resid ** 2)))
        ss_res, ss_tot = float(np.sum(resid ** 2)), float(np.sum((y_true[mask] - y_true[mask].mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        cv_rows.append({'vi_index': vi, 'model': name, 'n_predicted': int(mask.sum()),
                        'loso_rmse_days': rmse, 'loso_r2': r2})
        print(f"  {name:22s} n={int(mask.sum()):4d}  RMSE={rmse:6.2f}d  R2={r2:.3f}")

pd.DataFrame(contrast_rows).to_csv(OUTPUT_CONTRAST_CSV, index=False)
pd.DataFrame(comparison_rows).to_csv(OUTPUT_COMPARISON_CSV, index=False)
pd.DataFrame(cv_rows).to_csv(OUTPUT_CV_CSV, index=False)
print(f"\n\nWritten: '{OUTPUT_CONTRAST_CSV.name}', '{OUTPUT_COMPARISON_CSV.name}', '{OUTPUT_CV_CSV.name}'")
