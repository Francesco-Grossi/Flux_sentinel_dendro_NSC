"""
PIPELINE STEP 9 - Four-way model comparison for EOS90:

  A. environment only   : temperature_senescence_rate, photoperiod_senescence_rate,
                           mean_radiation_growing_season, mean_temperature_growing_season,
                           total_reco_growing_season
  B. environment + leaf_out_90
  C. environment + split GPP (gpp_sos10_to_solstice, gpp_solstice_to_eos90)
  D. environment + leaf_out_90 + split GPP

AIC/BIC (ML fit), pairwise likelihood-ratio tests between nested pairs, and
leave-one-site-out CV RMSE/R^2 - the out-of-sample step is what actually
disambiguates "fits better in-sample" from "generalizes to a new site",
and in this pipeline the two VI indices (NDVI, NIRv) turn out to pick up
different predictive signal (leaf_out for NDVI, split GPP for NIRv), which
only the CV step reveals - AIC alone favors adding predictors for both.

Inputs: data/phenology_flux_predictors_by_site_year_index.csv  (step 6)
        data/rate_of_change_predictors_by_site_year_index.csv  (step 7)

Outputs: data/eos90_four_model_comparison.csv
         data/eos90_four_model_lrt.csv
         data/eos90_four_model_cv.csv
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

OUTPUT_COMPARISON_CSV = DATA_DIR / "eos90_four_model_comparison.csv"
OUTPUT_LRT_CSV = DATA_DIR / "eos90_four_model_lrt.csv"
OUTPUT_CV_CSV = DATA_DIR / "eos90_four_model_cv.csv"

VI_INDICES = ['NDVI', 'NIRv']
MIN_OBS, MIN_SITES = 20, 5
TARGET = 'EOS90'

ENV_VARS = ['temperature_senescence_rate', 'photoperiod_senescence_rate',
            'mean_radiation_growing_season', 'mean_temperature_growing_season',
            'total_reco_growing_season']
LEAFOUT_VAR = ['leaf_out_90']
GPP_VARS = ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90']

MODEL_BLOCKS = {
    'A_environment_only': ENV_VARS,
    'B_environment_plus_leafout': ENV_VARS + LEAFOUT_VAR,
    'C_environment_plus_splitgpp': ENV_VARS + GPP_VARS,
    'D_environment_plus_leafout_plus_splitgpp': ENV_VARS + LEAFOUT_VAR + GPP_VARS,
}
LRT_PAIRS = [
    ('A_environment_only', 'B_environment_plus_leafout'),
    ('A_environment_only', 'C_environment_plus_splitgpp'),
    ('B_environment_plus_leafout', 'D_environment_plus_leafout_plus_splitgpp'),
    ('C_environment_plus_splitgpp', 'D_environment_plus_leafout_plus_splitgpp'),
    ('A_environment_only', 'D_environment_plus_leafout_plus_splitgpp'),
]

for p in (PREDICTORS_CSV, RATE_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run steps 1-7 first.")


def zscore(s):
    return (s - s.mean()) / s.std()


pred = pd.read_csv(PREDICTORS_CSV)
rate = pd.read_csv(RATE_CSV)[['site_id', 'year', 'vi_index',
                               'temperature_senescence_rate', 'photoperiod_senescence_rate']]
merged = pred.merge(rate, on=['site_id', 'year', 'vi_index'], how='inner')
print(f"{len(merged)} site-year-index rows available.\n")

comparison_rows, lrt_rows, cv_rows = [], [], []

for vi in VI_INDICES:
    print(f"\n{'=' * 70}\n{vi}\n{'=' * 70}")
    all_vars = sorted(set(ENV_VARS + LEAFOUT_VAR + GPP_VARS))
    d = merged[merged['vi_index'] == vi][[TARGET, 'site_id'] + all_vars].dropna().reset_index(drop=True)
    n_obs, n_sites = len(d), d['site_id'].nunique()
    if n_obs < MIN_OBS or n_sites < MIN_SITES:
        print(f"Skipping {vi}: only {n_obs} obs / {n_sites} sites.")
        continue
    print(f"n_obs={n_obs}, n_sites={n_sites}")

    for v in all_vars:
        d[v + '_z'] = zscore(d[v])

    fits_ml = {}
    for model_name, vars_ in MODEL_BLOCKS.items():
        formula = f"{TARGET} ~ " + " + ".join(v + '_z' for v in vars_)
        fit_ml = smf.mixedlm(formula, d, groups=d['site_id']).fit(reml=False)
        fits_ml[model_name] = fit_ml
        comparison_rows.append({'vi_index': vi, 'model': model_name, 'n_predictors': len(vars_),
                                'n_obs': n_obs, 'n_sites': n_sites,
                                'loglik': fit_ml.llf, 'aic': fit_ml.aic, 'bic': fit_ml.bic})

    print("\nModel fit (ML):")
    print(pd.DataFrame([r for r in comparison_rows if r['vi_index'] == vi])
          [['model', 'n_predictors', 'loglik', 'aic', 'bic']].to_string(index=False))

    print("\nPairwise likelihood-ratio tests:")
    for reduced_name, full_name in LRT_PAIRS:
        fit_reduced, fit_full = fits_ml[reduced_name], fits_ml[full_name]
        df_diff = len(MODEL_BLOCKS[full_name]) - len(MODEL_BLOCKS[reduced_name])
        stat = max(2 * (fit_full.llf - fit_reduced.llf), 0.0)
        p = chi2.sf(stat, df_diff)
        lrt_rows.append({'vi_index': vi, 'reduced_model': reduced_name, 'full_model': full_name,
                         'df': df_diff, 'chi2': stat, 'p_value': p})
        print(f"  {reduced_name:45s} -> {full_name:45s}  chi2={stat:6.2f}  df={df_diff}  p={p:.4g}")

    # Leave-one-site-out CV
    sites = d['site_id'].unique()
    preds = {name: np.full(len(d), np.nan) for name in MODEL_BLOCKS}
    for held_out_site in sites:
        train_mask = d['site_id'] != held_out_site
        test_mask = ~train_mask
        if test_mask.sum() == 0 or train_mask.sum() < MIN_OBS:
            continue
        for model_name, vars_ in MODEL_BLOCKS.items():
            train, test = d.loc[train_mask].copy(), d.loc[test_mask].copy()
            means = {v: train[v].mean() for v in vars_}
            stds = {v: train[v].std() for v in vars_}
            for v in vars_:
                train[v + '_zc'] = (train[v] - means[v]) / stds[v]
                test[v + '_zc'] = (test[v] - means[v]) / stds[v]
            formula_cv = f"{TARGET} ~ " + " + ".join(v + '_zc' for v in vars_)
            try:
                fit_cv = smf.mixedlm(formula_cv, train, groups=train['site_id']).fit(reml=True)
            except Exception:
                continue
            coef_names = ['Intercept'] + [v + '_zc' for v in vars_]
            X_test = np.column_stack([np.ones(len(test))] + [test[v + '_zc'].to_numpy() for v in vars_])
            coefs = np.array([fit_cv.params.get(c, 0.0) for c in coef_names])
            preds[model_name][test.index.to_numpy()] = X_test @ coefs

    y_true = d[TARGET].to_numpy()
    print("\nLeave-one-site-out CV:")
    for model_name, yhat in preds.items():
        mask = ~np.isnan(yhat)
        resid = y_true[mask] - yhat[mask]
        rmse = float(np.sqrt(np.mean(resid ** 2)))
        ss_res, ss_tot = float(np.sum(resid ** 2)), float(np.sum((y_true[mask] - y_true[mask].mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        cv_rows.append({'vi_index': vi, 'model': model_name, 'n_predicted': int(mask.sum()),
                        'loso_rmse_days': rmse, 'loso_r2': r2})
        print(f"  {model_name:45s} n={int(mask.sum()):4d}  RMSE={rmse:6.2f}d  R2={r2:.3f}")

pd.DataFrame(comparison_rows).to_csv(OUTPUT_COMPARISON_CSV, index=False)
pd.DataFrame(lrt_rows).to_csv(OUTPUT_LRT_CSV, index=False)
pd.DataFrame(cv_rows).to_csv(OUTPUT_CV_CSV, index=False)
print(f"\n\nWritten: '{OUTPUT_COMPARISON_CSV.name}', '{OUTPUT_LRT_CSV.name}', '{OUTPUT_CV_CSV.name}'")
