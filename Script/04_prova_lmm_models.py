import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Response variables: from prova_sigmoid_fit.py output. Each row there is one
# (site_id, year, vi_index) combination, so the same environmental/biological/
# biomass predictors get tested against phenology derived from every VI
# separately (NDVI, EVI, NIRv, ARI1, CRI1 - Hypotheses 1-4).
PHENOLOGY_CSV = DATA_DIR / "phenology_sos_eos_by_site_year_index.csv"

# Daily FLUXNET table (prova_fluxent_data.py output) - source of environmental
# predictors (Tair, VPD, ...). NOTE: precipitation (P_F) is not currently
# pulled by prova_fluxent_data.py's TARGET_MAPPING - add it there first if
# you want a precipitation term; it is silently skipped below if absent.
FLUX_CSV = DATA_DIR / "fluxnet_daily_selected_vars.csv"

# Site-year table of "biological" (canopy N, chlorophyll, PRI, water status)
# and "biomass" (GEDI/ICESat-2 height or S1 3D volume) predictors. This file
# does not yet exist in the pipeline - build it once those extractions are
# ready. Expected columns: site_id, year, plus whichever of BIO_PREDICTORS /
# BIOMASS_PREDICTORS below you have available. The script degrades
# gracefully (Model 1 only) if the file or specific columns are missing.
BIO_BIOMASS_CSV = DATA_DIR / "site_year_biological_biomass_predictors.csv"

OUTPUT_SUMMARY_CSV = DATA_DIR / "lmm_model_comparison_summary.csv"
OUTPUT_COEF_CSV = DATA_DIR / "lmm_fixed_effects_coefficients.csv"

RESPONSE_VARS = ['EOS10', 'growing_season_length']
GROUP_COL = 'site_id'

# Candidate predictors per tier - only the ones actually present (and not
# all-NaN) in the merged table are used, so this list can stay a superset.
ENV_CANDIDATES = ['TA_F_mean', 'VPD_F_mean', 'SW_IN_F_mean', 'P_F_sum']
BIO_CANDIDATES = ['PRI_mean', 'Chlorophyll_mean', 'CanopyN_mean', 'NDWI_mean']
BIOMASS_CANDIDATES = ['AGB', 'canopy_height_m', 'canopy_volume_m3']

MIN_GROUPS = 5           # minimum number of distinct sites for a random intercept to be meaningful
MIN_OBS_PER_MODEL = 20   # minimum rows for a given (response, vi_index, model) fit
STANDARDIZE_PREDICTORS = True  # z-score predictors so coefficients are comparable across tiers


# ---------------------------------------------------------------------------
# 1. Build environmental predictors from daily FLUXNET data
# ---------------------------------------------------------------------------
def aggregate_environmental(flux_csv):
    df = pd.read_csv(flux_csv)
    df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
    df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
    df['year'] = df['date'].dt.year

    agg_spec = {}
    if 'TA_F' in df.columns:
        agg_spec['TA_F'] = 'mean'
    if 'VPD_F' in df.columns:
        agg_spec['VPD_F'] = 'mean'
    if 'SW_IN_F' in df.columns:
        agg_spec['SW_IN_F'] = 'mean'
    if 'P_F' in df.columns:
        agg_spec['P_F'] = 'sum'
    else:
        print("Note: 'P_F' (precipitation) not found in FLUX_CSV - "
              "add it to TARGET_MAPPING in prova_fluxent_data.py to include it here.")

    env = df.groupby(['site_id', 'year']).agg(agg_spec).reset_index()
    rename_map = {'TA_F': 'TA_F_mean', 'VPD_F': 'VPD_F_mean', 'SW_IN_F': 'SW_IN_F_mean', 'P_F': 'P_F_sum'}
    env = env.rename(columns=rename_map)
    return env


# ---------------------------------------------------------------------------
# 2. Load biological / biomass predictors (optional)
# ---------------------------------------------------------------------------
def load_bio_biomass(path):
    if not os.path.exists(path):
        print(f"Note: '{path}' not found - only the environmental-only model "
              "will be fit. Add this file with the columns listed in "
              "BIO_CANDIDATES / BIOMASS_CANDIDATES to unlock Models 2 and 3.")
        return None
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# 3. Mixed model fitting + likelihood-ratio comparison
# ---------------------------------------------------------------------------
def fit_mixedlm(data, response, predictors, group_col):
    """Fit response ~ predictors + (1 | group_col) by ML (not REML, so
    likelihood-ratio tests between nested fixed-effects sets are valid)."""
    fixed = ' + '.join(predictors) if predictors else '1'
    formula = f"{response} ~ {fixed}"
    model = smf.mixedlm(formula, data=data, groups=data[group_col])

    # The default optimizer is more robust here than lbfgs, which tends to
    # hit a singular-matrix error when the random-intercept variance is
    # small; fall back to powell if the default also fails to converge.
    for method in (None, 'powell', 'cg'):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = model.fit(reml=False) if method is None else model.fit(reml=False, method=method)
            return result
        except Exception:
            continue

    print(f"   -> Fit failed for '{formula}' with all optimizers.")
    return None


def lrt_compare(simple_result, complex_result):
    """Likelihood-ratio (ANOVA-style) test between two nested mixed models
    fit on the SAME data with ML. Returns (lr_stat, df_diff, p_value)."""
    if simple_result is None or complex_result is None:
        return np.nan, np.nan, np.nan
    df_diff = complex_result.model.k_fe - simple_result.model.k_fe
    if df_diff <= 0:
        return np.nan, np.nan, np.nan
    lr_stat = 2 * (complex_result.llf - simple_result.llf)
    lr_stat = max(lr_stat, 0.0)  # guard against tiny negative values from optimizer noise
    p_value = chi2.sf(lr_stat, df_diff)
    return lr_stat, df_diff, p_value


def pseudo_r2(result):
    """Nakagawa & Schielzeth marginal (fixed effects only) and conditional
    (fixed + random) R^2 for a fitted MixedLM result."""
    try:
        fe_pred = result.model.exog.dot(result.fe_params)
        var_fixed = np.var(fe_pred)
        var_random = float(result.cov_re.iloc[0, 0]) if result.cov_re.size else 0.0
        var_resid = result.scale
        denom = var_fixed + var_random + var_resid
        r2_marginal = var_fixed / denom if denom > 0 else np.nan
        r2_conditional = (var_fixed + var_random) / denom if denom > 0 else np.nan
        return r2_marginal, r2_conditional
    except Exception:
        return np.nan, np.nan


def zscore(series):
    std = series.std()
    if not np.isfinite(std) or std == 0:
        return series - series.mean()
    return (series - series.mean()) / std


# ---------------------------------------------------------------------------
# 4. Build the merged site-year table
# ---------------------------------------------------------------------------
if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run prova_sigmoid_fit.py first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run prova_fluxent_data.py first.")

pheno = pd.read_csv(PHENOLOGY_CSV)
env = aggregate_environmental(FLUX_CSV)
bio_biomass = load_bio_biomass(BIO_BIOMASS_CSV)

merged = pheno.merge(env, on=['site_id', 'year'], how='left')
if bio_biomass is not None:
    merged = merged.merge(bio_biomass, on=['site_id', 'year'], how='left')

env_predictors = [c for c in ENV_CANDIDATES if c in merged.columns and merged[c].notna().any()]
bio_predictors = [c for c in BIO_CANDIDATES if c in merged.columns and merged[c].notna().any()]
biomass_predictors = [c for c in BIOMASS_CANDIDATES if c in merged.columns and merged[c].notna().any()]

print(f"Environmental predictors available: {env_predictors}")
print(f"Biological predictors available:    {bio_predictors}")
print(f"Biomass predictors available:       {biomass_predictors}\n")

if STANDARDIZE_PREDICTORS:
    for col in env_predictors + bio_predictors + biomass_predictors:
        merged[col] = zscore(merged[col])

# ---------------------------------------------------------------------------
# 5. Fit the three nested models per (response_var, vi_index)
# ---------------------------------------------------------------------------
summary_rows = []
coef_rows = []

vi_indices = sorted(merged['vi_index'].dropna().unique())

for response_var in RESPONSE_VARS:
    for vi_index in vi_indices:
        sub_all = merged[(merged['vi_index'] == vi_index) & merged[response_var].notna()].copy()
        if sub_all[GROUP_COL].nunique() < MIN_GROUPS:
            print(f"Skipping {response_var} / {vi_index}: fewer than {MIN_GROUPS} sites available.")
            continue

        tiers = {
            'env_only': env_predictors,
            'env_bio': env_predictors + bio_predictors,
            'env_bio_biomass': env_predictors + bio_predictors + biomass_predictors,
        }

        # For a fair, apples-to-apples LRT, fit all three tiers on the same
        # complete-case subset (rows with every predictor used by the most
        # complex model available). If that subset is too small, still fit
        # env_only on its own maximal data, but skip the nested comparison.
        full_predictor_set = tiers['env_bio_biomass']
        complete_cols = [response_var] + full_predictor_set
        complete_sub = sub_all.dropna(subset=complete_cols) if full_predictor_set else sub_all.dropna(subset=[response_var])

        results = {}
        if len(complete_sub) >= MIN_OBS_PER_MODEL and complete_sub[GROUP_COL].nunique() >= MIN_GROUPS:
            fit_data = complete_sub
            comparable = True
        else:
            fit_data = sub_all.dropna(subset=[response_var] + env_predictors)
            comparable = False
            print(f"Note: {response_var} / {vi_index} - not enough complete-case rows for a fair "
                  "3-model LRT comparison; fitting env_only on its own maximal subset instead.")

        for tier_name, predictors in tiers.items():
            if not predictors and tier_name != 'env_only':
                continue  # no bio/biomass predictors available - nothing new to fit
            data_for_tier = fit_data if comparable else (fit_data if tier_name == 'env_only' else None)
            if data_for_tier is None or len(data_for_tier) < MIN_OBS_PER_MODEL:
                continue
            result = fit_mixedlm(data_for_tier, response_var, predictors, GROUP_COL)
            results[tier_name] = result

            if result is not None:
                r2m, r2c = pseudo_r2(result)
                for pname, est, se, z, p in zip(result.params.index, result.params, result.bse, result.tvalues, result.pvalues):
                    coef_rows.append({
                        'response_var': response_var, 'vi_index': vi_index, 'model': tier_name,
                        'predictor': pname, 'estimate': est, 'se': se, 'z': z, 'p_value': p,
                        'n_obs': len(data_for_tier), 'n_sites': data_for_tier[GROUP_COL].nunique(),
                        'r2_marginal': r2m, 'r2_conditional': r2c
                    })

        lr1, df1, p1 = lrt_compare(results.get('env_only'), results.get('env_bio')) if comparable else (np.nan, np.nan, np.nan)
        lr2, df2, p2 = lrt_compare(results.get('env_bio'), results.get('env_bio_biomass')) if comparable else (np.nan, np.nan, np.nan)

        summary_rows.append({
            'response_var': response_var, 'vi_index': vi_index,
            'n_obs': len(fit_data), 'n_sites': fit_data[GROUP_COL].nunique(),
            'complete_case_comparison': comparable,
            'AIC_env_only': results['env_only'].aic if results.get('env_only') else np.nan,
            'AIC_env_bio': results['env_bio'].aic if results.get('env_bio') else np.nan,
            'AIC_env_bio_biomass': results['env_bio_biomass'].aic if results.get('env_bio_biomass') else np.nan,
            'LRT_env_vs_env+bio_stat': lr1, 'LRT_env_vs_env+bio_df': df1, 'LRT_env_vs_env+bio_p': p1,
            'LRT_env+bio_vs_env+bio+biomass_stat': lr2, 'LRT_env+bio_vs_env+bio+biomass_df': df2,
            'LRT_env+bio_vs_env+bio+biomass_p': p2,
        })

        print(f"{response_var} / {vi_index}: n={len(fit_data)}, sites={fit_data[GROUP_COL].nunique()}, "
              f"env-vs-env+bio p={p1:.4g}, env+bio-vs-+biomass p={p2:.4g}" if comparable else
              f"{response_var} / {vi_index}: n={len(fit_data)}, sites={fit_data[GROUP_COL].nunique()} (env_only only)")

# ---------------------------------------------------------------------------
# 6. Save outputs
# ---------------------------------------------------------------------------
pd.DataFrame(summary_rows).to_csv(OUTPUT_SUMMARY_CSV, index=False)
pd.DataFrame(coef_rows).to_csv(OUTPUT_COEF_CSV, index=False)

print(f"\nModel comparison summary (AIC + ANOVA/LRT p-values) written to '{OUTPUT_SUMMARY_CSV}'.")
print(f"Fixed-effects coefficient tables written to '{OUTPUT_COEF_CSV}'.")