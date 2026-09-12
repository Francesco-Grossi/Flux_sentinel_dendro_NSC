import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # script 5 output
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"  # script 4 output (full daily FLUXNET + sparse VI)

OUTPUT_SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"
OUTPUT_CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"

# Only NDVI and NIRv, as requested - EVI is intentionally excluded here.
VI_INDICES = ['NDVI', 'NIRv']

# The 4 parameters describing autumn/senescence phenology (script 5 output).
AUTUMN_PARAMS = ['senescence_10', 'senescence_50', 'senescence_90', 'senescence_kinetic_i']

# The 4 parameters describing spring/green-up phenology (script 5 output) -
# one of the predictor groups each autumn parameter is compared against.
SPRING_PARAMS = ['leaf_out_10', 'leaf_out_50', 'leaf_out_90', 'greenup_kinetic_i']

# Minimum number of non-NaN (predictor, autumn-parameter) pairs required
# before a correlation is computed/trusted - guards against a handful of
# site-years producing a meaningless "perfect" correlation.
MIN_PAIRS_FOR_CORR = 5

# FLUXNET columns used for the flux-based predictors below, and the
# physically plausible daily range (gC m-2 d-1) for each - anything outside
# this range is almost certainly a leftover fill/sentinel value (e.g.
# -9999) rather than real data, and gets excluded before summing so a
# single bad value can't blow up a whole site-year's total.
#
# GPP_NT_VUT_REF is GPP estimated via the NIGHTTIME partitioning method
# (Reichstein et al. 2005, using nighttime NEE to model respiration and
# subtracting it from NEE across the full day) - GPP_DT_VUT_REF (daytime
# partitioning, Lasslop et al. 2010) was tried here too, but produced
# site-years with total_gpp_growing_season in the millions (an unblanked
# -9999-style sentinel), so this stays on GPP_NT_VUT_REF.
#
# RECO_NT_VUT_REF is TOTAL ECOSYSTEM respiration (autotrophic + heterotrophic
# combined) from the same nighttime-partitioning method - it's what FLUXNET
# towers actually measure/partition. Autotrophic respiration (Ra) alone is
# NOT a standard FLUXNET variable (eddy covariance doesn't separate it from
# heterotrophic respiration without additional chamber/biometric data or a
# partitioning model), so it's intentionally not included here.
FLUX_VARS = {
    'gpp': {'column': 'GPP_NT_VUT_REF', 'plausible_range': (-5, 50)},
    'reco': {'column': 'RECO_NT_VUT_REF', 'plausible_range': (-5, 40)},
}

if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run the double-logistic phenology script first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run the FLUXNET + Landsat merge script first.")


# ---------------------------------------------------------------------------
# 1. Load phenology results - only successful fits, only NDVI/NIRv
# ---------------------------------------------------------------------------
# 'fit_failed' rows carry NaN for every phenology parameter (see the
# double-logistic script), so they can't contribute to a correlation and
# are dropped here rather than silently producing all-NaN pairs downstream.
pheno = pd.read_csv(PHENOLOGY_CSV)
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')].copy()
print(f"Phenology rows available for correlation (NDVI/NIRv, successful fits only): {len(pheno)}")

# ---------------------------------------------------------------------------
# 2. Load daily FLUXNET data (full year, not just growing season), clean
#    each flux variable, and index it by (site_id, year) for fast repeated
#    window-sum lookups.
# ---------------------------------------------------------------------------
flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year
flux['doy'] = flux['date'].dt.dayofyear

flux_lookups = {}
for var_name, spec in FLUX_VARS.items():
    col, (lo, hi) = spec['column'], spec['plausible_range']
    if col not in flux.columns:
        raise ValueError(f"'{col}' not found in FLUX_CSV - can't compute the {var_name.upper()}-based "
                          "predictors. Run the FLUXNET extraction script with this column included first.")

    implausible = (flux[col] < lo) | (flux[col] > hi)
    n_implausible = int(implausible.sum())
    if n_implausible:
        print(f"Excluding {n_implausible} implausible '{col}' value(s) outside [{lo}, {hi}] gC m-2 d-1 "
              "(likely leftover fill/sentinel values) before summing.")
        flux.loc[implausible, col] = np.nan

    flux_lookups[var_name] = {
        key: group[['doy', col]].dropna(subset=['doy'])
        for key, group in flux.groupby(['site_id', 'year'])
    }


def sum_flux_var(var_name, site_id, year, doy_start, doy_end):
    """Sum of a daily flux variable (GPP or RECO, per FLUX_VARS) over
    [doy_start, doy_end] (inclusive) for one site-year. NaN if either
    boundary is undefined, the window is inverted, or there's simply no
    flux record for that site-year."""
    if pd.isna(doy_start) or pd.isna(doy_end) or doy_start > doy_end:
        return np.nan
    g = flux_lookups[var_name].get((site_id, year))
    if g is None or g.empty:
        return np.nan
    col = FLUX_VARS[var_name]['column']
    mask = (g['doy'] >= doy_start) & (g['doy'] <= doy_end)
    vals = g.loc[mask, col]
    if vals.dropna().empty:
        return np.nan
    return float(vals.sum(skipna=True))


def solstice_doy(year):
    """DOY of June 21 for the given calendar year - a fixed-date proxy for
    the summer solstice (the exact astronomical moment shifts by up to a
    day depending on year/leap year, which a fixed calendar date already
    accounts for via dayofyear, unlike a hardcoded DOY like 172 would)."""
    return pd.Timestamp(year=year, month=6, day=21).dayofyear


# ---------------------------------------------------------------------------
# 3. Build the per-site-year-index predictor table
# ---------------------------------------------------------------------------
# Note: the window boundaries (SOS10, EOS90, etc.) come from each row's OWN
# vi_index fit, so the GPP/RECO sums below are computed on the same
# underlying daily FLUXNET series but with NDVI-derived and NIRv-derived
# windows for their respective rows - they are not expected to be identical
# between the two indices for the same site-year.
records = []
for _, row in pheno.iterrows():
    site_id, year = row['site_id'], int(row['year'])
    sos10, sos50 = row.get('leaf_out_10'), row.get('leaf_out_50')
    eos50, eos90 = row.get('senescence_50'), row.get('senescence_90')

    growing_season_length = (eos50 - sos50) if pd.notna(eos50) and pd.notna(sos50) else np.nan
    sol_doy = solstice_doy(year)

    rec = row.to_dict()
    rec['solstice_doy'] = sol_doy
    rec['growing_season_length'] = growing_season_length

    for var_name in FLUX_VARS:
        rec[f'total_{var_name}_growing_season'] = sum_flux_var(var_name, site_id, year, sos10, eos90)
        # Only GPP gets split into the pre-/post-solstice sub-windows for
        # now - RECO is added as a single growing-season total predictor.
        if var_name == 'gpp':
            rec[f'{var_name}_sos10_to_solstice'] = sum_flux_var(var_name, site_id, year, sos10, sol_doy)
            rec[f'{var_name}_solstice_to_eos90'] = sum_flux_var(var_name, site_id, year, sol_doy, eos90)

    records.append(rec)

analysis_df = pd.DataFrame(records)
analysis_df.to_csv(OUTPUT_SITEYEAR_CSV, index=False)
print(f"Site-year-index predictor table written to '{OUTPUT_SITEYEAR_CSV}' ({len(analysis_df)} rows).")

# ---------------------------------------------------------------------------
# 4. Correlate each autumn parameter against each predictor, per VI index
# ---------------------------------------------------------------------------
# Predictors: the 4 spring parameters, growing-season length (EOS50-SOS50),
# the 3 GPP-based sub-window metrics, and total RECO over the growing
# season as a single additional predictor.
FLUX_PREDICTORS = ['total_gpp_growing_season', 'gpp_sos10_to_solstice', 'gpp_solstice_to_eos90',
                    'total_reco_growing_season']

PREDICTORS = SPRING_PARAMS + ['growing_season_length'] + FLUX_PREDICTORS

# NOTE: these correlations pool site-years across ALL sites together, so
# they mix between-site differences (e.g. a boreal vs. temperate site
# simply senescing at different DOYs) with within-site year-to-year
# variation. That's the straightforward reading of "find a correlation
# between the parameters" - if what's actually wanted is each site's own
# year-to-year relationship (removing the between-site offset), that's a
# per-site version of the same loop, or a mixed model with a site random
# intercept (as used elsewhere in this pipeline).
corr_rows = []
for vi in VI_INDICES:
    sub = analysis_df[analysis_df['vi_index'] == vi]
    for autumn_p in AUTUMN_PARAMS:
        for pred in PREDICTORS:
            pair = sub[[autumn_p, pred]].dropna()
            if len(pair) < MIN_PAIRS_FOR_CORR or pair[autumn_p].std() == 0 or pair[pred].std() == 0:
                continue
            r, p = pearsonr(pair[pred], pair[autumn_p])
            corr_rows.append({
                'vi_index': vi, 'autumn_parameter': autumn_p, 'predictor': pred,
                'n': len(pair), 'pearson_r': r, 'p_value': p
            })

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUTPUT_CORR_CSV, index=False)
print(f"Correlation table written to '{OUTPUT_CORR_CSV}' ({len(corr_df)} rows).")

if not corr_df.empty:
    n_sig = int((corr_df['p_value'] < 0.05).sum())
    print(f"\n{n_sig}/{len(corr_df)} autumn-parameter x predictor correlations are significant at p<0.05.")
    top = corr_df.reindex(corr_df['pearson_r'].abs().sort_values(ascending=False).index).head(10)
    print("\nStrongest correlations (by |r|):")
    print(top[['vi_index', 'autumn_parameter', 'predictor', 'n', 'pearson_r', 'p_value']].to_string(index=False))