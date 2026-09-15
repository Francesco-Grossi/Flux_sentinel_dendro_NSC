import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

PREDICTORS_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # script 6 output
CLUSTERS_CSV = DATA_DIR / "site_year_climate_clusters.csv"                      # script 10 output (has is_drought)

OUTPUT_GRID_CSV = DATA_DIR / "drought_vs_nondrought_correlations.csv"       # every (group x pair) tested
OUTPUT_COMPARISON_CSV = DATA_DIR / "drought_vs_nondrought_comparison.csv"   # drought and non-drought r side by side
OUTPUT_SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_with_drought_flag.csv"  # site-year rows + is_drought

VI_INDICES = ['NDVI', 'NIRv']
AUTUMN_PARAMS = ['EOS90', 'EOS50', 'EOS10', 'senescence_kinetic_i']
SPRING_PARAMS = ['leaf_out_10', 'leaf_out_50', 'leaf_out_90', 'greenup_kinetic_i']
FLUX_PREDICTORS = ['total_gpp_growing_season', 'gpp_sos10_to_solstice', 'gpp_solstice_to_eos10',
                    'total_reco_growing_season', 'mean_radiation_growing_season',
                    'mean_temperature_growing_season', 'photoperiod_at_EOS90']
PREDICTORS = SPRING_PARAMS + ['growing_season_length'] + FLUX_PREDICTORS

# Smaller than script 6's pooled threshold, since splitting into 2 groups
# (and drought years are the minority by construction) leaves fewer
# site-years per group than the full pooled set.
MIN_PAIRS_FOR_CORR = 5

# Quality gate on the Landsat phenology curve fit - see script 6. Applied
# again here defensively, in case PREDICTORS_CSV predates that filter.
MIN_FIT_CORR = 0.8

if not os.path.exists(PREDICTORS_CSV):
    raise FileNotFoundError(f"Missing '{PREDICTORS_CSV}'. Run the autumn-phenology-correlation script first.")
if not os.path.exists(CLUSTERS_CSV):
    raise FileNotFoundError(f"Missing '{CLUSTERS_CSV}'. Run the climate-clustering script first.")


# ---------------------------------------------------------------------------
# 1. Load + merge - each (site_id, year, vi_index) row gets its drought
#    flag attached (from the climate-clustering script's within-site
#    precipitation/VPD anomaly). Site-years without a flag (e.g. dropped
#    upstream for insufficient climate-variable coverage) are excluded via
#    the inner join.
# ---------------------------------------------------------------------------
predictors_df = pd.read_csv(PREDICTORS_CSV)
predictors_df['year'] = predictors_df['year'].astype(int)
if 'corr' in predictors_df.columns:
    n_before = len(predictors_df)
    predictors_df = predictors_df[predictors_df['corr'] >= MIN_FIT_CORR].copy()
    if len(predictors_df) < n_before:
        print(f"Dropped {n_before - len(predictors_df)} row(s) with phenology fit corr < {MIN_FIT_CORR}.")

clusters_df = pd.read_csv(CLUSTERS_CSV)
if 'is_drought' not in clusters_df.columns:
    raise ValueError(f"'is_drought' not found in '{CLUSTERS_CSV}'. Re-run the climate-clustering script "
                      "(the drought-anomaly feature) first.")
clusters_df = clusters_df[['site_id', 'year', 'is_drought']].copy()
clusters_df['year'] = clusters_df['year'].astype(int)

merged = predictors_df.merge(clusters_df, on=['site_id', 'year'], how='inner')
n_dropped = len(predictors_df) - len(merged)
print(f"{len(merged)} of {len(predictors_df)} phenology rows have a drought-flag assignment "
      f"({n_dropped} dropped - no flag available for that site-year).")

n_drought = int((merged.drop_duplicates(['site_id', 'year'])['is_drought']).sum())
n_total_siteyears = merged.drop_duplicates(['site_id', 'year']).shape[0]
print(f"Drought site-years: {n_drought} / {n_total_siteyears}.")

# Save the merged site-year-index rows (predictors + autumn params +
# is_drought) for downstream plotting - e.g. scatter plots split by
# drought/non-drought, which need the actual data points, not just the
# aggregate r/p/n this script otherwise produces.
merged.to_csv(OUTPUT_SITEYEAR_CSV, index=False)
print(f"Merged site-year table (with drought flag) written to '{OUTPUT_SITEYEAR_CSV}'.")

GROUPS = {'drought': True, 'non_drought': False}

# ---------------------------------------------------------------------------
# 2. Correlate each autumn parameter against each predictor, SEPARATELY for
#    drought and non-drought site-years, per VI index.
# ---------------------------------------------------------------------------
results = []
for vi in VI_INDICES:
    vi_sub = merged[merged['vi_index'] == vi]
    for group_name, flag_value in GROUPS.items():
        group_sub = vi_sub[vi_sub['is_drought'] == flag_value]
        for autumn_p in AUTUMN_PARAMS:
            if autumn_p not in group_sub.columns:
                continue
            for pred in PREDICTORS:
                if pred not in group_sub.columns:
                    continue
                pair = group_sub[[pred, autumn_p]].dropna()
                if len(pair) < MIN_PAIRS_FOR_CORR or pair[pred].std() == 0 or pair[autumn_p].std() == 0:
                    continue
                r, p_value = pearsonr(pair[pred], pair[autumn_p])
                results.append({
                    'vi_index': vi, 'group': group_name,
                    'autumn_parameter': autumn_p, 'predictor': pred,
                    'n': len(pair), 'pearson_r': r, 'p_value': p_value,
                })

scan_df = pd.DataFrame(results)
scan_df.to_csv(OUTPUT_GRID_CSV, index=False)
print(f"\nDrought vs non-drought correlations ({len(scan_df):,} rows) written to '{OUTPUT_GRID_CSV}'.")

# ---------------------------------------------------------------------------
# 3. Side-by-side comparison table: drought r/p/n next to non-drought r/p/n
#    for the same (vi_index, autumn_parameter, predictor), plus the
#    difference in r - the quickest way to see where drought changes a
#    relationship's strength or even flips its sign.
# ---------------------------------------------------------------------------
pivot = scan_df.pivot_table(
    index=['vi_index', 'autumn_parameter', 'predictor'],
    columns='group', values=['pearson_r', 'p_value', 'n'], aggfunc='first'
)
pivot.columns = [f'{stat}_{group}' for stat, group in pivot.columns]
pivot = pivot.reset_index()

if 'pearson_r_drought' in pivot.columns and 'pearson_r_non_drought' in pivot.columns:
    pivot['r_difference'] = pivot['pearson_r_drought'] - pivot['pearson_r_non_drought']

pivot.to_csv(OUTPUT_COMPARISON_CSV, index=False)
print(f"Side-by-side drought/non-drought comparison written to '{OUTPUT_COMPARISON_CSV}'.")

if 'r_difference' in pivot.columns:
    top = pivot.reindex(pivot['r_difference'].abs().sort_values(ascending=False).index).head(10)
    print("\nPairs where drought changes the correlation most (by |r_drought - r_non_drought|):")
    cols = ['vi_index', 'autumn_parameter', 'predictor',
            'pearson_r_drought', 'pearson_r_non_drought', 'r_difference']
    print(top[[c for c in cols if c in top.columns]].to_string(index=False))