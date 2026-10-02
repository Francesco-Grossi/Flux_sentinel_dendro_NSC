"""
PIPELINE STEP 28 - Build ONE consolidated site-year-index predictor table
for every EOS source - satellite NDVI and NIRv (step 23), tower NDVI (step 24)
and PhenoCam GCC (step 25): spring/autumn phenology timing, cumulative-GPP
windows split around the summer solstice, and
environmental covariates (radiation, mean temperature, respiration). Also
writes an exploratory correlation table (autumn parameter x predictor, pooled
across sites). Step 40 plots every pair 1:1; the legacy hypothesis tests
(Script/legacy/08, 09) read the same table.

GPP windows computed:
    gpp_sos10_to_solstice   - spring green-up onset (SOS10) -> summer solstice
    gpp_solstice_to_eos10   - solstice -> EOS10 (near-dormant; "full autumn decline")
    gpp_solstice_to_eos90   - solstice -> EOS90 (still ~90% green; "early decline only")
                              This is the window the split-GPP (opposite-effect)
                              hypothesis is built on.
    gpp_solstice_to_eos10_fixed, gpp_solstice_to_eos90_fixed
                            - the same two windows, but ending at the SITE'S MEAN
                              EOS (same end date every year, >= 3 years).

Why both: a window that ends at the same year's EOS gets longer when EOS is
later, so its cumulative GPP rises with EOS by construction. Comparing the
1:1 plots of the year-anchored and the fixed window shows how much of the
positive post-solstice relation is that length effect and how much is real.

Growing-season totals and means (GPP, Reco, radiation, temperature) run from
the year's SOS10 to the SITE'S MEAN EOS10, for the same reason: ending them at
the same year's EOS10 made a late EOS add extra (cold, dark, low-GPP) autumn
days, which produced strong correlations with EOS by construction. The
same-year versions are kept in the table as <name>_same_year, but are not in
the correlation table or the 1:1 plots.

Output: data/phenology_flux_predictors_by_site_year_index.csv
        data/autumn_phenology_correlations.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # step 23 output
# tower NDVI (step 24) and PhenoCam (step 25) phenology, same columns; used when present
EXTRA_PHENOLOGY_CSVS = [DATA_DIR / "phenology_tower_by_site_year_index.csv",
                        DATA_DIR / "phenology_phenocam_by_site_year_index.csv"]
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                            # step 22 output

OUTPUT_SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"
OUTPUT_CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"

VI_INDICES = ['NDVI', 'NIRv', 'NDVI_tower', 'GCC']
AUTUMN_PARAMS = ['EOS90', 'EOS50', 'EOS10', 'senescence_kinetic_i']
SPRING_PARAMS = ['leaf_out_10', 'leaf_out_50', 'leaf_out_90']
MIN_PAIRS_FOR_CORR = 5
MIN_FIT_CORR = 0.8  # quality gate on the double-logistic curve fit itself

FLUX_VARS = {
    'gpp': {'column': 'GPP_NT_VUT_REF', 'plausible_range': (-5, 50), 'agg': 'sum'},
    'reco': {'column': 'RECO_NT_VUT_REF', 'plausible_range': (-5, 40), 'agg': 'sum'},
    'radiation': {'column': 'SW_IN_F', 'plausible_range': (0, 500), 'agg': 'mean'},
    'temperature': {'column': 'TA_F', 'plausible_range': (-60, 50), 'agg': 'mean'},
}

if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run 23_phenology_satellite.py first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run 22_merge_fluxnet_hls.py first.")

pheno = pd.concat([pd.read_csv(p) for p in [PHENOLOGY_CSV] + EXTRA_PHENOLOGY_CSVS if os.path.exists(p)],
                  ignore_index=True)
VI_INDICES = [v for v in VI_INDICES if v in set(pheno['vi_index'])]
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')
              & (pheno['corr'] >= MIN_FIT_CORR)].copy()
if 'qc_pass' in pheno.columns:  # step-23 QC flags (autumn coverage, R2, date order, ...)
    pheno = pheno[pheno['qc_pass'].astype(bool)].copy()
print(f"Phenology rows available (successful fits, corr >= {MIN_FIT_CORR}, QC passed): {len(pheno)}")
MIN_YEARS_ANCHOR = 3
for lvl in ('EOS90', 'EOS10'):
    grp = pheno.groupby(['site_id', 'vi_index'])[lvl]
    pheno[f'{lvl}_site_mean'] = grp.transform('mean').where(grp.transform('count') >= MIN_YEARS_ANCHOR)

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
        raise ValueError(f"'{col}' not found in FLUX_CSV.")
    implausible = (flux[col] < lo) | (flux[col] > hi)
    if implausible.any():
        print(f"Excluding {int(implausible.sum())} implausible '{col}' value(s) outside [{lo}, {hi}].")
        flux.loc[implausible, col] = np.nan
    flux_lookups[var_name] = {key: group[['doy', col]].dropna(subset=['doy'])
                               for key, group in flux.groupby(['site_id', 'year'])}


def agg_flux_var(var_name, site_id, year, doy_start, doy_end):
    if pd.isna(doy_start) or pd.isna(doy_end) or doy_start > doy_end:
        return np.nan
    g = flux_lookups[var_name].get((site_id, year))
    if g is None or g.empty:
        return np.nan
    col = FLUX_VARS[var_name]['column']
    mask = (g['doy'] >= doy_start) & (g['doy'] <= doy_end)
    vals = g.loc[mask, col].dropna()
    if vals.empty:
        return np.nan
    agg = FLUX_VARS[var_name]['agg']
    return float(vals.sum()) if agg == 'sum' else float(vals.mean())


def solstice_doy(year):
    return pd.Timestamp(year=year, month=6, day=21).dayofyear


def photoperiod_hours(doy, lat_deg):
    if pd.isna(doy) or pd.isna(lat_deg):
        return np.nan
    lat_rad = np.radians(lat_deg)
    declination_rad = np.radians(-23.44 * np.cos(np.radians(360.0 / 365.0 * (doy + 10))))
    cos_hour_angle = np.clip(-np.tan(lat_rad) * np.tan(declination_rad), -1.0, 1.0)
    return float((24.0 / np.pi) * np.arccos(cos_hour_angle))


site_lat = flux.groupby('site_id')['lat'].first().to_dict()

records = []
for _, row in pheno.iterrows():
    site_id, year = row['site_id'], int(row['year'])
    sos10, sos50 = row.get('leaf_out_10'), row.get('leaf_out_50')
    eos50, eos10, eos90 = row.get('EOS50'), row.get('EOS10'), row.get('EOS90')

    growing_season_length = (eos50 - sos50) if pd.notna(eos50) and pd.notna(sos50) else np.nan
    sol_doy = solstice_doy(year)

    rec = row.to_dict()
    rec['solstice_doy'] = sol_doy
    rec['growing_season_length'] = growing_season_length

    for var_name in FLUX_VARS:
        agg = FLUX_VARS[var_name]['agg']
        label = 'total' if agg == 'sum' else 'mean'
        rec[f'{label}_{var_name}_growing_season'] = agg_flux_var(var_name, site_id, year, sos10,
                                                                 row.get('EOS10_site_mean'))
        rec[f'{label}_{var_name}_growing_season_same_year'] = agg_flux_var(var_name, site_id, year, sos10, eos10)

    # Split GPP windows around the summer solstice - the core predictors for
    # the split-GPP (opposite-effect) hypothesis.
    rec['gpp_sos10_to_solstice'] = agg_flux_var('gpp', site_id, year, sos10, sol_doy)
    rec['gpp_solstice_to_eos10'] = agg_flux_var('gpp', site_id, year, sol_doy, eos10)
    rec['gpp_solstice_to_eos90'] = agg_flux_var('gpp', site_id, year, sol_doy, eos90)
    # same windows with a fixed end date (site mean EOS): no window-length effect
    rec['gpp_solstice_to_eos10_fixed'] = agg_flux_var('gpp', site_id, year, sol_doy, row.get('EOS10_site_mean'))
    rec['gpp_solstice_to_eos90_fixed'] = agg_flux_var('gpp', site_id, year, sol_doy, row.get('EOS90_site_mean'))

    rec['photoperiod_at_EOS90'] = photoperiod_hours(eos90, site_lat.get(site_id))
    records.append(rec)

analysis_df = pd.DataFrame(records)
analysis_df.to_csv(OUTPUT_SITEYEAR_CSV, index=False)
print(f"Predictor table -> '{OUTPUT_SITEYEAR_CSV}' ({len(analysis_df)} rows).")

# Exploratory pooled correlations (autumn parameter x predictor) - a first
# look before the site-controlled models (steps 33-38).
FLUX_PREDICTORS = ['total_gpp_growing_season', 'gpp_sos10_to_solstice', 'gpp_solstice_to_eos10',
                    'gpp_solstice_to_eos90', 'gpp_solstice_to_eos10_fixed', 'gpp_solstice_to_eos90_fixed',
                    'total_reco_growing_season', 'mean_radiation_growing_season',
                    'mean_temperature_growing_season', 'photoperiod_at_EOS90']
PREDICTORS = SPRING_PARAMS + ['growing_season_length'] + FLUX_PREDICTORS

corr_rows = []
for vi in VI_INDICES:
    sub = analysis_df[analysis_df['vi_index'] == vi]
    for autumn_p in AUTUMN_PARAMS:
        for pred in PREDICTORS:
            pair = sub[[autumn_p, pred]].dropna()
            if len(pair) < MIN_PAIRS_FOR_CORR or pair[autumn_p].std() == 0 or pair[pred].std() == 0:
                continue
            r, p = pearsonr(pair[pred], pair[autumn_p])
            corr_rows.append({'vi_index': vi, 'autumn_parameter': autumn_p, 'predictor': pred,
                               'n': len(pair), 'pearson_r': r, 'p_value': p})

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUTPUT_CORR_CSV, index=False)
print(f"Correlation table -> '{OUTPUT_CORR_CSV}' ({len(corr_df)} rows).")
