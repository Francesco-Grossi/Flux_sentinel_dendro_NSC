"""
PIPELINE STEP 6 - Build ONE consolidated site-year-index predictor table:
spring/autumn phenology timing (from step 5), cumulative-GPP windows split
around the summer solstice, and environmental covariates (radiation, mean
temperature, respiration) - everything the hypothesis-test scripts (08, 09)
need, in a single file. Also writes a general exploratory correlation table
(autumn parameter x predictor, pooled across sites).

GPP windows computed:
    gpp_sos10_to_solstice   - spring green-up onset (SOS10) -> summer solstice
    gpp_solstice_to_eos10   - solstice -> EOS10 (near-dormant; "full autumn decline")
    gpp_solstice_to_eos90   - solstice -> EOS90 (still ~90% green; "early decline only")
                              This is the window the opposite-effect hypothesis
                              (script 08/09) is built on.

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
PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # step 5 output
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                            # step 4 output

OUTPUT_SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"
OUTPUT_CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"

VI_INDICES = ['NDVI', 'NIRv']
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
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run 05_double_logistic_phenology.py first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run 04_merge_fluxnet_landsat.py first.")

pheno = pd.read_csv(PHENOLOGY_CSV)
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')
              & (pheno['corr'] >= MIN_FIT_CORR)].copy()
print(f"Phenology rows available (successful fits, corr >= {MIN_FIT_CORR}): {len(pheno)}")

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
        rec[f'{label}_{var_name}_growing_season'] = agg_flux_var(var_name, site_id, year, sos10, eos10)

    # Split GPP windows around the summer solstice - the core predictors for
    # the opposite-effect hypothesis (scripts 08/09).
    rec['gpp_sos10_to_solstice'] = agg_flux_var('gpp', site_id, year, sos10, sol_doy)
    rec['gpp_solstice_to_eos10'] = agg_flux_var('gpp', site_id, year, sol_doy, eos10)
    rec['gpp_solstice_to_eos90'] = agg_flux_var('gpp', site_id, year, sol_doy, eos90)

    rec['photoperiod_at_EOS90'] = photoperiod_hours(eos90, site_lat.get(site_id))
    records.append(rec)

analysis_df = pd.DataFrame(records)
analysis_df.to_csv(OUTPUT_SITEYEAR_CSV, index=False)
print(f"Predictor table -> '{OUTPUT_SITEYEAR_CSV}' ({len(analysis_df)} rows).")

# Exploratory pooled correlations (autumn parameter x predictor) - a first
# look before the site-controlled models in scripts 08/09.
FLUX_PREDICTORS = ['total_gpp_growing_season', 'gpp_sos10_to_solstice', 'gpp_solstice_to_eos10',
                    'gpp_solstice_to_eos90', 'total_reco_growing_season', 'mean_radiation_growing_season',
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
