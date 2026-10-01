"""
PIPELINE STEP 7 - Rate-of-change (slope, units/day) covariates over the
senescence window (EOS90 -> EOS10, onset-of-decline to near-dormancy):

    temperature_senescence_rate  - linear-regression slope of TA_F vs DOY
    photoperiod_senescence_rate  - exact two-point slope of day length vs DOY

These are the strongest independent predictors of autumn timing found in
this pipeline's modeling (scripts 08/09), so they're kept as covariates
there. (gpp_greenup_rate, a spring-side rate covariate from earlier
exploration, was dropped - it never showed a significant effect and isn't
part of the opposite-effect hypothesis.)

Output: data/rate_of_change_predictors_by_site_year_index.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # step 5 output
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                            # step 4 output
OUTPUT_CSV = DATA_DIR / "rate_of_change_predictors_by_site_year_index.csv"

VI_INDICES = ['NDVI', 'NIRv']
MIN_FIT_CORR = 0.8
MIN_POINTS_FOR_RATE = 5
TEMP_COL, TEMP_RANGE = 'TA_F', (-60, 50)

if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run 05_double_logistic_phenology.py first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run 04_merge_fluxnet_landsat.py first.")

pheno = pd.read_csv(PHENOLOGY_CSV)
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')
              & (pheno['corr'] >= MIN_FIT_CORR)].copy()

flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year
flux['doy'] = flux['date'].dt.dayofyear

if TEMP_COL not in flux.columns:
    raise ValueError(f"'{TEMP_COL}' not found in FLUX_CSV.")
lo, hi = TEMP_RANGE
implausible = (flux[TEMP_COL] < lo) | (flux[TEMP_COL] > hi)
if implausible.any():
    flux.loc[implausible, TEMP_COL] = np.nan
temp_lookup = {key: group[['doy', TEMP_COL]].dropna(subset=['doy'])
               for key, group in flux.groupby(['site_id', 'year'])}
site_lat = flux.groupby('site_id')['lat'].first().to_dict()


def rate_temperature(site_id, year, doy_start, doy_end):
    if pd.isna(doy_start) or pd.isna(doy_end) or doy_start >= doy_end:
        return np.nan
    g = temp_lookup.get((site_id, year))
    if g is None or g.empty:
        return np.nan
    mask = (g['doy'] >= doy_start) & (g['doy'] <= doy_end)
    sub = g.loc[mask, ['doy', TEMP_COL]].dropna()
    if len(sub) < MIN_POINTS_FOR_RATE or sub['doy'].nunique() < 2:
        return np.nan
    slope, *_ = linregress(sub['doy'], sub[TEMP_COL])
    return float(slope)


def photoperiod_hours(doy, lat_deg):
    if pd.isna(doy) or pd.isna(lat_deg):
        return np.nan
    lat_rad = np.radians(lat_deg)
    declination_rad = np.radians(-23.44 * np.cos(np.radians(360.0 / 365.0 * (doy + 10))))
    cos_hour_angle = np.clip(-np.tan(lat_rad) * np.tan(declination_rad), -1.0, 1.0)
    return float((24.0 / np.pi) * np.arccos(cos_hour_angle))


def photoperiod_rate(doy_start, doy_end, lat_deg):
    if pd.isna(doy_start) or pd.isna(doy_end) or pd.isna(lat_deg) or doy_start >= doy_end:
        return np.nan
    p_start, p_end = photoperiod_hours(doy_start, lat_deg), photoperiod_hours(doy_end, lat_deg)
    if pd.isna(p_start) or pd.isna(p_end):
        return np.nan
    return (p_end - p_start) / (doy_end - doy_start)


records = []
for _, row in pheno.iterrows():
    site_id, year = row['site_id'], int(row['year'])
    eos90, eos10 = row.get('EOS90'), row.get('EOS10')
    lat = site_lat.get(site_id)
    records.append({
        'site_id': site_id, 'year': year, 'vi_index': row['vi_index'],
        'EOS90': eos90, 'EOS50': row.get('EOS50'), 'EOS10': eos10,
        'senescence_kinetic_i': row.get('senescence_kinetic_i'), 'leaf_out_90': row.get('leaf_out_90'),
        'temperature_senescence_rate': rate_temperature(site_id, year, eos90, eos10),
        'photoperiod_senescence_rate': photoperiod_rate(eos90, eos10, lat),
    })

analysis_df = pd.DataFrame(records)
analysis_df.to_csv(OUTPUT_CSV, index=False)
print(f"Rate-of-change predictors -> '{OUTPUT_CSV}' ({len(analysis_df)} rows).")
