"""
PIPELINE STEP 23 - Fit a double-logistic phenology curve (Zhang/Beck/TIMESAT)
to each site-year x VI-index (NDVI, NIRv) series, and extract leaf_out
(green-up) and EOS (senescence) percentile-crossing DOYs plus the kinetic
(steepness) parameters for both phases.

The fitting itself (dormant-season background and snow fill, climatology
gap-fill, full-year transition dates, QC flags) lives in pheno_fit.py and is
shared with the tower (step 24) and PhenoCam (step 25) phenology - see that
module's docstring for the method and its thresholds.

Satellite-specific handling here:
  - observations with valid_pixel_frac < MIN_VALID_FRAC are dropped, and the
    rest are weighted by valid_pixel_frac
  - image dates with snow_frac >= SNOW_FRAC_MIN (from step 12) count as
    snow-covered and are filled with the dormant background

Output: data/phenology_double_logistic_by_site_year_index.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pheno_fit as pf

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
INPUT_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"
OUTPUT_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"

VI_COLUMNS = ['NDVI', 'NIRv']
PHYSICAL_BOUNDS = {'NDVI': (-1, 1), 'NIRv': (-1, 1)}
MIN_VALID_FRAC = 0.5
SNOW_FRAC_MIN = 0.5

if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run 22_merge_fluxnet_hls.py first.")

df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
df = df.dropna(subset=['date']).sort_values(['site_id', 'date']).copy()
df['year'] = df['date'].dt.year
df['doy'] = df['date'].dt.dayofyear

# frozen days from the full daily flux record (not only image dates)
df = pf.add_frozen_flag(df)

for vi_col, (lo, hi) in PHYSICAL_BOUNDS.items():
    if vi_col in df.columns:
        df.loc[(df[vi_col] < lo) | (df[vi_col] > hi), vi_col] = np.nan
if 'valid_pixel_frac' in df.columns:
    low_q = df['valid_pixel_frac'] < MIN_VALID_FRAC
    df.loc[low_q, [c for c in VI_COLUMNS if c in df.columns]] = np.nan
    df['fit_weight'] = df['valid_pixel_frac'].fillna(1.0).clip(lower=0.05, upper=1.0)
else:
    df['fit_weight'] = 1.0
has_snow_info = 'snow_frac' in df.columns
df['snowy'] = (df['snow_frac'] >= SNOW_FRAC_MIN) if has_snow_info else False
if not has_snow_info:
    print("NOTE: no snow_frac column (old step-12 output) - winter fill uses frozen days only.")

available_vi = [c for c in VI_COLUMNS if c in df.columns]
print(f"Fitting double-logistic phenology curves for: {available_vi}")

results = []
sites = df.groupby('site_id')
for si, (site_id, sdf) in enumerate(sites, start=1):
    for vi_col in available_vi:
        results.extend(pf.fit_site_index(site_id, sdf, vi_col))
    if si % 10 == 0 or si == sites.ngroups:
        print(f"  processed {si}/{sites.ngroups} sites")

pf.write_results(results, OUTPUT_CSV)
