"""
PIPELINE STEP 25 - Phenology (leaf-out / EOS) from PhenoCam canopy greenness
(GCC), for the flux sites matched and downloaded in step 13.

Uses gcc_90 (90th percentile of daytime GCC, the PhenoCam standard),
dropping days flagged as outliers, and fits it with the shared routine in
pheno_fit.py - NOT PhenoCam's own transition dates, which exist only for
10/25/50% of amplitude and use a different curve; fitting here gives
EOS90/EOS50/EOS10 defined exactly as for the satellite (step 23) and tower
(step 24) series.

Only QC-passing flux site-years are fitted (those are the ones with GPP
predictors), and the flux record supplies TA_F for the frozen-day flag.
Days with PhenoCam snow_flag set count as snow-covered.

Input : data/phenocam_site_matches.csv, data_raw/phenocam/<roi>_1day.csv (step 13);
        data/fluxnet_landsat_merged.csv (step 22)
Output: data/phenocam_gcc_daily.csv
        data/phenology_phenocam_by_site_year_index.csv   (same columns as step 23)
"""
import os
from pathlib import Path

import pandas as pd
import pheno_fit as pf

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
CACHE_DIR = ROOT / "data_raw" / "phenocam"
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"
MATCH_CSV = DATA_DIR / "phenocam_site_matches.csv"
DAILY_CSV = DATA_DIR / "phenocam_gcc_daily.csv"
OUTPUT_CSV = DATA_DIR / "phenology_phenocam_by_site_year_index.csv"
GCC_COL = 'gcc_90'

if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run steps 11-22 first.")
if not os.path.exists(MATCH_CSV):
    raise FileNotFoundError(f"Missing '{MATCH_CSV}'. Run 13_phenocam_download.py first.")

flux = pd.read_csv(FLUX_CSV, usecols=['site_id', 'TIMESTAMP', 'TA_F'])
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux = flux.dropna(subset=['date']).sort_values(['site_id', 'date'])
flux = pf.add_frozen_flag(flux)
matches = pd.read_csv(MATCH_CSV)
matches = matches[matches['site_id'].isin(flux['site_id'].unique())].reset_index(drop=True)
print(f"PhenoCam GCC for {len(matches)} flux sites with QC-passing flux years...")

daily_parts, rows = [], []
for i, m in matches.iterrows():
    cache = CACHE_DIR / f"{m['roi_name']}_1day.csv"
    if not cache.exists():
        print(f"  {m['site_id']}: '{cache.name}' not downloaded - skipped (rerun step 13)")
        continue
    g = pd.read_csv(cache, comment='#')
    g['date'] = pd.to_datetime(g['date'], errors='coerce')
    gcc = pd.to_numeric(g[GCC_COL], errors='coerce')
    flag_col = f'outlierflag_{GCC_COL}'
    if flag_col in g:
        gcc = gcc.where(pd.to_numeric(g[flag_col], errors='coerce').fillna(0) == 0)
    g['gcc'] = gcc
    g['snowy'] = pd.to_numeric(g['snow_flag'], errors='coerce').fillna(0) > 0 if 'snow_flag' in g else False

    sdf = flux[flux['site_id'] == m['site_id']].merge(g[['date', 'gcc', 'snowy']], on='date', how='left')
    sdf['snowy'] = sdf['snowy'].eq(True)
    sdf['year'], sdf['doy'] = sdf['date'].dt.year, sdf['date'].dt.dayofyear
    sdf['fit_weight'] = 1.0
    daily_parts.append(sdf.loc[sdf['gcc'].notna() | sdf['snowy'], ['site_id', 'date', 'gcc', 'snowy']]
                       .assign(roi_name=m['roi_name']))
    n_before = len(rows)
    rows.extend(pf.fit_site_index(m['site_id'], sdf, 'gcc', label='GCC'))
    print(f"  [{i + 1}/{len(matches)}] {m['site_id']} <- {m['roi_name']} ({m['dist_km']} km): "
          f"{int(sdf['gcc'].notna().sum())} GCC days in flux years, {len(rows) - n_before} site-year fits")

if not rows:
    raise SystemExit("PhenoCam records do not overlap the QC-passing flux years of any site.")
pd.concat(daily_parts, ignore_index=True).to_csv(DAILY_CSV, index=False)
print(f"\nDaily GCC -> '{DAILY_CSV}'")
pf.write_results(rows, OUTPUT_CSV)
