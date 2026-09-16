"""
PIPELINE STEP 4 - Merge QC-passing FLUXNET daily data with sparse Landsat VI
observations on (site_id, date). Both tables are first restricted to
QC-passing site-years explicitly.

Output: data/fluxnet_landsat_merged.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FLUX_CSV = DATA_DIR / "fluxnet_daily_all_vars_qc_filtered.csv"
QC_SUMMARY_CSV = DATA_DIR / "site_year_growing_season_qc_summary.csv"
LANDSAT_CSV = DATA_DIR / "fluxnet_all_highlat_landsat_indices.csv"
OUTPUT_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"
MERGE_HOW = 'left'

for p in (FLUX_CSV, QC_SUMMARY_CSV, LANDSAT_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run steps 1-3 first.")

flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year

landsat = pd.read_csv(LANDSAT_CSV)
landsat['date'] = pd.to_datetime(landsat['date'], errors='coerce')
landsat = landsat.dropna(subset=['date']).copy()
landsat['year'] = landsat['date'].dt.year

qc_summary = pd.read_csv(QC_SUMMARY_CSV)
passing = qc_summary.loc[qc_summary['passed'], ['site_id', 'year']].drop_duplicates()
flux = flux.merge(passing, on=['site_id', 'year'], how='inner')
landsat = landsat.merge(passing, on=['site_id', 'year'], how='inner')

vi_cols = [c for c in ['NDVI', 'EVI', 'NIRv', 'valid_pixel_frac'] if c in landsat.columns]
landsat_agg = landsat.groupby(['site_id', 'date'], as_index=False)[vi_cols].mean()

merged = flux.merge(landsat_agg, on=['site_id', 'date'], how=MERGE_HOW)
merged['has_landsat_obs'] = merged[vi_cols[0]].notna() if vi_cols else False

n_matched_days = int(merged['has_landsat_obs'].sum())
print(f"Merged table: {len(merged):,} rows, {n_matched_days:,} "
      f"({100 * n_matched_days / max(len(merged), 1):.1f}%) have a Landsat match.")

merged = merged.drop(columns=['date', 'year'])
merged.to_csv(OUTPUT_CSV, index=False)
print(f"Written -> '{OUTPUT_CSV}'.")
