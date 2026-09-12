import os
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

FLUX_CSV = DATA_DIR / "fluxnet_daily_all_vars_qc_filtered.csv"      # output of the growing-season QC filter
QC_SUMMARY_CSV = DATA_DIR / "site_year_growing_season_qc_summary.csv"  # per-site-year pass/fail, same script
LANDSAT_CSV = DATA_DIR / "fluxnet_all_highlat_landsat_indices.csv"   # output of the Landsat (HLSL30) script

OUTPUT_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"

# 'left'  - keep every (QC-passing) FLUXNET daily row, attach Landsat VI
#           columns where a clear-sky overpass exists for that site+date
#           (NaN elsewhere). Default: FLUXNET is the dense daily record,
#           Landsat is inherently sparse (16-day native revisit, further
#           thinned by cloud/shadow/snow masking), so most days simply
#           won't have a match and that's expected, not a bug.
# 'inner' - keep only site+date rows present in BOTH tables, i.e. only
#           actual flux-days that also have a usable Landsat observation.
#           Switch to this if what you need is matched pairs (e.g. directly
#           relating GPP to VI on the same day) rather than a continuous
#           daily series with VI as a sparse covariate.
MERGE_HOW = 'left'


# ---------------------------------------------------------------------------
# 1. Load + parse dates
# ---------------------------------------------------------------------------
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run the growing-season QC filter script first.")
if not os.path.exists(QC_SUMMARY_CSV):
    raise FileNotFoundError(f"Missing '{QC_SUMMARY_CSV}'. Run the growing-season QC filter script first.")
if not os.path.exists(LANDSAT_CSV):
    raise FileNotFoundError(f"Missing '{LANDSAT_CSV}'. Run the Landsat (HLSL30) extraction script first.")

flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year

landsat = pd.read_csv(LANDSAT_CSV)
landsat['date'] = pd.to_datetime(landsat['date'], errors='coerce')
landsat = landsat.dropna(subset=['date']).copy()
landsat['year'] = landsat['date'].dt.year

print(f"FLUXNET (QC-filtered): {len(flux):,} rows, {flux['site_id'].nunique()} sites.")
print(f"Landsat indices:       {len(landsat):,} rows, {landsat['site_id'].nunique()} sites.")

# ---------------------------------------------------------------------------
# 2. Restrict BOTH tables to QC-passing site-years, explicitly
# ---------------------------------------------------------------------------
# The Landsat download script pulls imagery for a site's full flux date
# range without knowing which of those years will later pass the
# growing-season QC filter - it can't, since that filter runs afterward.
# So the Landsat indices file can contain years that never passed QC (e.g.
# a year later dropped for >=50% low-quality or a >=50-day gap in GPP/NEE/
# RECO). FLUX_CSV is already QC-filtered so this is close to a no-op on the
# FLUXNET side, but the explicit filter is applied to both sides here so
# the restriction doesn't just happen implicitly as a side effect of which
# file FLUX_CSV happens to point to - if this script is ever pointed at an
# unfiltered FLUXNET file, the QC restriction still holds.
qc_summary = pd.read_csv(QC_SUMMARY_CSV)
passing = qc_summary.loc[qc_summary['passed'], ['site_id', 'year']].drop_duplicates()
print(f"QC-passing site-years: {len(passing):,} (out of {len(qc_summary):,} evaluated).")

n_flux_before = len(flux)
flux = flux.merge(passing, on=['site_id', 'year'], how='inner')
if len(flux) < n_flux_before:
    print(f"Dropped {n_flux_before - len(flux):,} FLUXNET row(s) belonging to non-passing site-years "
          "(should be ~0 if FLUX_CSV is already the QC-filtered output).")

n_landsat_before = len(landsat)
landsat = landsat.merge(passing, on=['site_id', 'year'], how='inner')
print(f"Dropped {n_landsat_before - len(landsat):,} Landsat row(s) belonging to non-passing site-years "
      f"({len(landsat):,} of {n_landsat_before:,} kept).")

# ---------------------------------------------------------------------------
# 3. Collapse Landsat to one row per site_id+date before merging
# ---------------------------------------------------------------------------
# The raw-band extraction step already averages same-day observations within
# a site (site_id, lat, lon, sensor, date), so duplicates at this point
# should be rare, but this is a defensive step in case the indices file was
# regenerated some other way (e.g. concatenated across multiple runs).
vi_cols = [c for c in ['NDVI', 'EVI', 'NIRv', 'valid_pixel_frac'] if c in landsat.columns]
n_before = len(landsat)
landsat_agg = (landsat
               .groupby(['site_id', 'date'], as_index=False)[vi_cols]
               .mean())
n_after = len(landsat_agg)
if n_after < n_before:
    print(f"Collapsed {n_before - n_after} duplicate site+date Landsat row(s) via mean.")

# ---------------------------------------------------------------------------
# 4. Merge on site_id + date
# ---------------------------------------------------------------------------
merged = flux.merge(landsat_agg, on=['site_id', 'date'], how=MERGE_HOW)

merged['has_landsat_obs'] = merged[vi_cols[0]].notna() if vi_cols else False

n_sites_flux = set(flux['site_id'].unique())
n_sites_landsat = set(landsat['site_id'].unique())
n_common = len(n_sites_flux & n_sites_landsat)
print(f"\nSites in both (QC-restricted) datasets: {n_common} "
      f"(FLUXNET-only: {len(n_sites_flux - n_sites_landsat)}, "
      f"Landsat-only: {len(n_sites_landsat - n_sites_flux)} - "
      "Landsat-only sites contribute no rows since FLUXNET drives the merge.)")

n_matched_days = int(merged['has_landsat_obs'].sum())
print(f"Merged table: {len(merged):,} rows ('{MERGE_HOW}' join on site_id+date), "
      f"{n_matched_days:,} rows ({100 * n_matched_days / max(len(merged), 1):.1f}%) "
      "have a matching Landsat observation.")

# TIMESTAMP already carries the date in the FLUXNET format (YYYYMMDD);
# drop the redundant datetime/year helper columns before saving.
merged = merged.drop(columns=['date', 'year'])
merged.to_csv(OUTPUT_CSV, index=False)
print(f"\nMerged FLUXNET + Landsat table (QC-passing site-years only) written to '{OUTPUT_CSV}'.")
