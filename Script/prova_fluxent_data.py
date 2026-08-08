import os
import re
import sys
import glob
import zipfile
from pathlib import Path
import pandas as pd

# ---------------------------------------------------------------------------
# 0. Repo-relative paths
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
DATA_DIR_RAW = Path(__file__).resolve().parent.parent / "Data_raw"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Import fluxnet_shuttle
# ---------------------------------------------------------------------------
try:
    from fluxnet_shuttle import listall, download
except ImportError:
    import subprocess
    print("Installing official 'fluxnet-shuttle' library from GitHub...")
    subprocess.check_call([
        sys.executable, "-m", "pip", "install",
        "git+https://github.com/fluxnet/shuttle.git"
    ])
    from fluxnet_shuttle import listall, download

# ---------------------------------------------------------------------------
# 2. Configuration & Variable Mapping
# ---------------------------------------------------------------------------
TARGET_MAPPING = {
    'TIMESTAMP': ['timestamp', 'date', 'time', 'date_str'],
    'GPP_NT_VUT_REF': ['gpp_nt_vut_ref', 'gpp_nt', 'gpp_nt_vut', 'gpp'],
    'GPP_DT_VUT_REF': ['gpp_dt_vut_ref', 'gpp_dt', 'gpp_dt_vut'],
    'NEE_VUT_REF': ['nee_vut_ref', 'nee_vut', 'nee'],
    'NEE_VUT_REF_QC': ['nee_vut_ref_qc', 'nee_vut_qc', 'nee_qc'],
    'RECO_NT_VUT_REF': ['reco_nt_vut_ref', 'reco_nt', 'reco_nt_vut', 'reco'],
    'LE_F_MDS': ['le_f_mds', 'le_f', 'le'],
    'LE_F_MDS_QC': ['le_f_mds_qc', 'le_f_qc'],
    'TA_F': ['ta_f', 'ta', 'temp'],
    'TA_F_QC': ['ta_f_qc'],
    'SW_IN_F': ['sw_in_f', 'sw_in', 'swin'],
    'SW_IN_F_QC': ['sw_in_f_qc'],
    'VPD_F': ['vpd_f', 'vpd'],
    'VPD_F_QC': ['vpd_f_qc']
}

# GPP_NT/DT and RECO are model outputs from the NEE partitioning step and don't
# carry their own QC flag in the FLUXNET daily product - they inherit NEE's QC.
QC_PAIRS = {
    'NEE_VUT_REF': 'NEE_VUT_REF_QC',
    'LE_F_MDS': 'LE_F_MDS_QC',
    'TA_F': 'TA_F_QC',
    'SW_IN_F': 'SW_IN_F_QC',
    'VPD_F': 'VPD_F_QC',
}
# FLUXNET QC convention: 0=measured, 1=good gap-fill, 2=medium, 3=poor.
# Values with QC > QC_THRESHOLD are set to NaN rather than dropped, so the
# QC columns stay available downstream if you want a looser/stricter cut later.
QC_THRESHOLD = 2

FINAL_COLUMNS = [
    'site_id', 'lat', 'lon', 'igbp', 'TIMESTAMP',
    'GPP_NT_VUT_REF', 'GPP_DT_VUT_REF', 'NEE_VUT_REF', 'NEE_VUT_REF_QC',
    'RECO_NT_VUT_REF', 'LE_F_MDS', 'LE_F_MDS_QC',
    'TA_F', 'TA_F_QC', 'SW_IN_F', 'SW_IN_F_QC', 'VPD_F', 'VPD_F_QC'
]

# "High latitude" cutoff - 30N is really just "Northern Hemisphere extratropical".
# Bump this up (e.g. 50-55) if you want to restrict to boreal/subarctic sites.
MIN_LAT = 30.0

# A site is kept only if its daily record extends to at least CUTOFF_YEAR.
# This used to be 2018 because Sentinel-2 (the only optical source in the
# pipeline at the time) starts in 2015/2017 and needed a couple of years of
# overlap to be worth including a site for. Now that HLS is in the pipeline
# (HLSL30 draws on Landsat 8, which starts Feb 2013), that constraint is
# gone - a site whose record extends to 2013 already has usable HLS overlap,
# so the cutoff is loosened to match. Note this only affects which sites
# make it into fluxnet_daily_selected_vars.csv; it does not change how far
# back the daily record for each kept site extends (that's the full FLUXNET
# record either way).
CUTOFF_YEAR = 2013

# All downloaded archives and derived outputs live under data/ (processed data).
DOWNLOAD_DIR = DATA_DIR_RAW 
OUTPUT_CSV = DATA_DIR / "fluxnet_daily_selected_vars.csv"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 3. Query Master Catalog & Filter Sites > 30°N
# ---------------------------------------------------------------------------
print("Generating global FLUXNET site snapshot...")
latest_snapshot = listall(output_dir=str(DATA_DIR))
print(f"Catalog generated: {latest_snapshot}")

df_catalog = pd.read_csv(latest_snapshot)

cols = {str(c).lower().strip(): c for c in df_catalog.columns}
lat_col = next(cols[k] for k in cols if k in ['lat', 'latitude', 'location_lat'])
lon_col = next(cols[k] for k in cols if k in ['lon', 'longitude', 'location_long'])
site_col = next(cols[k] for k in cols if k in ['site_id', 'site', 'id'])
igbp_col = next(
    (cols[k] for k in cols if k in ['igbp', 'igbp_land_cover', 'igbp_class', 'vegetation_type', 'land_cover']),
    None
)
if igbp_col is None:
    print("WARNING: No IGBP/land-cover column found in the site catalog - "
          "'igbp' will be saved as blank. Check df_catalog.columns and add the "
          "right alias to the igbp_col lookup above if the metadata uses a different name.")

df_catalog['lat_num'] = pd.to_numeric(df_catalog[lat_col], errors='coerce')
df_catalog['lon_num'] = pd.to_numeric(df_catalog[lon_col], errors='coerce')

target_sites_df = df_catalog[
    (df_catalog['lat_num'] > MIN_LAT) &
    df_catalog['lat_num'].notnull()
].drop_duplicates(subset=[site_col]).reset_index(drop=True)

target_site_ids = target_sites_df[site_col].tolist()
print(f"Found {len(target_site_ids)} global FLUXNET sites with Latitude > {MIN_LAT}°N.")

site_meta = {
    str(row[site_col]): {
        'lat': float(row['lat_num']),
        'lon': float(row['lon_num']),
        'igbp': row[igbp_col] if igbp_col and pd.notna(row[igbp_col]) else None,
    }
    for _, row in target_sites_df.iterrows()
}

# ---------------------------------------------------------------------------
# 4. Batch Download Daily Flux Archives
# ---------------------------------------------------------------------------
print("\nChecking for already-downloaded sites...")

existing_files = glob.glob(os.path.join(DOWNLOAD_DIR, "*.zip")) + glob.glob(os.path.join(DOWNLOAD_DIR, "*.csv"))
existing_filenames = [os.path.basename(f) for f in existing_files]


def _site_already_downloaded(sid, filenames):
    pattern = re.compile(rf'(^|[_-]){re.escape(sid)}([_-]|$)', re.IGNORECASE)
    return any(pattern.search(fname) for fname in filenames)


already_downloaded = [sid for sid in target_site_ids if _site_already_downloaded(sid, existing_filenames)]
sites_to_download = [sid for sid in target_site_ids if sid not in already_downloaded]

print(f"  {len(already_downloaded)} site(s) already present in '{DOWNLOAD_DIR}' - skipping re-download.")
print(f"  {len(sites_to_download)} site(s) remaining to download.")

if sites_to_download:
    download(
        site_ids=sites_to_download,
        snapshot_file=latest_snapshot,
        output_dir=str(DOWNLOAD_DIR)
    )
else:
    print("All target sites already downloaded - skipping download step entirely.")

# ---------------------------------------------------------------------------
# 5. Extract, Filter CUTOFF_YEAR+ & Target Variables
# ---------------------------------------------------------------------------
print(f"\nExtracting, filtering {CUTOFF_YEAR}+ cutoff, and processing target variables...")

if os.path.exists(OUTPUT_CSV):
    os.remove(OUTPUT_CSV)

zip_files = glob.glob(os.path.join(DOWNLOAD_DIR, "*.zip")) + glob.glob(os.path.join(DOWNLOAD_DIR, "*.csv"))
processed_count = 0

for filepath in zip_files:
    try:
        filename = os.path.basename(filepath)
        site_id = None
        for sid in target_site_ids:
            if re.search(rf'(^|[_-]){re.escape(sid)}([_-]|$)', filename, re.IGNORECASE):
                site_id = sid
                break

        if not site_id:
            continue

        meta = site_meta.get(site_id, {})
        lat, lon, igbp = meta.get('lat'), meta.get('lon'), meta.get('igbp')
        df_raw = None

        if filepath.endswith('.zip'):
            with zipfile.ZipFile(filepath, 'r') as z:
                dd_files = [f for f in z.namelist() if ('_DD_' in f or '_DD.' in f or 'daily' in f.lower()) and f.endswith('.csv')]
                if dd_files:
                    # Prefer FULLSET over SUBSET if a zip happens to contain both
                    fullset_files = [f for f in dd_files if 'FULLSET' in f.upper()]
                    chosen_file = fullset_files[0] if fullset_files else dd_files[0]
                    if len(dd_files) > 1:
                        print(f"   -> Note: {site_id} zip has {len(dd_files)} daily-file candidates, using '{chosen_file}'")
                    with z.open(chosen_file) as csv_file:
                        df_raw = pd.read_csv(csv_file)
        elif filepath.endswith('.csv'):
            df_raw = pd.read_csv(filepath)

        if df_raw is not None and not df_raw.empty:
            existing_cols = {str(c).lower().strip(): c for c in df_raw.columns}

            # Find TIMESTAMP column
            ts_col_name = None
            for alias in TARGET_MAPPING['TIMESTAMP']:
                if alias in existing_cols:
                    ts_col_name = existing_cols[alias]
                    break

            if not ts_col_name:
                print(f"   -> Skipped {site_id}: missing timestamp column")
                continue

            # Check if record reaches at least CUTOFF_YEAR
            ts_dates = pd.to_datetime(df_raw[ts_col_name].astype(str), format='%Y%m%d', errors='coerce')
            if ts_dates.isna().all():
                ts_dates = pd.to_datetime(df_raw[ts_col_name].astype(str), errors='coerce')

            max_year = ts_dates.dt.year.max()
            if pd.isna(max_year) or max_year < CUTOFF_YEAR:
                print(f"   -> Skipped {site_id}: Record ends in {int(max_year) if pd.notna(max_year) else 'N/A'} (< {CUTOFF_YEAR})")
                continue

            clean_df = pd.DataFrame()
            clean_df['site_id'] = [site_id] * len(df_raw)
            clean_df['lat'] = [lat] * len(df_raw)
            clean_df['lon'] = [lon] * len(df_raw)
            clean_df['igbp'] = [igbp] * len(df_raw)

            for target_var, aliases in TARGET_MAPPING.items():
                found_col = None
                for alias in aliases:
                    if alias in existing_cols:
                        found_col = existing_cols[alias]
                        break
                clean_df[target_var] = df_raw[found_col] if found_col else None

            # Blank out poor-quality (heavily gap-filled) values rather than dropping rows,
            # so downstream z-scores/cumulative sums aren't biased by low-confidence data.
            for var, qc_col in QC_PAIRS.items():
                if qc_col in clean_df.columns:
                    qc_numeric = pd.to_numeric(clean_df[qc_col], errors='coerce')
                    clean_df.loc[qc_numeric > QC_THRESHOLD, var] = None

            clean_df = clean_df[FINAL_COLUMNS]

            file_exists = os.path.exists(OUTPUT_CSV)
            clean_df.to_csv(OUTPUT_CSV, mode='a', header=not file_exists, index=False)
            processed_count += 1
            print(f"   -> PASSED {site_id} (Max year: {int(max_year)}): {len(clean_df):,} daily records saved")

    except Exception as e:
        print(f"   -> Skipped {filepath}: {e}")

print("\n" + "=" * 60)
print(f"Extraction complete! Saved {processed_count} sites with records extending to at least {CUTOFF_YEAR} to '{OUTPUT_CSV}'.")