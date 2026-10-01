"""
PIPELINE STEP 1 - Download and extract daily FLUXNET data.

Queries the global FLUXNET site catalog, keeps sites > 30 deg N with a
natural (non-wetland/urban/cropland) land cover, downloads each site's daily
archive, and keeps only the longest gap-free >= 10-year run per site.
Writes a single combined CSV with GPP/NEE/RECO/LE/TA/SW_IN/VPD/P and their
QC flags.

Output: data/fluxnet_daily_all_vars.csv
"""
import os
import re
import sys
import glob
import zipfile
from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR_RAW = Path(__file__).resolve().parent.parent / "data_raw"
DATA_DIR.mkdir(parents=True, exist_ok=True)

try:
    from fluxnet_shuttle import listall, download
except ImportError:
    import subprocess
    print("Installing 'fluxnet-shuttle'...")
    subprocess.check_call([sys.executable, "-m", "pip", "install",
                            "git+https://github.com/fluxnet/shuttle.git"])
    from fluxnet_shuttle import listall, download

TARGET_MAPPING = {
    'TIMESTAMP': ['timestamp', 'date', 'time', 'date_str'],
    'GPP_NT_VUT_REF': ['gpp_nt_vut_ref', 'gpp_nt', 'gpp_nt_vut', 'gpp'],
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
    'VPD_F_QC': ['vpd_f_qc'],
    'P_F': ['p_f', 'precip', 'precipitation', 'p'],
    'P_F_QC': ['p_f_qc', 'precip_qc'],
}
QC_PAIRS = {
    'NEE_VUT_REF': 'NEE_VUT_REF_QC',
    'GPP_NT_VUT_REF': 'NEE_VUT_REF_QC',
    'RECO_NT_VUT_REF': 'NEE_VUT_REF_QC',
    'LE_F_MDS': 'LE_F_MDS_QC', 'TA_F': 'TA_F_QC',
    'SW_IN_F': 'SW_IN_F_QC', 'VPD_F': 'VPD_F_QC', 'P_F': 'P_F_QC',
}
QC_THRESHOLD = 2
FINAL_COLUMNS = [
    'site_id', 'lat', 'lon', 'igbp', 'TIMESTAMP',
    'GPP_NT_VUT_REF', 'GPP_NT_VUT_REF_QC', 'NEE_VUT_REF', 'NEE_VUT_REF_QC',
    'RECO_NT_VUT_REF', 'RECO_NT_VUT_REF_QC', 'LE_F_MDS', 'LE_F_MDS_QC',
    'TA_F', 'TA_F_QC', 'SW_IN_F', 'SW_IN_F_QC', 'VPD_F', 'VPD_F_QC', 'P_F', 'P_F_QC',
]

MIN_LAT = 30.0
MIN_CONSECUTIVE_YEARS = 10
MIN_CONSECUTIVE_DAYS = round(365.25 * MIN_CONSECUTIVE_YEARS)
EXCLUDED_LC = {'WET', 'URB', 'CRO', 'CVM'}
DOWNLOAD_DIR = DATA_DIR_RAW
OUTPUT_CSV = DATA_DIR / "fluxnet_daily_all_vars.csv"
os.makedirs(DOWNLOAD_DIR, exist_ok=True)


def classify_land_cover(igbp_raw):
    if igbp_raw is None or pd.isna(igbp_raw):
        return None
    s = str(igbp_raw).strip().upper()
    aliases = {
        'CROPLAND': 'CRO', 'CROPLANDS': 'CRO',
        'CROPLAND/NATURAL VEGETATION MOSAIC': 'CVM', 'CROPLAND/NATURAL VEGETATION MOSAICS': 'CVM',
        'PERMANENT WETLANDS': 'WET', 'WETLAND': 'WET', 'WETLANDS': 'WET',
        'EVERGREEN NEEDLELEAF FOREST': 'ENF', 'EVERGREEN NEEDLELEAF FORESTS': 'ENF',
        'EVERGREEN BROADLEAF FOREST': 'EBF', 'EVERGREEN BROADLEAF FORESTS': 'EBF',
        'GRASSLAND': 'GRA', 'GRASSLANDS': 'GRA',
        'URBAN': 'URB', 'URBAN AND BUILT-UP': 'URB', 'URBAN AND BUILT-UP LANDS': 'URB',
        'BUILT-UP': 'URB', 'BUILT-UP LANDS': 'URB',
    }
    return aliases.get(s, s[:3])


def parse_year_range_from_filename(filename):
    match = re.search(r'(19|20)\d{2}-(19|20)\d{2}', filename)
    if not match:
        return None, None
    start_str, end_str = match.group(0).split('-')
    return int(start_str), int(end_str)


def longest_consecutive_window(sorted_unique_dates):
    if len(sorted_unique_dates) == 0:
        return None, None, 0
    run_start = sorted_unique_dates[0]
    best_start, best_end, best_len = run_start, run_start, 1
    prev, cur_len = run_start, 1
    for d in sorted_unique_dates[1:]:
        if (d - prev).days == 1:
            cur_len += 1
        else:
            if cur_len > best_len:
                best_start, best_end, best_len = run_start, prev, cur_len
            run_start, cur_len = d, 1
        prev = d
    if cur_len > best_len:
        best_start, best_end, best_len = run_start, prev, cur_len
    return best_start, best_end, best_len


print("Generating global FLUXNET site snapshot...")
latest_snapshot = listall(output_dir=str(DATA_DIR))
print(f"Snapshot written to '{latest_snapshot}'.", flush=True)
df_catalog = pd.read_csv(latest_snapshot)
cols = {str(c).lower().strip(): c for c in df_catalog.columns}
lat_col = next(cols[k] for k in cols if k in ['lat', 'latitude', 'location_lat'])
lon_col = next(cols[k] for k in cols if k in ['lon', 'longitude', 'location_long'])
site_col = next(cols[k] for k in cols if k in ['site_id', 'site', 'id'])
igbp_col = next((cols[k] for k in cols if k in
                  ['igbp', 'igbp_land_cover', 'igbp_class', 'vegetation_type', 'land_cover']), None)

df_catalog['lat_num'] = pd.to_numeric(df_catalog[lat_col], errors='coerce')
df_catalog['lon_num'] = pd.to_numeric(df_catalog[lon_col], errors='coerce')
candidate = df_catalog[(df_catalog['lat_num'] > MIN_LAT) & df_catalog['lat_num'].notnull()] \
    .drop_duplicates(subset=[site_col]).reset_index(drop=True)

if igbp_col is not None:
    lc_codes = candidate[igbp_col].apply(classify_land_cover)
    target_sites_df = candidate[~lc_codes.isin(EXCLUDED_LC)].reset_index(drop=True)
else:
    target_sites_df = candidate

target_site_ids = target_sites_df[site_col].tolist()
print(f"Found {len(target_site_ids)} candidate sites (Lat > {MIN_LAT}N, natural land cover).", flush=True)

site_meta = {
    str(row[site_col]): {'lat': float(row['lat_num']), 'lon': float(row['lon_num']),
                          'igbp': row[igbp_col] if igbp_col and pd.notna(row[igbp_col]) else None}
    for _, row in target_sites_df.iterrows()
}

existing_files = glob.glob(os.path.join(DOWNLOAD_DIR, "*.zip")) + glob.glob(os.path.join(DOWNLOAD_DIR, "*.csv"))
existing_filenames = [os.path.basename(f) for f in existing_files]


def _already_downloaded(sid, filenames):
    pattern = re.compile(rf'(^|[_-]){re.escape(sid)}([_-]|$)', re.IGNORECASE)
    return any(pattern.search(fname) for fname in filenames)


sites_to_download = [sid for sid in target_site_ids if not _already_downloaded(sid, existing_filenames)]
print(f"{len(existing_filenames)} file(s) already present in '{DOWNLOAD_DIR}'.", flush=True)
print(f"{len(target_site_ids) - len(sites_to_download)} already downloaded, {len(sites_to_download)} to fetch.",
      flush=True)
if sites_to_download:
    print(f"Calling download() for {len(sites_to_download)} site(s) - this step talks to the network "
          "and prints its own progress, which may take a while with no output here in between.", flush=True)
    download(site_ids=sites_to_download, snapshot_file=latest_snapshot, output_dir=str(DOWNLOAD_DIR))
else:
    print("Nothing to download - every candidate site already has a local file. Moving on to extraction.",
          flush=True)

print(f"\nExtracting, requiring >= {MIN_CONSECUTIVE_YEARS}y consecutive daily record...")
if os.path.exists(OUTPUT_CSV):
    os.remove(OUTPUT_CSV)

zip_files = glob.glob(os.path.join(DOWNLOAD_DIR, "*.zip")) + glob.glob(os.path.join(DOWNLOAD_DIR, "*.csv"))
processed_count = 0

for filepath in zip_files:
    try:
        filename = os.path.basename(filepath)
        site_id = next((sid for sid in target_site_ids
                         if re.search(rf'(^|[_-]){re.escape(sid)}([_-]|$)', filename, re.IGNORECASE)), None)
        if not site_id:
            continue

        span_start, span_end = parse_year_range_from_filename(filename)
        if span_start is not None and (span_end - span_start + 1) < MIN_CONSECUTIVE_YEARS:
            continue

        meta = site_meta.get(site_id, {})
        lat, lon, igbp = meta.get('lat'), meta.get('lon'), meta.get('igbp')
        lc_code = classify_land_cover(igbp)
        df_raw = None

        if filepath.endswith('.zip'):
            with zipfile.ZipFile(filepath, 'r') as z:
                dd_files = [f for f in z.namelist() if ('_DD_' in f or '_DD.' in f or 'daily' in f.lower())
                            and f.endswith('.csv')]
                if dd_files:
                    fullset = [f for f in dd_files if 'FULLSET' in f.upper()]
                    chosen = fullset[0] if fullset else dd_files[0]
                    with z.open(chosen) as csv_file:
                        df_raw = pd.read_csv(csv_file)
        elif filepath.endswith('.csv'):
            df_raw = pd.read_csv(filepath)

        if df_raw is None or df_raw.empty:
            continue

        existing_cols = {str(c).lower().strip(): c for c in df_raw.columns}
        ts_col_name = next((existing_cols[a] for a in TARGET_MAPPING['TIMESTAMP'] if a in existing_cols), None)
        if not ts_col_name:
            continue

        ts_dates = pd.to_datetime(df_raw[ts_col_name].astype(str), format='%Y%m%d', errors='coerce')
        if ts_dates.isna().all():
            ts_dates = pd.to_datetime(df_raw[ts_col_name].astype(str), errors='coerce')
        valid_mask = ts_dates.notna()
        if not valid_mask.any():
            continue
        df_valid, ts_valid = df_raw.loc[valid_mask].copy(), ts_dates.loc[valid_mask]

        unique_sorted_dates = sorted(pd.Timestamp(d) for d in ts_valid.dropna().unique())
        run_start, run_end, run_len = longest_consecutive_window(unique_sorted_dates)
        if run_len < MIN_CONSECUTIVE_DAYS or lc_code in EXCLUDED_LC:
            continue

        window_mask = (ts_valid >= run_start) & (ts_valid <= run_end)
        df_window, ts_window = df_valid.loc[window_mask], ts_valid.loc[window_mask]

        clean_df = pd.DataFrame()
        clean_df['site_id'] = [site_id] * len(df_window)
        clean_df['lat'] = [lat] * len(df_window)
        clean_df['lon'] = [lon] * len(df_window)
        clean_df['igbp'] = [igbp] * len(df_window)

        for target_var, aliases in TARGET_MAPPING.items():
            if target_var == 'TIMESTAMP':
                clean_df['TIMESTAMP'] = df_window[ts_col_name].to_numpy()
                continue
            found_col = next((existing_cols[a] for a in aliases if a in existing_cols), None)
            clean_df[target_var] = df_window[found_col].to_numpy() if found_col else None

        if 'NEE_VUT_REF_QC' in clean_df.columns:
            clean_df['GPP_NT_VUT_REF_QC'] = clean_df['NEE_VUT_REF_QC']
            clean_df['RECO_NT_VUT_REF_QC'] = clean_df['NEE_VUT_REF_QC']
        else:
            clean_df['GPP_NT_VUT_REF_QC'] = None
            clean_df['RECO_NT_VUT_REF_QC'] = None

        for var, qc_col in QC_PAIRS.items():
            qc_source = 'NEE_VUT_REF_QC' if var in ('GPP_NT_VUT_REF', 'RECO_NT_VUT_REF') else qc_col
            if qc_source in clean_df.columns:
                qc_numeric = pd.to_numeric(clean_df[qc_source], errors='coerce')
                clean_df.loc[qc_numeric > QC_THRESHOLD, var] = None

        clean_df = clean_df[FINAL_COLUMNS]
        file_exists = os.path.exists(OUTPUT_CSV)
        clean_df.to_csv(OUTPUT_CSV, mode='a', header=not file_exists, index=False)
        processed_count += 1
        print(f"   -> PASSED {site_id} ({run_start.date()} to {run_end.date()}, {run_len}d): "
              f"{len(clean_df):,} records")

    except Exception as e:
        print(f"   -> Skipped {filepath}: {e}")

print(f"\nDone. {processed_count} sites saved to '{OUTPUT_CSV}'.")