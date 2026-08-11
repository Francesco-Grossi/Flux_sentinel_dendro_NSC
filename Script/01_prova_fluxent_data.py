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
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR_RAW = Path(__file__).resolve().parent.parent / "data_raw"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Per-land-cover-class output folders. Every site that passes all filters
# still lands in the single combined OUTPUT_CSV below (as before); sites
# classified as evergreen forest or grassland ALSO get their own per-site
# CSV copied into one of these two folders, since those are the classes
# planned for separate analysis.
EVERGREEN_DIR = DATA_DIR / "evergreen"
GRASSLAND_DIR = DATA_DIR / "grassland"
EVERGREEN_DIR.mkdir(parents=True, exist_ok=True)
GRASSLAND_DIR.mkdir(parents=True, exist_ok=True)

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
    'VPD_F_QC': ['vpd_f_qc'],
    # Precipitation - previously missing (prova_lmm_models.py already had a
    # note flagging this: "add it to TARGET_MAPPING ... to include it here").
    'P_F': ['p_f', 'precip', 'precipitation', 'p'],
    'P_F_QC': ['p_f_qc', 'precip_qc'],
}

# GPP_NT/DT and RECO are model outputs from the NEE partitioning step and don't
# carry their own QC flag in the FLUXNET daily product - they inherit NEE's QC.
QC_PAIRS = {
    'NEE_VUT_REF': 'NEE_VUT_REF_QC',
    'LE_F_MDS': 'LE_F_MDS_QC',
    'TA_F': 'TA_F_QC',
    'SW_IN_F': 'SW_IN_F_QC',
    'VPD_F': 'VPD_F_QC',
    'P_F': 'P_F_QC',
}
# FLUXNET QC convention: 0=measured, 1=good gap-fill, 2=medium, 3=poor.
# Values with QC > QC_THRESHOLD are set to NaN rather than dropped, so the
# QC columns stay available downstream if you want a looser/stricter cut later.
QC_THRESHOLD = 2

FINAL_COLUMNS = [
    'site_id', 'lat', 'lon', 'igbp', 'TIMESTAMP',
    'GPP_NT_VUT_REF', 'GPP_DT_VUT_REF', 'NEE_VUT_REF', 'NEE_VUT_REF_QC',
    'RECO_NT_VUT_REF', 'LE_F_MDS', 'LE_F_MDS_QC',
    'TA_F', 'TA_F_QC', 'SW_IN_F', 'SW_IN_F_QC', 'VPD_F', 'VPD_F_QC',
    'P_F', 'P_F_QC',
]

# "High latitude" cutoff - 30N is really just "Northern Hemisphere extratropical".
# Bump this up (e.g. 50-55) if you want to restrict to boreal/subarctic sites.
MIN_LAT = 30.0

# Every site's daily record is cropped to CUTOFF_YEAR-01-01 onwards before
# anything else (the consecutive-run check, land-cover filter don't care,
# but the length check below does). HLSL30 draws on Landsat 8, which starts
# Feb 2013, so pre-2013 flux data has no HLS imagery to pair with anyway -
# no point carrying it through the pipeline.
CUTOFF_YEAR = 2013
CUTOFF_DATE = pd.Timestamp(f"{CUTOFF_YEAR}-01-01")

# A site is kept only if, within its post-CUTOFF_YEAR record, it has at
# least this many consecutive calendar days with a daily timestamp (i.e. no
# gap in the date sequence) - not just "some data reaching back 10 years
# with holes in between". 365.25*10 rounds to 3653 to allow for leap years.
MIN_CONSECUTIVE_YEARS = 10
MIN_CONSECUTIVE_DAYS = round(365.25 * MIN_CONSECUTIVE_YEARS)

# IGBP land-cover classes to exclude outright (croplands are managed/
# harvested, and wetlands have hydrology-driven phenology - both would
# confound the GPP-vs-senescence analysis this pipeline is built around).
# CVM ("cropland/natural vegetation mosaic") is excluded alongside pure
# cropland (CRO) for the same reason - it's still cropland-dominated.
EXCLUDED_LC = {'CRO', 'CVM', 'WET'}
EVERGREEN_LC = {'ENF', 'EBF'}  # evergreen needleleaf / broadleaf forest
GRASSLAND_LC = {'GRA'}

# If the site catalog from listall() happens to expose each site's overall
# data-year range (some FLUXNET/AmeriFlux listings do, some don't), sites
# that obviously can't satisfy the 10y/2013+ requirement are dropped BEFORE
# download() is ever called on them - avoiding the bandwidth/disk cost
# entirely, rather than paying it and rejecting the site afterward. This is
# a coarse pre-filter only (it can't see gaps within the span, so a site
# that passes here can still get rejected later in Section 5 once its
# actual daily record is read) - it only ever removes sites that are
# already impossible on their catalog-listed span alone.
CATALOG_YEAR_COL_ALIASES = {
    'start_year': ['start_year', 'startyear', 'first_year', 'data_start_year', 'yr_start'],
    'end_year': ['end_year', 'endyear', 'last_year', 'data_end_year', 'yr_end'],
}

# Set True to delete a raw archive from data_raw/ once Section 5's fast
# filename-based pre-check confirms it can't possibly qualify - reclaims
# disk space, at the cost of re-downloading it if CUTOFF_YEAR or
# MIN_CONSECUTIVE_YEARS are ever loosened later. Off by default so this
# script never deletes something from disk without an explicit opt-in.
DELETE_REJECTED_RAW_FILES = False

# All downloaded archives and derived outputs live under data/ (processed data).
DOWNLOAD_DIR = DATA_DIR_RAW
OUTPUT_CSV = DATA_DIR / "fluxnet_daily_selected_vars.csv"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)


def find_year_range_column(df, sample_size=30, match_threshold=0.5):
    """Scan every column in the catalog for one whose values embed a
    year-range pattern like '1997-2010' - the same pattern already proven
    to exist in the downloaded archive filenames (see the screenshot that
    prompted this). Doesn't rely on knowing the exact column name (e.g.
    'start_year'/'end_year', which this library's public docs don't
    document) - if the snapshot has a download-URL or filename-preview
    column that already carries this pattern, it'll be found here too.
    Returns the column name, or None if nothing matches well enough."""
    for col in df.columns:
        sample = df[col].dropna().astype(str).head(sample_size)
        if sample.empty:
            continue
        hit_rate = sample.apply(lambda s: parse_year_range_from_filename(s)[0] is not None).mean()
        if hit_rate >= match_threshold:
            return col
    return None


def classify_land_cover(igbp_raw):
    """Normalize whatever the site catalog stores (3-letter IGBP code or a
    full name) down to the standard 3-letter IGBP code. Falls back to the
    first 3 letters uppercased for anything unrecognized, so an unexpected
    label doesn't silently vanish - it just won't match EXCLUDED_LC /
    EVERGREEN_LC / GRASSLAND_LC and the site is treated as "keep, unknown
    class" rather than being dropped."""
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
    }
    return aliases.get(s, s[:3])


def parse_year_range_from_filename(filename):
    """FLUXNET/AmeriFlux archive names embed the site's overall data-year
    range, e.g. 'AMF_CA-Ca1_FLUXNET_1997-2010_v1.3_r1.zip' -> (1997, 2010).
    Used to reject an archive BEFORE opening/unzipping it: if the archive's
    own overall span can't possibly contain a >= MIN_CONSECUTIVE_YEARS
    gap-free run reaching CUTOFF_YEAR+, no amount of parsing the daily data
    inside it will change that - so there's no reason to pay the unzip/read
    cost. This can only produce safe rejections (a gap-free run is always a
    subset of the overall span), never a false rejection of a usable site."""
    match = re.search(r'(19|20)\d{2}-(19|20)\d{2}', filename)
    if not match:
        return None, None
    start_str, end_str = match.group(0).split('-')
    return int(start_str), int(end_str)


def longest_consecutive_window(sorted_unique_dates):
    """Given a sorted array/list of unique daily Timestamps, return
    (start, end, length_days) of the longest run with no gap (every
    consecutive pair exactly 1 day apart)."""
    if len(sorted_unique_dates) == 0:
        return None, None, 0

    run_start = sorted_unique_dates[0]
    best_start, best_end, best_len = run_start, run_start, 1
    prev = run_start
    cur_len = 1

    for d in sorted_unique_dates[1:]:
        if (d - prev).days == 1:
            cur_len += 1
        else:
            if cur_len > best_len:
                best_start, best_end, best_len = run_start, prev, cur_len
            run_start = d
            cur_len = 1
        prev = d

    if cur_len > best_len:
        best_start, best_end, best_len = run_start, prev, cur_len

    return best_start, best_end, best_len


# ---------------------------------------------------------------------------
# 3. Query Master Catalog & Filter Sites > 30°N, drop cropland/wetland
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
          "cropland/wetland exclusion and the evergreen/grassland folders "
          "cannot be applied. 'igbp' will be saved as blank. Check "
          "df_catalog.columns and add the right alias to igbp_col above.")

df_catalog['lat_num'] = pd.to_numeric(df_catalog[lat_col], errors='coerce')
df_catalog['lon_num'] = pd.to_numeric(df_catalog[lon_col], errors='coerce')

candidate_sites_df = df_catalog[
    (df_catalog['lat_num'] > MIN_LAT) &
    df_catalog['lat_num'].notnull()
].drop_duplicates(subset=[site_col]).reset_index(drop=True)

if igbp_col is not None:
    lc_codes = candidate_sites_df[igbp_col].apply(classify_land_cover)
    excluded_mask = lc_codes.isin(EXCLUDED_LC)
    n_excluded = int(excluded_mask.sum())
    target_sites_df = candidate_sites_df[~excluded_mask].reset_index(drop=True)
    print(f"Excluded {n_excluded} cropland/wetland site(s) ({sorted(EXCLUDED_LC)}) from the candidate list.")
else:
    target_sites_df = candidate_sites_df

# Coarse pre-filter on catalog-listed data years, BEFORE downloading -
# tries an explicit column name first (varies by what listall() returns;
# not something this library's public docs specify), then falls back to
# scanning every column for the embedded year-range pattern itself (see
# find_year_range_column). If neither finds anything, this is a no-op and
# every site still gets downloaded, with the archive-filename pre-check in
# Section 5 catching disqualified ones before they're opened.
start_year_col = next((cols[k] for k in cols if k in CATALOG_YEAR_COL_ALIASES['start_year']), None)
end_year_col = next((cols[k] for k in cols if k in CATALOG_YEAR_COL_ALIASES['end_year']), None)

if start_year_col and end_year_col:
    start_yr = pd.to_numeric(target_sites_df[start_year_col], errors='coerce')
    end_yr = pd.to_numeric(target_sites_df[end_year_col], errors='coerce')
    span_years = end_yr - start_yr + 1
    catalog_ok = (end_yr >= CUTOFF_YEAR) & (span_years >= MIN_CONSECUTIVE_YEARS)
    catalog_ok = catalog_ok | start_yr.isna() | end_yr.isna()
    n_dropped = int((~catalog_ok).sum())
    target_sites_df = target_sites_df[catalog_ok].reset_index(drop=True)
    print(f"Catalog exposes data-year range ('{start_year_col}'/'{end_year_col}') - "
          f"dropped {n_dropped} site(s) BEFORE download whose listed span can't reach "
          f"{CUTOFF_YEAR}+ with >= {MIN_CONSECUTIVE_YEARS}y.")
else:
    range_col = find_year_range_column(target_sites_df)
    if range_col:
        parsed = target_sites_df[range_col].astype(str).apply(parse_year_range_from_filename)
        start_yr = parsed.apply(lambda t: t[0])
        end_yr = parsed.apply(lambda t: t[1])
        span_years = end_yr - start_yr + 1
        catalog_ok = (end_yr >= CUTOFF_YEAR) & (span_years >= MIN_CONSECUTIVE_YEARS)
        catalog_ok = catalog_ok | start_yr.isna() | end_yr.isna()
        n_dropped = int((~catalog_ok).sum())
        target_sites_df = target_sites_df[catalog_ok].reset_index(drop=True)
        print(f"Found an embedded year-range pattern in catalog column '{range_col}' - "
              f"dropped {n_dropped} site(s) BEFORE download whose listed span can't reach "
              f"{CUTOFF_YEAR}+ with >= {MIN_CONSECUTIVE_YEARS}y.")
    else:
        print("Catalog has no detectable data-year range (no named column, and no column "
              "with an embedded year-range pattern) - duration can't be pre-checked before "
              "download; short/old sites will still be downloaded and rejected in Section 5 "
              "(see the fast filename-based pre-check there, and DELETE_REJECTED_RAW_FILES "
              "if you want those archives cleaned up automatically). Run "
              "`print(df_catalog.columns.tolist())` and inspect a row or two to check "
              "whether a duration field exists under a name not covered here, then add it "
              "to CATALOG_YEAR_COL_ALIASES above.")

target_site_ids = target_sites_df[site_col].tolist()
print(f"Found {len(target_site_ids)} global FLUXNET sites with Latitude > {MIN_LAT}°N "
      f"(cropland/wetland excluded).")

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
# 5. Extract, crop to CUTOFF_YEAR+, require 10y consecutive record,
#    process target variables (incl. precipitation)
# ---------------------------------------------------------------------------
print(f"\nExtracting, cropping to {CUTOFF_YEAR}+, requiring >= {MIN_CONSECUTIVE_YEARS}y "
      "consecutive daily record, and processing target variables...")

if os.path.exists(OUTPUT_CSV):
    os.remove(OUTPUT_CSV)

zip_files = glob.glob(os.path.join(DOWNLOAD_DIR, "*.zip")) + glob.glob(os.path.join(DOWNLOAD_DIR, "*.csv"))
processed_count = 0
n_evergreen, n_grassland = 0, 0

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

        # Fast pre-check using the year range already embedded in the
        # archive's own filename - rejects sites whose overall record
        # can't possibly satisfy CUTOFF_YEAR/MIN_CONSECUTIVE_YEARS WITHOUT
        # opening or unzipping the file at all. A gap-free run is always a
        # subset of the full listed span, so this never rejects a site
        # that could otherwise have passed.
        span_start, span_end = parse_year_range_from_filename(filename)
        if span_start is not None:
            span_years = span_end - span_start + 1
            if span_end < CUTOFF_YEAR or span_years < MIN_CONSECUTIVE_YEARS:
                print(f"   -> Skipped {site_id} early (archive spans {span_start}-{span_end}, "
                      f"{span_years}y): can't reach {CUTOFF_YEAR}+ with >= {MIN_CONSECUTIVE_YEARS}y "
                      "- not opening the archive.")
                if DELETE_REJECTED_RAW_FILES:
                    os.remove(filepath)
                    print(f"      (deleted '{filename}' - DELETE_REJECTED_RAW_FILES is on)")
                continue

        meta = site_meta.get(site_id, {})
        lat, lon, igbp = meta.get('lat'), meta.get('lon'), meta.get('igbp')
        lc_code = classify_land_cover(igbp)
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

            ts_dates = pd.to_datetime(df_raw[ts_col_name].astype(str), format='%Y%m%d', errors='coerce')
            if ts_dates.isna().all():
                ts_dates = pd.to_datetime(df_raw[ts_col_name].astype(str), errors='coerce')

            # Crop to CUTOFF_YEAR onwards before anything else - pre-2013
            # flux data has no HLS imagery to pair with, so it's dropped
            # rather than carried through and trimmed later.
            keep_mask = ts_dates >= CUTOFF_DATE
            if not keep_mask.any():
                print(f"   -> Skipped {site_id}: no records at/after {CUTOFF_YEAR}")
                continue
            df_cropped = df_raw.loc[keep_mask].copy()
            ts_cropped = ts_dates.loc[keep_mask]

            # Require >= MIN_CONSECUTIVE_YEARS of gap-free daily timestamps
            # within the cropped (2013+) record - not just "the record
            # spans that many years somewhere with holes in it".
            unique_sorted_dates = sorted(ts_cropped.dropna().unique())
            unique_sorted_dates = [pd.Timestamp(d) for d in unique_sorted_dates]
            run_start, run_end, run_len = longest_consecutive_window(unique_sorted_dates)

            if run_len < MIN_CONSECUTIVE_DAYS:
                print(f"   -> Skipped {site_id}: longest consecutive run in {CUTOFF_YEAR}+ record "
                      f"is {run_len} days (< {MIN_CONSECUTIVE_DAYS} needed for {MIN_CONSECUTIVE_YEARS}y)")
                continue

            # Land-cover exclusion is already applied at the catalog stage
            # (Section 3), but re-check here too in case site_meta['igbp']
            # differs from what was in the snapshot used to build it.
            if lc_code in EXCLUDED_LC:
                print(f"   -> Skipped {site_id}: excluded land cover ({lc_code})")
                continue

            # Keep only the longest consecutive window itself.
            window_mask = (ts_cropped >= run_start) & (ts_cropped <= run_end)
            df_window = df_cropped.loc[window_mask]
            ts_window = ts_cropped.loc[window_mask]

            clean_df = pd.DataFrame()
            clean_df['site_id'] = [site_id] * len(df_window)
            clean_df['lat'] = [lat] * len(df_window)
            clean_df['lon'] = [lon] * len(df_window)
            clean_df['igbp'] = [igbp] * len(df_window)

            for target_var, aliases in TARGET_MAPPING.items():
                if target_var == 'TIMESTAMP':
                    clean_df['TIMESTAMP'] = df_window[ts_col_name].to_numpy()
                    continue
                found_col = None
                for alias in aliases:
                    if alias in existing_cols:
                        found_col = existing_cols[alias]
                        break
                clean_df[target_var] = df_window[found_col].to_numpy() if found_col else None

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
            print(f"   -> PASSED {site_id} (lc={lc_code}, {run_start.date()} to {run_end.date()}, "
                  f"{run_len} consecutive days): {len(clean_df):,} daily records saved")

            # Per-land-cover-class folders for evergreen forest / grassland.
            if lc_code in EVERGREEN_LC:
                clean_df.to_csv(EVERGREEN_DIR / f"{site_id}.csv", index=False)
                n_evergreen += 1
            elif lc_code in GRASSLAND_LC:
                clean_df.to_csv(GRASSLAND_DIR / f"{site_id}.csv", index=False)
                n_grassland += 1

    except Exception as e:
        print(f"   -> Skipped {filepath}: {e}")

print("\n" + "=" * 60)
print(f"Extraction complete! Saved {processed_count} sites "
      f"(>= {MIN_CONSECUTIVE_YEARS}y consecutive record within {CUTOFF_YEAR}+, "
      f"cropland/wetland excluded) to '{OUTPUT_CSV}'.")
print(f"Evergreen forest (ENF/EBF) sites: {n_evergreen} -> '{EVERGREEN_DIR}'")
print(f"Grassland (GRA) sites: {n_grassland} -> '{GRASSLAND_DIR}'")