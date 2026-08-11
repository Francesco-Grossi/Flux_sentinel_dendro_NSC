import ee
import numpy as np
import pandas as pd
import time
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# 0. Repo-relative paths
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

RAW_BANDS_CSV = DATA_DIR / "fluxnet_all_highlat_hls_raw_bands.csv"
INDICES_CSV = DATA_DIR / "fluxnet_all_highlat_hls_indices.csv"

# ---------------------------------------------------------------------------
# 1. Initialize Earth Engine
# ---------------------------------------------------------------------------
PROJECT_ID = 'prova-sentinel'  # Replace with your GEE project ID

try:
    ee.Initialize(project=PROJECT_ID)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=PROJECT_ID)

# ---------------------------------------------------------------------------
# 2. Filter Sites: > 30°N, >= 5 Years Span, Last Year >= CUTOFF_YEAR
# ---------------------------------------------------------------------------
CUTOFF_YEAR = 2013

flux_file = DATA_DIR / "fluxnet_daily_selected_vars.csv"
if not os.path.exists(flux_file):
    raise FileNotFoundError(f"Missing '{flux_file}'. Run 1_download_fluxnet_daily.py first!")

print(f"Filtering FLUXNET sites for > 30°N, >= 10 years duration, and data reaching >= {CUTOFF_YEAR}...")
df_flux = pd.read_csv(flux_file)

df_flux['date'] = pd.to_datetime(df_flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df_flux['date'] = df_flux['date'].fillna(pd.to_datetime(df_flux['TIMESTAMP'].astype(str), errors='coerce'))

valid_sites = []
for site_id, group in df_flux.groupby('site_id'):
    lat = group['lat'].dropna().iloc[0] if not group['lat'].dropna().empty else None
    lon = group['lon'].dropna().iloc[0] if not group['lon'].dropna().empty else None

    if lat is None or lon is None or lat <= 30.0:
        continue

    dates = group['date'].dropna()
    if dates.empty:
        continue

    min_date = dates.min()
    max_date = dates.max()
    duration_days = (max_date - min_date).days

    if duration_days >= 1825 and max_date.year >= CUTOFF_YEAR:
        valid_sites.append({
            'site_id': site_id,
            'lat': lat,
            'lon': lon,
            'start_year': min_date.year,
            'end_year': max_date.year,
            'duration_years': round(duration_days / 365.25, 1)
        })

target_sites = pd.DataFrame(valid_sites)
print(f"Found {len(target_sites)} qualifying sites (Lat > 30°N, >= 5 yrs data, reaching >= {CUTOFF_YEAR}).\n")

site_date_range = df_flux.dropna(subset=['date']).groupby('site_id')['date'].agg(['min', 'max'])

# ---------------------------------------------------------------------------
# 3. HLS band mapping - RAW BANDS ONLY, no indices computed here
# ---------------------------------------------------------------------------
# HLSL30 (Landsat 8/9 OLI, harmonized) has NO red-edge band, so it only ever
# gives blue/green/red/nir. HLSS30 (Sentinel-2 MSI, harmonized) additionally
# has the red-edge band (re1) needed for the pigment indices (ARI1, CRI1).
L30_BANDS = {'blue': 'B2', 'green': 'B3', 'red': 'B4', 'nir': 'B5'}
S30_BANDS = {'blue': 'B2', 'green': 'B3', 'red': 'B4', 'nir': 'B8A', 're1': 'B5'}

# Fmask quality band bit layout (HLS v2.0):
#   bit1 = cloud, bit2 = adjacent to cloud/shadow, bit3 = cloud shadow, bit4 = snow/ice
CLOUD_BIT, ADJ_BIT, SHADOW_BIT, SNOW_BIT = 1, 2, 3, 4


def fmask_valid(img):
    """Boolean mask: True where pixel is clear of cloud/shadow/adjacency/snow."""
    fmask = img.select('Fmask')
    bad = (fmask.bitwiseAnd(1 << CLOUD_BIT).neq(0)
           .Or(fmask.bitwiseAnd(1 << ADJ_BIT).neq(0))
           .Or(fmask.bitwiseAnd(1 << SHADOW_BIT).neq(0))
           .Or(fmask.bitwiseAnd(1 << SNOW_BIT).neq(0)))
    return bad.Not()


# NOTE ON SCALING: earlier versions of this pipeline multiplied the raw HLS
# DN values by 0.0001 here, at the Earth Engine step, on the assumption that
# HLS bands are always stored as reflectance*10000. That assumption held for
# the old 'v2.0' asset IDs, but doesn't necessarily hold for whatever GEE
# collection is actually being read (asset re-registrations, like the
# 'v2.0' -> 'v002' rename already noted below, are exactly the kind of thing
# that can silently change this). Downstream, EVI computed from those
# "scaled" bands came out ~10,000x too small while NDVI (which is a ratio,
# so scale-invariant) looked perfectly normal - a classic symptom of the
# bands actually having already been reflectance-scaled upstream, and this
# script scaling them down by *another* 0.0001.
#
# To avoid guessing again, this version downloads the completely RAW,
# UNSCALED band values and lets the processing step (Section 6 below)
# inspect their actual magnitude and decide the scale factor from the data
# itself, rather than assuming it blind at extraction time.
def make_l30_mapper(roi):
    def _map(img):
        mask = fmask_valid(img)
        refl = img.select(list(L30_BANDS.values())).updateMask(mask)
        bands = {b: refl.select(L30_BANDS[b]) for b in ['blue', 'green', 'red', 'nir']}

        stats = ee.Image.cat(list(bands.values())).rename(list(bands.keys())) \
            .reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9)
        valid_frac = mask.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9).get('Fmask')

        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'),
            'sensor': 'L30',
            'blue': stats.get('blue'), 'green': stats.get('green'),
            'red': stats.get('red'), 'nir': stats.get('nir'), 're1': None,
            'valid_pixel_frac': valid_frac
        })
    return _map


def make_s30_mapper(roi):
    def _map(img):
        mask = fmask_valid(img)
        refl = img.select(list(S30_BANDS.values())).updateMask(mask)
        bands = {b: refl.select(S30_BANDS[b]) for b in ['blue', 'green', 'red', 'nir', 're1']}

        stats = ee.Image.cat(list(bands.values())).rename(list(bands.keys())) \
            .reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9)
        valid_frac = mask.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9).get('Fmask')

        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'),
            'sensor': 'S30',
            'blue': stats.get('blue'), 'green': stats.get('green'),
            'red': stats.get('red'), 'nir': stats.get('nir'), 're1': stats.get('re1'),
            'valid_pixel_frac': valid_frac
        })
    return _map


# ---------------------------------------------------------------------------
# 4. Raw-band extraction loop
# ---------------------------------------------------------------------------
if os.path.exists(RAW_BANDS_CSV):
    os.remove(RAW_BANDS_CSV)

CLOUD_COVERAGE_MAX = 80  # scene-level pre-filter; per-pixel Fmask does the real work

# S30 (Sentinel-2) is preferred over L30 (Landsat) whenever available, since
# it carries the extra red-edge band needed for ARI1/CRI1/PRI. L30 only
# fills temporal gaps where no S30 observation falls within this many days.
GAP_FILL_WINDOW_DAYS = 5


def _fetch_features(collection_id, mapper_factory, roi, start_date, end_date, min_chunk_days=7):
    """getInfo() on a mapped ImageCollection aborts once the query
    accumulates over 5000 elements. Query in date-range chunks, halving on
    failure, until it fits (or min_chunk_days is hit, at which point the
    underlying error is something else and gets raised)."""
    try:
        col = (ee.ImageCollection(collection_id)
               .filterBounds(roi)
               .filterDate(start_date, end_date)
               .filter(ee.Filter.lt('CLOUD_COVERAGE', CLOUD_COVERAGE_MAX))
               .map(mapper_factory(roi)))
        return col.getInfo()['features']
    except Exception as e:
        if '5000 elements' not in str(e):
            raise
        start_ts, end_ts = pd.Timestamp(start_date), pd.Timestamp(end_date)
        span_days = (end_ts - start_ts).days
        if span_days <= min_chunk_days:
            raise
        mid = (start_ts + pd.Timedelta(days=span_days // 2)).strftime('%Y-%m-%d')
        first_half = _fetch_features(collection_id, mapper_factory, roi, start_date, mid, min_chunk_days)
        second_half = _fetch_features(collection_id, mapper_factory, roi, mid, end_date, min_chunk_days)
        return first_half + second_half


def extract_raw_bands_site(site_id, lat, lon, radius_m=500, start_date='2013-01-01', end_date='2021-01-01'):
    point = ee.Geometry.Point([lon, lat])
    roi = point.buffer(radius_m)

    # NOTE: Earth Engine re-registered these collections under the 'v002'
    # suffix (the old 'v2.0' asset IDs no longer resolve).
    s30_features = _fetch_features('NASA/HLS/HLSS30/v002', make_s30_mapper, roi, start_date, end_date)
    l30_features = _fetch_features('NASA/HLS/HLSL30/v002', make_l30_mapper, roi, start_date, end_date)

    # getInfo() drops feature properties entirely when their server-side
    # value is None (e.g. re1 on L30 scenes) rather than keeping the key
    # with a null value - read with .get(), not direct indexing.
    def _to_records(features):
        recs = []
        for f in features:
            p = f['properties']
            if p.get('nir') is None or p.get('red') is None:
                continue
            recs.append({
                'site_id': site_id, 'lat': lat, 'lon': lon,
                'date': p.get('date'), 'sensor': p.get('sensor'),
                'blue': p.get('blue'), 'green': p.get('green'),
                'red': p.get('red'), 'nir': p.get('nir'), 're1': p.get('re1'),
                'valid_pixel_frac': p.get('valid_pixel_frac')
            })
        return recs

    s30_records = _to_records(s30_features)
    l30_records = _to_records(l30_features)

    if not s30_records:
        return l30_records
    if not l30_records:
        return s30_records

    s30_dates = pd.to_datetime([r['date'] for r in s30_records])
    kept_l30 = []
    for r in l30_records:
        d = pd.to_datetime(r['date'])
        gaps = abs((s30_dates - d).days.to_numpy())
        min_gap = gaps.min() if len(gaps) else np.inf
        if min_gap > GAP_FILL_WINDOW_DAYS:
            kept_l30.append(r)

    return s30_records + kept_l30


total = len(target_sites)
for idx, row in target_sites.iterrows():
    s_id, lat, lon = str(row['site_id']), float(row['lat']), float(row['lon'])

    if s_id not in site_date_range.index:
        print(f"[{idx+1}/{total}] Skipping HLS: {s_id} - no valid flux dates found")
        continue

    d_start = site_date_range.loc[s_id, 'min'].strftime('%Y-%m-%d')
    d_end = (site_date_range.loc[s_id, 'max'] + pd.Timedelta(days=1)).strftime('%Y-%m-%d')

    print(f"[{idx+1}/{total}] Fetching raw bands: {s_id} (Flux span used: {d_start} to {d_end})")

    try:
        recs = extract_raw_bands_site(s_id, lat, lon, start_date=d_start, end_date=d_end)
        if recs:
            df_s = (pd.DataFrame(recs)
                    .groupby(['site_id', 'lat', 'lon', 'sensor', 'date']).mean().reset_index())
            file_exists = os.path.exists(RAW_BANDS_CSV)
            df_s.to_csv(RAW_BANDS_CSV, mode='a', header=not file_exists, index=False)
            n_l30 = (df_s['sensor'] == 'L30').sum()
            n_s30 = (df_s['sensor'] == 'S30').sum()
            print(f"   -> Saved {len(df_s)} raw-band records (L30: {n_l30}, S30: {n_s30})")
        else:
            print("   -> No usable HLS imagery found")
    except Exception as e:
        print(f"   -> Skipped {s_id}: {e}")

    time.sleep(0.2)

print(f"\nRaw-band extraction complete! Saved to '{RAW_BANDS_CSV}'.")


# ---------------------------------------------------------------------------
# 5. Determine the reflectance scale factor FROM THE DATA
# ---------------------------------------------------------------------------
# HLS surface reflectance is nominally stored as DN = reflectance * 10000
# (scale factor 0.0001), but that's an assumption about the asset, not a
# physical constant - and the 'v2.0' -> 'v002' rename above is a reminder
# that asset internals can change under us. Rather than hardcoding
# `* 0.0001` again, look at how big the raw NIR values actually are and
# infer the right factor: physical surface reflectance is in [0, ~1.2] for
# a well-behaved pixel, so:
#   - raw values in the thousands  -> DN encoding, scale by 0.0001
#   - raw values already in [0, ~2] -> already reflectance, scale by 1.0
def infer_reflectance_scale(raw_band_values):
    finite = raw_band_values.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return 1.0
    typical_magnitude = finite.abs().median()
    return 0.0001 if typical_magnitude > 10 else 1.0


# ---------------------------------------------------------------------------
# 6. Compute indices from the raw bands (pandas, not Earth Engine)
# ---------------------------------------------------------------------------
print("\nComputing indices from raw bands...")

raw = pd.read_csv(RAW_BANDS_CSV)
band_cols = ['blue', 'green', 'red', 'nir', 're1']

scale_factor = infer_reflectance_scale(raw['nir'])
print(f"Inferred reflectance scale factor: {scale_factor} "
      f"(median raw NIR = {raw['nir'].dropna().abs().median():.4g})")

for col in band_cols:
    raw[col] = raw[col] * scale_factor

# Physically implausible reflectance after scaling (allow a little negative
# slack for atmospheric-correction noise, same convention as similar
# filters elsewhere in this pipeline) -> NaN rather than feeding the index
# formulas garbage.
for col in band_cols:
    raw.loc[(raw[col] < -0.05) | (raw[col] > 1.2), col] = np.nan

blue, green, red, nir, re1 = (raw[c] for c in band_cols)

raw['NDVI'] = (nir - red) / (nir + red)
raw['EVI'] = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
raw['NIRv'] = raw['NDVI'] * nir

# ARI1, CRI1, PRI only defined where re1 exists (S30 rows) - reciprocal
# terms blow up as any reflectance band approaches 0, which is expected
# behavior for this index family, not a bug; downstream fitting scripts
# already handle the resulting outliers (percentile clip + rescale).
has_re1 = re1.notna()
raw['ARI1'] = np.where(has_re1, (1.0 / green) - (1.0 / re1), np.nan)
raw['CRI1'] = np.where(has_re1, (1.0 / blue) - (1.0 / green), np.nan)
# PRI approximation: HLS/S2 has no true 531nm band, so blue (~492nm) and
# green (~560nm) are used as the closest available proxy.
raw['PRI'] = np.where(has_re1, (blue - green) / (blue + green), np.nan)

indices_df = raw[['site_id', 'lat', 'lon', 'sensor', 'date',
                   'NDVI', 'EVI', 'NIRv', 'ARI1', 'CRI1', 'PRI', 'valid_pixel_frac']]
indices_df.to_csv(INDICES_CSV, index=False)

print(f"Indices written to '{INDICES_CSV}' ({len(indices_df)} rows).")
print("Note: ARI1/CRI1/PRI are populated only for S30 rows - L30 (Landsat) has no red-edge band.")