"""
PIPELINE STEP 2 - HLS (Landsat 8/9 HLSL30 + Sentinel-2 HLSS30) raw-band
extraction + NDVI/EVI/NIRv.

For each qualifying FLUXNET site (>30N, >=5y span - re-derived from the
daily file here so this step can run independently of step 1's exact site
list), pulls cloud/shadow/snow-masked surface reflectance from BOTH HLS
collections over the site's own flux date range. HLS L30 and S30 are
already harmonized (common grid, BRDF- and bandpass-adjusted), so the two
sensors can be pooled into one time series; S30 roughly triples the number
of clear observations per season from 2016 onwards.

Footprint: mean over a RADIUS_M buffer around the tower, restricted to the
pixels whose ESA WorldCover class matches the site's IGBP vegetation type
(e.g. tree cover for forests, grassland for GRA), so that roads, fields,
water or buildings inside the buffer do not contaminate the signal. If
fewer than MIN_LC_MATCH_FRAC of the buffer matches, all pixels are used
and lc_masked=False is recorded.

Snow: snow pixels are masked out of the reflectance means, but the per-image
fraction of snow-covered (clear, non-cloud) pixels is kept as snow_frac,
even for images with no usable reflectance. Step 5 uses these dates as
dormant-season (background) observations.

Outputs: data/fluxnet_all_highlat_landsat_raw_bands.csv
         data/fluxnet_all_highlat_landsat_indices.csv   (NDVI, EVI, NIRv, snow_frac)
Resumable: rerunning skips sites already in the raw-band file (delete the
file, or set HLS_RESTART=1, to start from scratch).
"""
import ee
import numpy as np
import pandas as pd
import time
import os
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RAW_BANDS_CSV = DATA_DIR / "fluxnet_all_highlat_landsat_raw_bands.csv"
INDICES_CSV = DATA_DIR / "fluxnet_all_highlat_landsat_indices.csv"

PROJECT_ID = 'prova-sentinel'  # replace with your GEE project ID
try:
    ee.Initialize(project=PROJECT_ID)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=PROJECT_ID)

RADIUS_M = 1000             # buffer radius around the tower (was 500 m, all pixels)
MIN_LC_MATCH_FRAC = 0.2     # below this share of matching land cover, fall back to all pixels
CLOUD_COVERAGE_MAX = 80
SCALE_M = 30

# HLS band names per collection. S30 uses the narrow NIR (B8A), which is the
# band HLS harmonizes to Landsat's B5.
COLLECTIONS = {
    'L30': ('NASA/HLS/HLSL30/v002', {'blue': 'B2', 'green': 'B3', 'red': 'B4', 'nir': 'B5'}),
    'S30': ('NASA/HLS/HLSS30/v002', {'blue': 'B2', 'green': 'B3', 'red': 'B4', 'nir': 'B8A'}),
}
CLOUD_BIT, ADJ_BIT, SHADOW_BIT, SNOW_BIT = 1, 2, 3, 4

# IGBP class -> ESA WorldCover v200 classes that represent the same vegetation
# (10 tree cover, 20 shrubland, 30 grassland, 100 moss & lichen).
IGBP_TO_WORLDCOVER = {
    'ENF': [10], 'EBF': [10], 'DNF': [10], 'DBF': [10], 'MF': [10],
    'CSH': [20], 'OSH': [20, 100],
    'WSA': [10, 20, 30], 'SAV': [10, 20, 30],
    'GRA': [30, 100],
}
WORLDCOVER = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map')

flux_file = DATA_DIR / "fluxnet_daily_all_vars.csv"
if not os.path.exists(flux_file):
    raise FileNotFoundError(f"Missing '{flux_file}'. Run 01_fluxnet_download.py first.")

print("Filtering FLUXNET sites for > 30 deg N and >= 5 years duration...")
df_flux = pd.read_csv(flux_file, usecols=['site_id', 'lat', 'lon', 'igbp', 'TIMESTAMP'])
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
    duration_days = (dates.max() - dates.min()).days
    if duration_days >= 1825:
        igbp = group['igbp'].dropna().iloc[0] if not group['igbp'].dropna().empty else None
        valid_sites.append({'site_id': site_id, 'lat': lat, 'lon': lon, 'igbp': igbp})

target_sites = pd.DataFrame(valid_sites).reset_index(drop=True)
print(f"Found {len(target_sites)} qualifying sites.\n")
site_date_range = df_flux.dropna(subset=['date']).groupby('site_id')['date'].agg(['min', 'max'])


def landcover_mask(roi, igbp):
    """Pixels of the site's own vegetation type, or None if too few match."""
    classes = IGBP_TO_WORLDCOVER.get(str(igbp))
    if not classes:
        return None, np.nan
    match = WORLDCOVER.remap(classes, [1] * len(classes), 0).rename('lc')
    frac = match.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=SCALE_M,
                              maxPixels=1e9).get('lc').getInfo()
    if frac is None or frac < MIN_LC_MATCH_FRAC:
        return None, frac
    return match.selfMask(), frac


def make_mapper(sensor, band_map, roi, lc_mask):
    def _map(img):
        fmask = img.select('Fmask')
        cloudy = (fmask.bitwiseAnd(1 << CLOUD_BIT).neq(0)
                  .Or(fmask.bitwiseAnd(1 << ADJ_BIT).neq(0))
                  .Or(fmask.bitwiseAnd(1 << SHADOW_BIT).neq(0)))
        snow = fmask.bitwiseAnd(1 << SNOW_BIT).neq(0)
        valid = cloudy.Or(snow).Not()
        if lc_mask is not None:
            valid, cloudy, snow = valid.updateMask(lc_mask), cloudy.updateMask(lc_mask), snow.updateMask(lc_mask)
        refl = img.select(list(band_map.values()), list(band_map.keys())).updateMask(valid.selfMask())
        stats = refl.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=SCALE_M, maxPixels=1e9)
        valid_frac = valid.rename('v').reduceRegion(ee.Reducer.mean(), roi, SCALE_M, maxPixels=1e9).get('v')
        # snow share among the clear (non-cloud) pixels: tells step 5 the ground is snow-covered
        snow_frac = snow.updateMask(cloudy.Not()).rename('s') \
            .reduceRegion(ee.Reducer.mean(), roi, SCALE_M, maxPixels=1e9).get('s')
        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'), 'sensor': sensor,
            'blue': stats.get('blue'), 'green': stats.get('green'),
            'red': stats.get('red'), 'nir': stats.get('nir'),
            'valid_pixel_frac': valid_frac, 'snow_frac': snow_frac,
        })
    return _map


def _fetch_features(collection_id, mapper, roi, start_date, end_date, min_chunk_days=7):
    try:
        col = (ee.ImageCollection(collection_id).filterBounds(roi).filterDate(start_date, end_date)
               .filter(ee.Filter.lt('CLOUD_COVERAGE', CLOUD_COVERAGE_MAX)).map(mapper))
        return col.getInfo()['features']
    except Exception as e:
        if '5000 elements' not in str(e) and 'memory' not in str(e).lower():
            raise
        start_ts, end_ts = pd.Timestamp(start_date), pd.Timestamp(end_date)
        span_days = (end_ts - start_ts).days
        if span_days <= min_chunk_days:
            raise
        mid = (start_ts + pd.Timedelta(days=span_days // 2)).strftime('%Y-%m-%d')
        return (_fetch_features(collection_id, mapper, roi, start_date, mid, min_chunk_days)
                + _fetch_features(collection_id, mapper, roi, mid, end_date, min_chunk_days))


def extract_raw_bands_site(site_id, lat, lon, igbp, start_date, end_date):
    roi = ee.Geometry.Point([lon, lat]).buffer(RADIUS_M)
    lc_mask, lc_frac = landcover_mask(roi, igbp)
    recs = []
    for sensor, (collection_id, band_map) in COLLECTIONS.items():
        mapper = make_mapper(sensor, band_map, roi, lc_mask)
        for f in _fetch_features(collection_id, mapper, roi, start_date, end_date):
            p = f['properties']
            has_refl = p.get('nir') is not None and p.get('red') is not None
            if not has_refl and not p.get('snow_frac'):
                continue  # fully cloudy image: carries no information
            recs.append({'site_id': site_id, 'lat': lat, 'lon': lon, 'igbp': igbp,
                         'date': p.get('date'), 'sensor': p.get('sensor'),
                         'blue': p.get('blue'), 'green': p.get('green'),
                         'red': p.get('red'), 'nir': p.get('nir'),
                         'valid_pixel_frac': p.get('valid_pixel_frac'), 'snow_frac': p.get('snow_frac'),
                         'radius_m': RADIUS_M, 'lc_masked': lc_mask is not None, 'lc_match_frac': lc_frac})
    return recs


done_sites = set()
if os.path.exists(RAW_BANDS_CSV):
    existing_cols = pd.read_csv(RAW_BANDS_CSV, nrows=0).columns
    if os.environ.get('HLS_RESTART') == '1' or 'snow_frac' not in existing_cols:
        os.remove(RAW_BANDS_CSV)  # old-format (L30-only) file or explicit restart
    else:
        done_sites = set(pd.read_csv(RAW_BANDS_CSV, usecols=['site_id'])['site_id'].unique())
        print(f"Resuming: {len(done_sites)} sites already extracted.")

total = len(target_sites)
for idx, row in target_sites.iterrows():
    s_id, lat, lon, igbp = str(row['site_id']), float(row['lat']), float(row['lon']), row['igbp']
    if s_id not in site_date_range.index or s_id in done_sites:
        continue
    d_start = site_date_range.loc[s_id, 'min'].strftime('%Y-%m-%d')
    d_end = (site_date_range.loc[s_id, 'max'] + pd.Timedelta(days=1)).strftime('%Y-%m-%d')
    if d_end < '2013-04-01':
        continue  # flux record ends before HLS starts
    print(f"[{idx+1}/{total}] Fetching HLS L30+S30: {s_id} [{igbp}] ({d_start} to {d_end})")
    try:
        recs = extract_raw_bands_site(s_id, lat, lon, igbp, d_start, d_end)
        if recs:
            df_s = pd.DataFrame(recs)
            file_exists = os.path.exists(RAW_BANDS_CSV)
            df_s.to_csv(RAW_BANDS_CSV, mode='a', header=not file_exists, index=False)
            print(f"   -> Saved {len(df_s)} records "
                  f"(L30 {(df_s.sensor == 'L30').sum()}, S30 {(df_s.sensor == 'S30').sum()}; "
                  f"land-cover mask: {bool(df_s.lc_masked.iloc[0])}, match {df_s.lc_match_frac.iloc[0]})")
    except Exception as e:
        print(f"   -> Skipped {s_id}: {e}")
    time.sleep(0.2)

print(f"\nRaw-band extraction complete -> '{RAW_BANDS_CSV}'.")


def infer_reflectance_scale(raw_band_values):
    finite = raw_band_values.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return 1.0
    return 0.0001 if finite.abs().median() > 10 else 1.0


print("\nComputing indices (NDVI, EVI, NIRv) from raw bands...")
raw = pd.read_csv(RAW_BANDS_CSV)
band_cols = ['blue', 'green', 'red', 'nir']
scale_factor = infer_reflectance_scale(raw['nir'])
print(f"Inferred reflectance scale factor: {scale_factor}")
for col in band_cols:
    raw[col] = raw[col] * scale_factor
for col in band_cols:
    raw.loc[(raw[col] < -0.05) | (raw[col] > 1.2), col] = np.nan

blue, green, red, nir = (raw[c] for c in band_cols)
raw['NDVI'] = (nir - red) / (nir + red)
raw['EVI'] = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
raw['NIRv'] = raw['NDVI'] * nir

indices_df = raw[['site_id', 'lat', 'lon', 'igbp', 'sensor', 'date', 'NDVI', 'EVI', 'NIRv',
                  'valid_pixel_frac', 'snow_frac', 'radius_m', 'lc_masked', 'lc_match_frac']]
indices_df.to_csv(INDICES_CSV, index=False)
print(f"Indices written to '{INDICES_CSV}' ({len(indices_df)} rows; "
      f"L30 {(indices_df.sensor == 'L30').sum()}, S30 {(indices_df.sensor == 'S30').sum()}).")
