"""
PIPELINE STEP 2 - Landsat HLS (HLSL30) raw-band extraction + NDVI/EVI/NIRv.

For each qualifying FLUXNET site (>30N, >=5y span - re-derived from the
daily file here so this step can run independently of step 1's exact site
list), pulls cloud/shadow/snow-masked Landsat 8/9 surface reflectance over
the site's own flux date range, then computes indices from the RAW
(unscaled) bands with a data-inferred reflectance scale factor.

Output: data/fluxnet_all_highlat_landsat_raw_bands.csv
        data/fluxnet_all_highlat_landsat_indices.csv   (NDVI, EVI, NIRv)
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

flux_file = DATA_DIR / "fluxnet_daily_all_vars.csv"
if not os.path.exists(flux_file):
    raise FileNotFoundError(f"Missing '{flux_file}'. Run 01_fluxnet_download.py first.")

print("Filtering FLUXNET sites for > 30 deg N and >= 5 years duration...")
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
    duration_days = (dates.max() - dates.min()).days
    if duration_days >= 1825:
        valid_sites.append({'site_id': site_id, 'lat': lat, 'lon': lon})

target_sites = pd.DataFrame(valid_sites)
print(f"Found {len(target_sites)} qualifying sites.\n")
site_date_range = df_flux.dropna(subset=['date']).groupby('site_id')['date'].agg(['min', 'max'])

L30_BANDS = {'blue': 'B2', 'green': 'B3', 'red': 'B4', 'nir': 'B5'}
CLOUD_BIT, ADJ_BIT, SHADOW_BIT, SNOW_BIT = 1, 2, 3, 4


def fmask_valid(img):
    fmask = img.select('Fmask')
    bad = (fmask.bitwiseAnd(1 << CLOUD_BIT).neq(0)
           .Or(fmask.bitwiseAnd(1 << ADJ_BIT).neq(0))
           .Or(fmask.bitwiseAnd(1 << SHADOW_BIT).neq(0))
           .Or(fmask.bitwiseAnd(1 << SNOW_BIT).neq(0)))
    return bad.Not()


def make_l30_mapper(roi):
    def _map(img):
        mask = fmask_valid(img)
        refl = img.select(list(L30_BANDS.values())).updateMask(mask)
        bands = {b: refl.select(L30_BANDS[b]) for b in ['blue', 'green', 'red', 'nir']}
        stats = ee.Image.cat(list(bands.values())).rename(list(bands.keys())) \
            .reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9)
        valid_frac = mask.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9).get('Fmask')
        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'), 'sensor': 'L30',
            'blue': stats.get('blue'), 'green': stats.get('green'),
            'red': stats.get('red'), 'nir': stats.get('nir'), 'valid_pixel_frac': valid_frac
        })
    return _map


if os.path.exists(RAW_BANDS_CSV):
    os.remove(RAW_BANDS_CSV)
CLOUD_COVERAGE_MAX = 80


def _fetch_features(collection_id, mapper_factory, roi, start_date, end_date, min_chunk_days=7):
    try:
        col = (ee.ImageCollection(collection_id).filterBounds(roi).filterDate(start_date, end_date)
               .filter(ee.Filter.lt('CLOUD_COVERAGE', CLOUD_COVERAGE_MAX)).map(mapper_factory(roi)))
        return col.getInfo()['features']
    except Exception as e:
        if '5000 elements' not in str(e):
            raise
        start_ts, end_ts = pd.Timestamp(start_date), pd.Timestamp(end_date)
        span_days = (end_ts - start_ts).days
        if span_days <= min_chunk_days:
            raise
        mid = (start_ts + pd.Timedelta(days=span_days // 2)).strftime('%Y-%m-%d')
        return (_fetch_features(collection_id, mapper_factory, roi, start_date, mid, min_chunk_days)
                + _fetch_features(collection_id, mapper_factory, roi, mid, end_date, min_chunk_days))


def extract_raw_bands_site(site_id, lat, lon, radius_m=500, start_date='2000-01-01', end_date='2025-01-01'):
    roi = ee.Geometry.Point([lon, lat]).buffer(radius_m)
    l30_features = _fetch_features('NASA/HLS/HLSL30/v002', make_l30_mapper, roi, start_date, end_date)

    def _to_records(features):
        recs = []
        for f in features:
            p = f['properties']
            if p.get('nir') is None or p.get('red') is None:
                continue
            recs.append({'site_id': site_id, 'lat': lat, 'lon': lon, 'date': p.get('date'),
                         'sensor': p.get('sensor'), 'blue': p.get('blue'), 'green': p.get('green'),
                         'red': p.get('red'), 'nir': p.get('nir'), 'valid_pixel_frac': p.get('valid_pixel_frac')})
        return recs
    return _to_records(l30_features)


total = len(target_sites)
for idx, row in target_sites.iterrows():
    s_id, lat, lon = str(row['site_id']), float(row['lat']), float(row['lon'])
    if s_id not in site_date_range.index:
        continue
    d_start = site_date_range.loc[s_id, 'min'].strftime('%Y-%m-%d')
    d_end = (site_date_range.loc[s_id, 'max'] + pd.Timedelta(days=1)).strftime('%Y-%m-%d')
    print(f"[{idx+1}/{total}] Fetching Landsat: {s_id} ({d_start} to {d_end})")
    try:
        recs = extract_raw_bands_site(s_id, lat, lon, start_date=d_start, end_date=d_end)
        if recs:
            df_s = pd.DataFrame(recs).groupby(['site_id', 'lat', 'lon', 'sensor', 'date']).mean().reset_index()
            file_exists = os.path.exists(RAW_BANDS_CSV)
            df_s.to_csv(RAW_BANDS_CSV, mode='a', header=not file_exists, index=False)
            print(f"   -> Saved {len(df_s)} records")
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

indices_df = raw[['site_id', 'lat', 'lon', 'sensor', 'date', 'NDVI', 'EVI', 'NIRv', 'valid_pixel_frac']]
indices_df.to_csv(INDICES_CSV, index=False)
print(f"Indices written to '{INDICES_CSV}' ({len(indices_df)} rows).")
