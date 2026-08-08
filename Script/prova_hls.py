import ee
import pandas as pd
import time
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# 0. Repo-relative paths
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

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
#    Same site-level filter as prova_sentinel_1.py / prova_sentinel_2.py,
#    except CUTOFF_YEAR is loosened from 2018 to 2013 here: HLSL30 draws on
#    Landsat 8, which starts Feb 2013, so a site only needs to reach 2013
#    (not the old Sentinel-2-only 2018 constraint) to have usable HLS
#    overlap. Keep this in sync with CUTOFF_YEAR in prova_fluxent_data.py.
# ---------------------------------------------------------------------------
CUTOFF_YEAR = 2013

flux_file = DATA_DIR / "fluxnet_daily_selected_vars.csv"
if not os.path.exists(flux_file):
    raise FileNotFoundError(f"Missing '{flux_file}'. Run 1_download_fluxnet_daily.py first!")

print(f"Filtering FLUXNET sites for > 30°N, >= 5 years duration, and data reaching >= {CUTOFF_YEAR}...")
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
# 3. HLS band mapping and index definitions
# ---------------------------------------------------------------------------
# HLSL30 (Landsat 8/9 OLI, harmonized) has NO red-edge bands, so it can only
# support the "greenness" indices (NDVI, EVI, NIRv). HLSS30 (Sentinel-2
# MSI, harmonized) additionally has the green/red-edge bands needed for the
# pigment indices (ARI1, CRI1, PRI) - see project notes on Hypothesis 4.
# Reflectance bands in both collections use a 0.0001 scale factor.
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


def add_greenness_indices(refl, blue, red, nir):
    """NDVI, EVI, NIRv from harmonized reflectance bands common to L30 and S30."""
    ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI')
    evi = nir.subtract(red).multiply(2.5).divide(
        nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1)
    ).rename('EVI')
    nirv = ndvi.multiply(nir).rename('NIRv')
    return ndvi.addBands([evi, nirv])


def add_pigment_indices(green, blue, re1):
    """ARI1, CRI1, PRI - only computable from S30 (needs green/red-edge bands).
    PRI is an approximation: S2/HLS-S30 has no true 531 nm band, so blue (B2,
    ~492 nm) and green (B3, ~560 nm) are used as the closest available proxy,
    consistent with common practice in the S2-phenology literature."""
    ari1 = ee.Image(1).divide(green).subtract(ee.Image(1).divide(re1)).rename('ARI1')
    cri1 = ee.Image(1).divide(blue).subtract(ee.Image(1).divide(green)).rename('CRI1')
    pri = blue.subtract(green).divide(blue.add(green)).rename('PRI')
    return ari1.addBands([cri1, pri])


def make_l30_mapper(roi):
    def _map(img):
        mask = fmask_valid(img)
        refl = img.select(list(L30_BANDS.values())).updateMask(mask).multiply(0.0001)
        blue, green, red, nir = [refl.select(L30_BANDS[b]) for b in ['blue', 'green', 'red', 'nir']]
        indices = add_greenness_indices(refl, blue, red, nir)

        stats = indices.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9)
        valid_frac = mask.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9).get('Fmask')

        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'),
            'sensor': 'L30',
            'NDVI': stats.get('NDVI'), 'EVI': stats.get('EVI'), 'NIRv': stats.get('NIRv'),
            'ARI1': None, 'CRI1': None, 'PRI': None,
            'valid_pixel_frac': valid_frac
        })
    return _map


def make_s30_mapper(roi):
    def _map(img):
        mask = fmask_valid(img)
        refl = img.select(list(S30_BANDS.values())).updateMask(mask).multiply(0.0001)
        blue, green, red, nir, re1 = [refl.select(S30_BANDS[b]) for b in ['blue', 'green', 'red', 'nir', 're1']]
        indices = add_greenness_indices(refl, blue, red, nir).addBands(add_pigment_indices(green, blue, re1))

        stats = indices.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9)
        valid_frac = mask.reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=30, maxPixels=1e9).get('Fmask')

        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'),
            'sensor': 'S30',
            'NDVI': stats.get('NDVI'), 'EVI': stats.get('EVI'), 'NIRv': stats.get('NIRv'),
            'ARI1': stats.get('ARI1'), 'CRI1': stats.get('CRI1'), 'PRI': stats.get('PRI'),
            'valid_pixel_frac': valid_frac
        })
    return _map


# ---------------------------------------------------------------------------
# 4. HLS Extraction Loop
# ---------------------------------------------------------------------------
output_csv = DATA_DIR / "fluxnet_all_highlat_hls_indices.csv"
if os.path.exists(output_csv):
    os.remove(output_csv)

CLOUD_COVERAGE_MAX = 80  # scene-level pre-filter; per-pixel Fmask does the real work


def extract_hls_site(site_id, lat, lon, radius_m=500, start_date='2013-01-01', end_date='2021-01-01'):
    point = ee.Geometry.Point([lon, lat])
    roi = point.buffer(radius_m)

    # NOTE: Earth Engine re-registered these collections under the 'v002'
    # suffix (the old 'v2.0' asset IDs no longer resolve and raise
    # "ImageCollection asset ... not found").
    l30_col = (ee.ImageCollection('NASA/HLS/HLSL30/v002')
               .filterBounds(roi)
               .filterDate(start_date, end_date)
               .filter(ee.Filter.lt('CLOUD_COVERAGE', CLOUD_COVERAGE_MAX))
               .map(make_l30_mapper(roi)))

    s30_col = (ee.ImageCollection('NASA/HLS/HLSS30/v002')
               .filterBounds(roi)
               .filterDate(start_date, end_date)
               .filter(ee.Filter.lt('CLOUD_COVERAGE', CLOUD_COVERAGE_MAX))
               .map(make_s30_mapper(roi)))

    merged = l30_col.merge(s30_col)
    features = merged.getInfo()['features']

    records = []
    for f in features:
        p = f['properties']
        if p.get('NDVI') is None:
            continue
        records.append({
            'site_id': site_id, 'lat': lat, 'lon': lon,
            'date': p['date'], 'sensor': p['sensor'],
            'NDVI': p['NDVI'], 'EVI': p['EVI'], 'NIRv': p['NIRv'],
            'ARI1': p['ARI1'], 'CRI1': p['CRI1'], 'PRI': p['PRI'],
            'valid_pixel_frac': p.get('valid_pixel_frac')
        })
    return records


total = len(target_sites)
for idx, row in target_sites.iterrows():
    s_id, lat, lon = str(row['site_id']), float(row['lat']), float(row['lon'])

    if s_id not in site_date_range.index:
        print(f"[{idx+1}/{total}] Skipping HLS: {s_id} - no valid flux dates found")
        continue

    d_start = site_date_range.loc[s_id, 'min'].strftime('%Y-%m-%d')
    d_end = (site_date_range.loc[s_id, 'max'] + pd.Timedelta(days=1)).strftime('%Y-%m-%d')

    print(f"[{idx+1}/{total}] Processing HLS: {s_id} (Flux span used: {d_start} to {d_end})")

    try:
        recs = extract_hls_site(s_id, lat, lon, start_date=d_start, end_date=d_end)
        if recs:
            df_s = (pd.DataFrame(recs)
                    .groupby(['site_id', 'lat', 'lon', 'sensor', 'date']).mean().reset_index())
            file_exists = os.path.exists(output_csv)
            df_s.to_csv(output_csv, mode='a', header=not file_exists, index=False)
            n_l30 = (df_s['sensor'] == 'L30').sum()
            n_s30 = (df_s['sensor'] == 'S30').sum()
            print(f"   -> Saved {len(df_s)} HLS records (L30: {n_l30}, S30: {n_s30})")
        else:
            print("   -> No usable HLS imagery found")
    except Exception as e:
        print(f"   -> Skipped {s_id}: {e}")

    time.sleep(0.2)

print(f"\nHLS extraction complete! Saved output to '{output_csv}'.")
print("Note: ARI1/CRI1/PRI are populated only for S30 rows - L30 (Landsat) has no red-edge bands.")