"""
PIPELINE STEP 20 - PhenoCam canopy greenness (GCC) at the flux sites, and the
phenology (leaf-out / EOS) fitted to it. PhenoCam is the standard ground
reference for validating satellite land-surface phenology: a camera on or
next to the tower, daily, unaffected by clouds, independent of GPP.

1. Match cameras to flux sites via the PhenoCam API (phenocam.nau.edu/api):
   a camera matches when its metadata lists the FLUXNET site ID
   (flux_sitenames) or it stands within MATCH_KM of the tower.
2. Pick one region of interest (ROI) per site: the vegetation type that
   corresponds to the site's IGBP class (e.g. DB for DBF, EN for ENF, GR for
   GRA), then the longest record.
3. Download that ROI's 1-day summary file (cached in data_raw/phenocam/) and
   take gcc_90 (90th percentile of daytime GCC, the PhenoCam standard),
   dropping days flagged as outliers.
4. Fit with the shared routine in pheno_fit.py - NOT PhenoCam's own
   transition dates, which exist only for 10/25/50% of amplitude and use a
   different curve; fitting here gives EOS90/EOS50/EOS10 defined exactly as
   for the satellite (step 5) and tower (step 19) series.

Only QC-passing flux site-years are fitted (those are the ones with GPP
predictors), and the flux record supplies TA_F for the frozen-day flag.
Days with PhenoCam snow_flag set count as snow-covered.

Input : data/fluxnet_landsat_merged.csv; PhenoCam API (network)
Output: data/phenocam_site_matches.csv
        data/phenocam_gcc_daily.csv
        data/phenology_phenocam_by_site_year_index.csv   (same columns as step 5)
"""
import io
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import pheno_fit as pf

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
CACHE_DIR = ROOT / "data_raw" / "phenocam"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"
SITE_CSV = DATA_DIR / "fluxnet_daily_all_vars.csv"
MATCH_CSV = DATA_DIR / "phenocam_site_matches.csv"
DAILY_CSV = DATA_DIR / "phenocam_gcc_daily.csv"
OUTPUT_CSV = DATA_DIR / "phenology_phenocam_by_site_year_index.csv"

API = "https://phenocam.nau.edu/api"
MATCH_KM = 1.0
GCC_COL = 'gcc_90'
# IGBP class -> PhenoCam ROI vegetation types, in order of preference
IGBP_TO_ROI = {
    'DBF': ['DB', 'MX'], 'DNF': ['DN', 'DB'], 'MF': ['DB', 'MX', 'EN'],
    'ENF': ['EN', 'MX'], 'EBF': ['EB', 'MX'],
    'GRA': ['GR'], 'CSH': ['SH', 'GR'], 'OSH': ['SH', 'GR'],
    'WSA': ['DB', 'EN', 'EB', 'GR', 'SH'], 'SAV': ['GR', 'DB', 'EB', 'SH'],
}
NON_VEG_ROI = {'AG', 'NV', 'RF', 'UN', 'WL', 'XX'}

if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run steps 1-4 first.")


def api_all(endpoint):
    url, out = f"{API}/{endpoint}/?format=json&limit=500", []
    while url:
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        page = r.json()
        out.extend(page['results'])
        url = page.get('next')
    return out


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


# ------------------------------------------------------------------ matching
flux = pd.read_csv(FLUX_CSV, usecols=['site_id', 'TIMESTAMP', 'TA_F'])
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux = flux.dropna(subset=['date']).sort_values(['site_id', 'date'])
flux = pf.add_frozen_flag(flux)
site_meta = pd.read_csv(SITE_CSV, usecols=['site_id', 'lat', 'lon', 'igbp']).drop_duplicates('site_id')
site_meta = site_meta[site_meta['site_id'].isin(flux['site_id'].unique())].reset_index(drop=True)

print("Reading the PhenoCam camera and ROI lists...")
cameras = api_all('cameras')
rois = pd.DataFrame(api_all('roilists'))
rois['site_years'] = pd.to_numeric(rois['site_years'], errors='coerce')
print(f"  {len(cameras)} cameras, {len(rois)} ROIs")

cam = pd.DataFrame([{
    'camera': c['Sitename'], 'cam_lat': c['Lat'], 'cam_lon': c['Lon'],
    'flux_ids': set(re.split(r'[,;\s]+', ((c.get('sitemetadata') or {}).get('flux_sitenames') or '').strip())) - {''},
} for c in cameras]).dropna(subset=['cam_lat', 'cam_lon'])

matches = []
for _, s in site_meta.iterrows():
    dist = haversine_km(s['lat'], s['lon'], cam['cam_lat'].to_numpy(float), cam['cam_lon'].to_numpy(float))
    named = cam['flux_ids'].apply(lambda ids: s['site_id'] in ids).to_numpy()
    cand = cam[named | (dist <= MATCH_KM)].assign(dist_km=dist[named | (dist <= MATCH_KM)], by_name=named[named | (dist <= MATCH_KM)])
    r = rois[rois['site'].isin(cand['camera']) & ~rois['roitype'].isin(NON_VEG_ROI)].merge(
        cand[['camera', 'dist_km', 'by_name']], left_on='site', right_on='camera')
    if r.empty:
        continue
    pref = IGBP_TO_ROI.get(str(s['igbp']), [])
    r['type_rank'] = r['roitype'].apply(lambda t: pref.index(t) if t in pref else len(pref))
    best = r.sort_values(['type_rank', 'site_years', 'dist_km'], ascending=[True, False, True]).iloc[0]
    matches.append({'site_id': s['site_id'], 'igbp': s['igbp'], 'camera': best['site'], 'roi_name': best['roi_name'],
                    'roitype': best['roitype'], 'roi_matches_igbp': best['roitype'] in pref,
                    'dist_km': round(float(best['dist_km']), 3), 'matched_by_flux_id': bool(best['by_name']),
                    'first_date': best['first_date'], 'last_date': best['last_date'],
                    'site_years': best['site_years'], 'one_day_summary': best['one_day_summary']})

matches = pd.DataFrame(matches)
if matches.empty:
    raise SystemExit("No flux site has a PhenoCam camera.")
matches.to_csv(MATCH_CSV, index=False)
print(f"{len(matches)}/{len(site_meta)} flux sites have a PhenoCam ROI -> '{MATCH_CSV}'")

# ----------------------------------------------------- download + fit GCC
daily_parts, rows = [], []
for i, m in matches.iterrows():
    cache = CACHE_DIR / f"{m['roi_name']}_1day.csv"
    if not cache.exists():
        r = requests.get(m['one_day_summary'], timeout=300)
        if r.status_code != 200:
            print(f"  {m['site_id']}: could not download {m['roi_name']} (HTTP {r.status_code})")
            continue
        cache.write_bytes(r.content)
        time.sleep(0.5)
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
