"""
PIPELINE STEP 13 - Download PhenoCam canopy greenness (GCC) for the flux sites.
PhenoCam is the standard ground reference for validating satellite
land-surface phenology: a camera on or next to the tower, daily, unaffected
by clouds, independent of GPP. The phenology is fitted in step 25.

1. Match cameras to flux sites via the PhenoCam API (phenocam.nau.edu/api):
   a camera matches when its metadata lists the FLUXNET site ID
   (flux_sitenames), or it stands within MATCH_KM of the tower, or within
   NEAR_KM and it has an ROI of the site's own vegetation type. The best
   of these three tiers that applies is used.
2. Pick the regions of interest (ROIs): the vegetation type that best
   corresponds to the site's IGBP class (e.g. DB for DBF, EN for ENF, GR for
   GRA), and ALL ROIs of that type. PhenoCam starts a new ROI file
   (..._1000, _2000, _3000) whenever a camera's field of view shifts, and a
   site can have several cameras, so one ROI often covers only a few of the
   flux years. Step 25 picks, year by year, the ROI with the most data.
3. Download each ROI's 1-day summary file to data_raw/phenocam/. Files
   already there are kept (delete one to refresh it).

PhenoCam covers mainly North America; European flux sites with their own
cameras (ICOS) are not in this archive.

Input : data/fluxnet_daily_all_vars.csv (step 11: site list, coordinates, IGBP);
        PhenoCam API (network)
Output: data/phenocam_site_matches.csv
        data_raw/phenocam/<roi>_1day.csv
"""
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
CACHE_DIR = ROOT / "data_raw" / "phenocam"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
SITE_CSV = DATA_DIR / "fluxnet_daily_all_vars.csv"
MATCH_CSV = DATA_DIR / "phenocam_site_matches.csv"

API = "https://phenocam.nau.edu/api"
MATCH_KM = 1.0              # any vegetation ROI of a camera this close is accepted
NEAR_KM = 2.0               # up to here only ROIs of the site's own vegetation type
# IGBP class -> PhenoCam ROI vegetation types, in order of preference
IGBP_TO_ROI = {
    'DBF': ['DB', 'MX'], 'DNF': ['DN', 'DB'], 'MF': ['DB', 'MX', 'EN'],
    'ENF': ['EN', 'MX'], 'EBF': ['EB', 'MX'],
    'GRA': ['GR'], 'CSH': ['SH', 'GR'], 'OSH': ['SH', 'GR'],
    'WSA': ['DB', 'EN', 'EB', 'GR', 'SH'], 'SAV': ['GR', 'DB', 'EB', 'SH'],
}
NON_VEG_ROI = {'AG', 'NV', 'RF', 'UN', 'WL', 'XX'}

if not os.path.exists(SITE_CSV):
    raise FileNotFoundError(f"Missing '{SITE_CSV}'. Run 11_fluxnet_download.py first.")


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


site_meta = pd.read_csv(SITE_CSV, usecols=['site_id', 'lat', 'lon', 'igbp']).drop_duplicates('site_id') \
    .reset_index(drop=True)

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
    sel = named | (dist <= NEAR_KM)
    cand = cam[sel].assign(dist_km=dist[sel], by_name=named[sel])
    r = rois[rois['site'].isin(cand['camera']) & ~rois['roitype'].isin(NON_VEG_ROI)].merge(
        cand[['camera', 'dist_km', 'by_name']], left_on='site', right_on='camera')
    pref = IGBP_TO_ROI.get(str(s['igbp']), [])
    r['type_rank'] = r['roitype'].apply(lambda t: pref.index(t) if t in pref else len(pref))
    # cameras between MATCH_KM and NEAR_KM count only if they look at the site's own vegetation type
    r = r[r['by_name'] | (r['dist_km'] <= MATCH_KM) | (r['type_rank'] < len(pref))]
    if r.empty:
        continue
    # the tier a camera was matched on, then the vegetation type; ALL ROIs of the best tier and type are kept
    r['tier'] = np.where(r['by_name'], 0, np.where(r['dist_km'] <= MATCH_KM, 1, 2))
    r = r[r['tier'] == r['tier'].min()]
    r = r[r['type_rank'] == r['type_rank'].min()].sort_values(['site_years', 'dist_km'], ascending=[False, True])
    for priority, (_, b) in enumerate(r.iterrows()):
        matches.append({'site_id': s['site_id'], 'igbp': s['igbp'], 'camera': b['site'], 'roi_name': b['roi_name'],
                        'roitype': b['roitype'], 'roi_matches_igbp': b['roitype'] in pref, 'priority': priority,
                        'dist_km': round(float(b['dist_km']), 3), 'matched_by_flux_id': bool(b['by_name']),
                        'match_tier': ['flux site ID', f'within {MATCH_KM:g} km',
                                       f'within {NEAR_KM:g} km, same vegetation type'][int(b['tier'])],
                        'first_date': b['first_date'], 'last_date': b['last_date'],
                        'site_years': b['site_years'], 'one_day_summary': b['one_day_summary']})

matches = pd.DataFrame(matches)
if matches.empty:
    raise SystemExit("No flux site has a PhenoCam camera.")

n_new = 0
for _, m in matches.iterrows():
    cache = CACHE_DIR / f"{m['roi_name']}_1day.csv"
    if cache.exists():
        continue
    r = requests.get(m['one_day_summary'], timeout=300)
    if r.status_code != 200:
        print(f"  {m['site_id']}: could not download {m['roi_name']} (HTTP {r.status_code})")
        continue
    cache.write_bytes(r.content)
    n_new += 1
    print(f"  {m['site_id']} <- {m['roi_name']} ({m['dist_km']} km)")
    time.sleep(0.5)

matches['downloaded'] = [(CACHE_DIR / f"{r}_1day.csv").exists() for r in matches['roi_name']]
matches.to_csv(MATCH_CSV, index=False)
print(f"{matches['site_id'].nunique()}/{len(site_meta)} flux sites have a PhenoCam ROI "
      f"({len(matches)} ROIs in total) -> '{MATCH_CSV}'")
print(matches.drop_duplicates('site_id')['match_tier'].value_counts().to_string())
print(f"{n_new} file(s) downloaded, {int(matches['downloaded'].sum())} available in '{CACHE_DIR}'.")
