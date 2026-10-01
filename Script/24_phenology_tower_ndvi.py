"""
PIPELINE STEP 19 - Tower-based greenness: broadband NDVI from the flux
tower's own radiation sensors, and the phenology (leaf-out / EOS) fitted
to them. Independent of GPP and of satellites, daily, cloud-free, and from
the tower footprint itself - a reference for the satellite EOS of step 5.

Source: the half-hourly (HH) or hourly (HR) FLUXMET file inside each site's
FLUXNET zip in data_raw/ (already downloaded by step 1). No new download.

Daily values use midday records only (MIDDAY_HOURS, local standard time),
where the sun is high and the sensors' cosine response is reliable:
    rho_PAR = sum(PPFD_OUT) / sum(PPFD_IN)
    rho_NIR = sum(SW_OUT - k*PPFD_OUT) / sum(SW_IN - k*PPFD_IN)
    NDVI    = (rho_NIR - rho_PAR) / (rho_NIR + rho_PAR)
with k = PAR_W_PER_UMOL converting PPFD (umol m-2 s-1) to W m-2, so that
shortwave minus PAR approximates the near-infrared part (broadband NDVI:
Huemmrich et al. 1999; Wang et al. 2004; Jenkins et al. 2007). Only measured
SW_IN is used (SW_IN_F with SW_IN_F_QC == 0), never gap-filled values.

Snow: a day is 'snowy' when the shortwave reflectance sum(SW_OUT)/sum(SW_IN)
exceeds the site's summer (JJA) median by SNOW_REFL_EXCESS while air
temperature is below SNOW_MAX_TA_C. Snowy days are excluded from the fit and
filled with the dormant background (pheno_fit). This ratio is used ONLY as the
snow flag; it is not fitted or written out.

NDVI_tower is fitted with the shared routine in pheno_fit.py. It needs SW_OUT,
PPFD_IN and PPFD_OUT (~half of the sites).

Input : data_raw/*_<site>_FLUXNET_*.zip, data/fluxnet_landsat_merged.csv (TA_F,
        and the QC-passing site-years to fit)
Output: data/tower_ndvi_daily.csv
        data/phenology_tower_by_site_year_index.csv   (same columns as step 5)
"""
import glob
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pheno_fit as pf

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
RAW_DIR = ROOT / "data_raw"
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"
DAILY_CSV = DATA_DIR / "tower_ndvi_daily.csv"
OUTPUT_CSV = DATA_DIR / "phenology_tower_by_site_year_index.csv"

MIDDAY_HOURS = (10, 14)        # [start, end), local standard time
MIN_MIDDAY_FRAC = 0.5          # share of midday records that must be valid
MIN_SW_IN = 150.0              # W m-2
MIN_PPFD_IN = 300.0            # umol m-2 s-1
PAR_W_PER_UMOL = 0.219         # 1 / 4.57 (McCree 1972)
SNOW_REFL_EXCESS = 0.10
SNOW_MAX_TA_C = 5.0
SMOOTH_DAYS = 5                # centred rolling median on the daily series before fitting
RAD_COLS = ['SW_IN_F', 'SW_IN_F_QC', 'SW_OUT', 'PPFD_IN', 'PPFD_OUT']

if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run steps 1-4 first.")


def read_subdaily(zip_path):
    """Half-hourly/hourly radiation columns of one FLUXNET zip, or None."""
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if 'FLUXMET' in n and re.search(r'_(HH|HR)_', n) and n.endswith('.csv')]
        if not names:
            return None
        with z.open(names[0]) as f:
            header = pd.read_csv(f, nrows=0).columns
        use = ['TIMESTAMP_START'] + [c for c in RAD_COLS if c in header]
        if not {'SW_OUT', 'SW_IN_F', 'PPFD_IN', 'PPFD_OUT'} <= set(use):
            return None
        with z.open(names[0]) as f:
            d = pd.read_csv(f, usecols=use, na_values=[-9999, -9999.0])
    return d


def daily_from_subdaily(d):
    ts = d['TIMESTAMP_START'].astype('int64')
    d['date'] = pd.to_datetime((ts // 10000).astype(str), format='%Y%m%d', errors='coerce')
    hour = (ts // 100) % 100
    d = d[(hour >= MIDDAY_HOURS[0]) & (hour < MIDDAY_HOURS[1])].dropna(subset=['date']).copy()
    n_midday = d.groupby('date').size().median()   # 8 for half-hourly, 4 for hourly files

    sw_in = d['SW_IN_F'].where(d['SW_IN_F_QC'] == 0) if 'SW_IN_F_QC' in d else d['SW_IN_F']
    ok_sw = (sw_in >= MIN_SW_IN) & (d['SW_OUT'] > 0) & (d['SW_OUT'] < sw_in)
    a = pd.DataFrame({'date': d['date'], 'sw_in': sw_in.where(ok_sw), 'sw_out': d['SW_OUT'].where(ok_sw),
                      'n_sw': ok_sw.astype(int)})
    nir_in = sw_in - PAR_W_PER_UMOL * d['PPFD_IN']
    nir_out = d['SW_OUT'] - PAR_W_PER_UMOL * d['PPFD_OUT']
    ok_par = ok_sw & (d['PPFD_IN'] >= MIN_PPFD_IN) & (d['PPFD_OUT'] > 0) & (d['PPFD_OUT'] < d['PPFD_IN']) \
        & (nir_in > 0) & (nir_out > 0)
    a['par_in'], a['par_out'] = d['PPFD_IN'].where(ok_par), d['PPFD_OUT'].where(ok_par)
    a['nir_in'], a['nir_out'] = nir_in.where(ok_par), nir_out.where(ok_par)
    a['n_ndvi'] = ok_par.astype(int)

    g = a.groupby('date').sum(min_count=1)
    out = pd.DataFrame(index=g.index)
    # shortwave reflectance: snow flag only
    out['sw_refl'] = (g['sw_out'] / g['sw_in']).where(g['n_sw'] >= MIN_MIDDAY_FRAC * n_midday)
    enough = g['n_ndvi'] >= MIN_MIDDAY_FRAC * n_midday
    rho_par = (g['par_out'] / g['par_in']).where(enough)
    rho_nir = (g['nir_out'] / g['nir_in']).where(enough)
    ok = (rho_par > 0) & (rho_par < 1) & (rho_nir > 0) & (rho_nir < 1)
    out['rho_par'], out['rho_nir'] = rho_par.where(ok), rho_nir.where(ok)
    out['ndvi_tower'] = ((rho_nir - rho_par) / (rho_nir + rho_par)).where(ok)
    return out.reset_index()


flux = pd.read_csv(FLUX_CSV, usecols=['site_id', 'TIMESTAMP', 'TA_F'])
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux = flux.dropna(subset=['date']).sort_values(['site_id', 'date'])
flux = pf.add_frozen_flag(flux)
sites = sorted(flux['site_id'].unique())
zips = glob.glob(str(RAW_DIR / "*.zip"))
print(f"Tower radiation for {len(sites)} sites with QC-passing flux years...")

daily_parts, rows, skipped = [], [], {}
for si, site_id in enumerate(sites, start=1):
    match = [z for z in zips if re.search(rf'(^|[_-]){re.escape(site_id)}([_-]|$)', os.path.basename(z))]
    if not match:
        skipped[site_id] = 'no zip in data_raw'
        continue
    try:
        sub = read_subdaily(match[0])
    except zipfile.BadZipFile:
        skipped[site_id] = 'corrupt zip (re-download with step 1)'
        continue
    if sub is None:
        skipped[site_id] = 'no SW_OUT / PPFD_OUT'
        continue
    day = daily_from_subdaily(sub)
    if day['ndvi_tower'].notna().sum() == 0:
        skipped[site_id] = 'no valid SW_OUT / PPFD_OUT data'
        continue

    sdf = flux[flux['site_id'] == site_id].merge(day, on='date', how='left')
    sdf['year'], sdf['doy'] = sdf['date'].dt.year, sdf['date'].dt.dayofyear
    summer_refl = sdf.loc[sdf['date'].dt.month.isin([6, 7, 8]), 'sw_refl'].median()
    sdf['snowy'] = (sdf['sw_refl'] > summer_refl + SNOW_REFL_EXCESS) & (sdf['TA_F'] < SNOW_MAX_TA_C)
    sdf['fit_weight'] = 1.0
    clear = sdf['ndvi_tower'].where(~sdf['snowy'])
    sdf['ndvi_tower_smooth'] = clear.rolling(SMOOTH_DAYS, center=True, min_periods=2).median().where(clear.notna())

    daily_parts.append(sdf.loc[sdf['ndvi_tower'].notna() | sdf['snowy'],
                               ['site_id', 'date', 'rho_par', 'rho_nir', 'ndvi_tower', 'snowy']])
    n_before = len(rows)
    rows.extend(pf.fit_site_index(site_id, sdf, 'ndvi_tower_smooth', label='NDVI_tower'))
    print(f"  [{si}/{len(sites)}] {site_id}: {int(sdf['ndvi_tower'].notna().sum())} NDVI days, "
          f"{int(sdf['snowy'].sum())} snowy, {len(rows) - n_before} site-year fits")

if not rows:
    raise SystemExit("No site had usable SW_OUT / PPFD_OUT data.")
pd.concat(daily_parts, ignore_index=True).to_csv(DAILY_CSV, index=False)
print(f"\nDaily tower NDVI -> '{DAILY_CSV}'")
pf.write_results(rows, OUTPUT_CSV)
if skipped:
    reasons = pd.Series(skipped).value_counts().to_dict()
    print(f"\nSkipped {len(skipped)} sites: {reasons}")
    corrupt = [s for s, r in skipped.items() if r.startswith('corrupt')]
    if corrupt:
        print(f"  corrupt zips: {corrupt}")
