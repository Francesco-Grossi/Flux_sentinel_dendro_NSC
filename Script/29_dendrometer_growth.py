"""
PIPELINE STEP 29 - Stem growth per flux site and year from dendrometers: how
much the trees grew, how fast, and when (the carbon SINK, measured directly).

Five datasets (step 14) are brought to one form. A "unit" is one tree, or
one plot x species where only plot values exist:
    DenDrought2018   daily growth per tree (zero-growth concept, um)
    Hyytiala         30-minute stem radius of two pines -> daily maximum
    Harvard Forest   band readings of every tree, 4-19 per year -> basal area
    Zoebelboden      monthly stem growth per plot and species (kg m-2)
    Tonzi (Rao)      hourly point dendrometers on five oaks -> daily maximum
                     (used if the Dryad zip was put in data_raw/dendro/rao_oaks_dryad/)

1. Cumulative growth within the year, per unit. Growth is irreversible, so
   the series is made non-decreasing (running maximum): stems also shrink
   and swell with their water content, and counting only new maxima as
   growth is the zero-growth concept of Zweifel et al. (2016). Where the
   record runs on across years the maximum is carried over, so a stem that
   shrank in a winter frost or a drought must first regain its earlier size.
   Increases before DOY COUNT_FROM or after DOY COUNT_TO are not counted:
   outside the growing season stems swell with water, they do not grow.
2. Each unit is divided by its own mean annual growth over all its years, so
   large and small trees, and datasets in different units, weigh the same.
   A value of 1 is "a normal year's growth of this tree".
3. Site-year curve = mean over units, on a daily grid (linear interpolation
   between readings). A unit-year is used only if its readings start by
   DOY FIRST_BY and run to DOY LAST_FROM.
4. From the site-year curve:
       growth_annual        total growth (1 = the site's normal year)
       growth_pre, growth_post   growth before / from the summer solstice
       frac_pre             share of the year's growth done by the solstice
       rate_pre, rate_post  mean growth rate (per day) from DOY SEASON_START
                            to the solstice / from the solstice to SEASON_END
       rate_max, doy_rate_max    fastest growth (RATE_WINDOW-day mean) and its date
       doy_g10, doy_g50, doy_g90 day when 10 / 50 / 90 % of the year's growth
                            is reached: onset, midpoint and cessation of growth
   Timing values (doy_*, rate_max) need readings no more than
   MAX_GAP_FOR_TIMING days apart in the growing season; with sparser
   readings they are left empty and only the totals and the pre / post split
   are kept (flag timing_ok).

Input : data_raw/dendro/... (step 14)
Output: data/dendro_growth_daily.csv            site-year curves, daily
        data/dendro_growth_by_site_year.csv     the metrics above
        figure/dendro_growth/<site>.png
"""
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data_raw" / "dendro"
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figure" / "dendro_growth"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DAILY = DATA_DIR / "dendro_growth_daily.csv"
OUT_YEAR = DATA_DIR / "dendro_growth_by_site_year.csv"

FIRST_BY, LAST_FROM = 135, 260          # a unit-year must have readings by / up to these DOYs
SEASON_START, SEASON_END = 100, 300     # window for the mean growth rates
RATE_WINDOW = 21                        # days, for the fastest growth
MAX_GAP_FOR_TIMING = 31                 # days between readings, DOY SEASON_START - SEASON_END
COUNT_FROM, COUNT_TO = 90, 305          # increases outside these DOYs are not counted as growth
CONTINUOUS = {'Harvard Forest HF069'}   # size series that run on across years (the others restart each year)
MIN_UNITS = 3                           # trees (or plot x species) per site-year; drops Hyytiala (two trees,
                                        # whose yearly totals differ tenfold between years)
DENDROUGHT_SITES = {83: 'CH-Dav', 35: 'CZ-BK1', 31: 'CZ-RAJ', 32: 'CZ-RAJ', 36: 'CZ-Stn', 7: 'FR-Fon'}
GRID = np.arange(1, 366)


def solstice_doy(year):
    return pd.Timestamp(year=int(year), month=6, day=21).dayofyear


# ---------------------------------------------------------------- readers: site_id, unit, date, size
def read_dendrought():
    """Daily growth per tree; gr_dmax restarts at zero each year."""
    src, cache = RAW / "dendrought2018" / "Data - Daily aggregation.xlsx", RAW / "dendrought2018" / "_daily_growth_cache.csv"
    if not src.exists():
        return None
    if cache.exists():
        d = pd.read_csv(cache, parse_dates=['date'])
    else:
        print("  reading the DenDrought2018 workbook (slow, done once) ...")
        d = pd.read_excel(src, sheet_name='Data', usecols=['site_id', 'tree_id', 'tree_species', 'data_timestamp', 'gr_dmax'])
        d = d.rename(columns={'data_timestamp': 'date'})
        d['date'] = pd.to_datetime(d['date'])
        d.to_csv(cache, index=False)
    d = d[d['site_id'].isin(DENDROUGHT_SITES)].dropna(subset=['gr_dmax'])
    return pd.DataFrame({'site_id': d['site_id'].map(DENDROUGHT_SITES), 'unit': 'tree' + d['tree_id'].astype(str),
                         'date': d['date'], 'size': d['gr_dmax'], 'dataset': 'DenDrought2018'})


def read_hyytiala():
    src = RAW / "hyytiala_liu2023" / "Hyy_30min.xlsx"
    if not src.exists():
        return None
    parts = []
    for tree in ('Pentti', 'Sylvi'):
        d = pd.read_excel(src, sheet_name=f'{tree}_30min', usecols=['Year', 'Month', 'Day', 'RBH_mm'])
        d['date'] = pd.to_datetime(d[['Year', 'Month', 'Day']].rename(columns=str.lower))
        day = d.groupby('date')['RBH_mm'].max().dropna().reset_index()          # daily maximum radius
        parts.append(pd.DataFrame({'site_id': 'FI-Hyy', 'unit': tree, 'date': day['date'], 'size': day['RBH_mm'],
                                   'dataset': 'Hyytiala (Liu et al.)'}))
    return pd.concat(parts, ignore_index=True)


def read_harvard():
    src = RAW / "harvard_forest" / "hf069-09-ems-trees.csv"
    if not src.exists():
        return None
    d = pd.read_csv(src, parse_dates=['date'])
    d = d[(d['tree_type'] == 'live') & d['basal_area_cm2'].notna()]
    return pd.DataFrame({'site_id': 'US-Ha1', 'unit': d['plottag'].astype(str), 'date': d['date'],
                         'size': d['basal_area_cm2'], 'dataset': 'Harvard Forest HF069'})


def read_zoebelboden():
    """Monthly stem growth per plot and species -> a cumulative series, dated at month end."""
    files = list((RAW / "zoebelboden_lter").glob("LTER_EU_AT_003_ZOEBELBODEN_TREEGROWTH*.txt"))
    if not files:
        return None
    d = pd.read_csv(files[0], sep='\t', encoding='latin1', decimal=',')
    d = d[d['VARIABLE'] == 'Stemgrowth'].copy()
    d['date'] = pd.to_datetime(d['TIME'] + '-01') + pd.offsets.MonthEnd(0)
    d['unit'] = 'plot' + d['STATION_CODE'].astype(str) + ' ' + d['TAXA']
    d = d.sort_values(['unit', 'date'])
    d['yr'] = d['date'].dt.year
    d['size'] = d.groupby(['unit', 'yr'])['VALUE'].cumsum()
    start = d[['unit', 'yr']].drop_duplicates()                                  # zero at 1 January
    start['date'] = pd.to_datetime(start['yr'].astype(str) + '-01-01')
    d = pd.concat([d[['unit', 'date', 'size']], start.assign(size=0.0)[['unit', 'date', 'size']]], ignore_index=True)
    return d.assign(site_id='AT-Zoe', dataset='LTER Zoebelboden')


def read_tonzi():
    """Hourly point dendrometers on five blue oaks (Rao et al., Dryad; downloaded by hand, see
    Output/dendrometer_datasets.md). 'gr' is the authors' zero-growth series."""
    name = "r_markdown/data/US_CA_Tonzi/dendrometers/US_CA_TON_20240125_clean.csv"
    folder = RAW / "rao_oaks_dryad"
    src = folder / name
    if not src.exists() and (folder / "r_markdown.zip").exists():
        with zipfile.ZipFile(folder / "r_markdown.zip") as z:
            z.extract(name, folder)
    if not src.exists():
        return None
    d = pd.read_csv(src, usecols=['time', 'tree', 'gr'])
    d['time'] = pd.to_datetime(d['time'], errors='coerce', format='mixed')
    d = d.dropna(subset=['time'])
    day = d.groupby(['tree', d['time'].dt.normalize()])['gr'].max().dropna().reset_index()
    return pd.DataFrame({'site_id': 'US-Ton', 'unit': day['tree'], 'date': day['time'], 'size': day['gr'],
                         'dataset': 'Rao et al. (Dryad)'})


parts = [p for p in (read_dendrought(), read_hyytiala(), read_harvard(), read_zoebelboden(), read_tonzi())
         if p is not None and len(p)]
if not parts:
    raise SystemExit(f"No dendrometer data in '{RAW}'. Run 14_dendrometer_download.py first.")
obs = pd.concat(parts, ignore_index=True)
obs['year'], obs['doy'] = obs['date'].dt.year, obs['date'].dt.dayofyear
print(f"{len(obs):,} readings, {obs.groupby('site_id')['unit'].nunique().to_dict()} units per site")

# ---------------------------------------------------------------- 1.-2. cumulative growth per unit-year, on the daily grid
curves, meta = [], []
for (site, unit), g in obs.groupby(['site_id', 'unit'], sort=False):
    g = g.sort_values('date')
    if g['dataset'].iloc[0] in CONTINUOUS:        # zero-growth over the whole record: a stem that shrank in a
        g = g.assign(size=g['size'].cummax())     # frost or drought first has to regain its earlier maximum
    per_year = {}
    for year, gy in g.groupby('year'):
        gy = gy.sort_values('doy').drop_duplicates('doy')
        if len(gy) < 3 or gy['doy'].iloc[0] > FIRST_BY or gy['doy'].iloc[-1] < LAST_FROM:
            continue
        cum = np.maximum.accumulate(gy['size'].to_numpy(float) - gy['size'].iloc[0])     # zero-growth: new maxima only
        doy = gy['doy'].to_numpy(float)
        # only increases between COUNT_FROM and COUNT_TO are growth; outside, stems swell with water
        lo, hi = np.interp(COUNT_FROM, doy, cum), np.interp(COUNT_TO, doy, cum)
        cum = np.clip(cum, lo, hi) - lo
        in_season = doy[(doy >= SEASON_START) & (doy <= SEASON_END)]
        gap = float(np.max(np.diff(np.concatenate([[SEASON_START], in_season, [SEASON_END]])))) if len(in_season) else np.inf
        per_year[year] = (np.interp(GRID, doy, cum), gap)       # flat before the first and after the last reading
    totals = {y: c[-1] for y, (c, _) in per_year.items()}
    norm = np.mean(list(totals.values())) if totals else 0.0
    if norm <= 0:
        continue
    for year, (c, gap) in per_year.items():
        curves.append(pd.DataFrame({'site_id': site, 'unit': unit, 'year': year, 'doy': GRID, 'growth': c / norm}))
        meta.append({'site_id': site, 'unit': unit, 'year': year, 'max_gap_days': gap})
if not curves:
    raise SystemExit("No unit-year had readings over the whole growing season.")
curves, meta = pd.concat(curves, ignore_index=True), pd.DataFrame(meta)

# ---------------------------------------------------------------- 3. site-year curves
daily = curves.groupby(['site_id', 'year', 'doy'], as_index=False).agg(growth=('growth', 'mean'), n_units=('unit', 'nunique'),
                                                                       growth_sd=('growth', 'std'))
info = meta.groupby(['site_id', 'year'], as_index=False).agg(n_units=('unit', 'nunique'), max_gap_days=('max_gap_days', 'median'))
info = info[info['n_units'] >= MIN_UNITS]
daily = daily.merge(info[['site_id', 'year']], on=['site_id', 'year'])
daily.to_csv(OUT_DAILY, index=False)
dataset = obs.drop_duplicates('site_id').set_index('site_id')['dataset']

# ---------------------------------------------------------------- 4. metrics
rows = []
for (site, year), g in daily.groupby(['site_id', 'year']):
    c = g.sort_values('doy')['growth'].to_numpy(float)
    sol, total = solstice_doy(year), float(c[-1])
    i = info[(info['site_id'] == site) & (info['year'] == year)].iloc[0]
    r = {'site_id': site, 'year': int(year), 'dataset': dataset[site], 'n_units': int(i['n_units']),
         'max_gap_days': float(i['max_gap_days']), 'timing_ok': bool(i['max_gap_days'] <= MAX_GAP_FOR_TIMING),
         'growth_annual': total}
    if total > 0:
        at = lambda d: float(c[int(d) - 1])
        r.update({'growth_pre': at(sol), 'growth_post': total - at(sol), 'frac_pre': at(sol) / total,
                  'rate_pre': (at(sol) - at(SEASON_START)) / (sol - SEASON_START),
                  'rate_post': (at(SEASON_END) - at(sol)) / (SEASON_END - sol)})
        if r['timing_ok']:
            rate = pd.Series(np.gradient(c)).rolling(RATE_WINDOW, center=True, min_periods=RATE_WINDOW // 2).mean()
            r.update({'rate_max': float(rate.max()), 'doy_rate_max': int(GRID[int(rate.idxmax())])})
            for q in (10, 50, 90):
                r[f'doy_g{q}'] = int(GRID[np.argmax(c >= q / 100 * total)])
    rows.append(r)
res = pd.DataFrame(rows)
res.to_csv(OUT_YEAR, index=False)
pd.set_option('display.width', 220)
print(f"\nSite-year growth -> '{OUT_YEAR}' ({len(res)} site-years at {res['site_id'].nunique()} sites)")
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    summ = res.groupby('site_id').agg(dataset=('dataset', 'first'), years=('year', 'size'), first=('year', 'min'),
                                      last=('year', 'max'), units=('n_units', 'median'), timing_years=('timing_ok', 'sum'),
                                      frac_pre=('frac_pre', 'mean'), onset=('doy_g10', 'mean'), cessation=('doy_g90', 'mean'),
                                      peak_rate_doy=('doy_rate_max', 'mean'))
print(summ.round(2).to_string())

# ---------------------------------------------------------------- figures
for site, g in daily.groupby('site_id'):
    years = sorted(g['year'].unique())
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cmap = plt.get_cmap('viridis')
    for k, year in enumerate(years):
        s = g[g['year'] == year].sort_values('doy')
        ax.plot(s['doy'], s['growth'], color=cmap(k / max(len(years) - 1, 1)), lw=1.2,
                label=str(year) if len(years) <= 12 else None)
    ax.axvline(172, color='gray', ls='--', lw=0.8)
    ax.set_xlabel('day of year (dashed: summer solstice)')
    ax.set_ylabel("cumulative stem growth (1 = a normal year's growth)")
    ax.set_title(f"{site} - {dataset[site]} ({years[0]}-{years[-1]}, {len(years)} years)", fontsize=10)
    if len(years) <= 12:
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{site}.png", dpi=150)
    plt.close(fig)
print(f"Figures -> '{FIG_DIR}'")
