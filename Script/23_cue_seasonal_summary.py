"""
PIPELINE STEP 23 - Describe the CUE estimates of step 22 at three time scales
and check how they relate to each other:

  annual    one CUE per site-year (the Luo et al. method)
  seasonal  mean of the daily CUE before / after the summer solstice
            (CUE_pre: DOY 80 - solstice, CUE_post: solstice - DOY 300)
  daily     the daily CUE series interpolated from the sliding windows

Questions:
  1. What does the seasonal course of CUE look like, per leaf habit?
  2. Do the seasonal means carry information beyond the annual value?
     (within-site correlation between annual CUE and CUE_pre / CUE_post, and
     between CUE_pre and CUE_post)
  3. How different is NPP from the daily CUE (NPPd) from NPP from the annual
     CUE (NPP), for the pre- and post-solstice sums?
The EOS relationships themselves are tested by steps 12-17, which pick up
NPP, NPPd and CUEd (window-mean daily CUE) automatically.

Input : data/cue_luo2025_site_year.csv, cue_luo2025_seasonal.csv, npp_luo2025_daily.csv
Output: data/cue_seasonal_curve_by_leaf_habit.csv
        data/cue_seasonal_means_by_site_year.csv
        data/cue_scale_agreement.csv
        figure/cue_seasonal/seasonal_curve.png, annual_vs_seasonal.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

CUE_CSV = ec.DATA_DIR / "cue_luo2025_site_year.csv"
SEASON_CSV = ec.DATA_DIR / "cue_luo2025_seasonal.csv"
OUT_CURVE = ec.DATA_DIR / "cue_seasonal_curve_by_leaf_habit.csv"
OUT_MEANS = ec.DATA_DIR / "cue_seasonal_means_by_site_year.csv"
OUT_AGREE = ec.DATA_DIR / "cue_scale_agreement.csv"
FIG_DIR = ec.FIGURE_DIR / "cue_seasonal"
FIG_DIR.mkdir(parents=True, exist_ok=True)
ec.require(CUE_CSV, SEASON_CSV, ec.NPP_CSV, hint="Run 22_cue_luo2025.py first.")

PRE_START, POST_END = 80, 300
MIN_DAYS = 30                 # days with a value needed for a seasonal mean / sum
MIN_SITES_PER_DOY = 5
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous',
              'ENF': 'evergreen', 'EBF': 'evergreen'}   # everything else: grass/shrub

igbp = pd.read_csv(ec.DATA_DIR / "fluxnet_daily_all_vars.csv", usecols=['site_id', 'igbp']).drop_duplicates('site_id')
habit = igbp.set_index('site_id')['igbp'].map(LEAF_HABIT).fillna('grass/shrub')

# ---------------------------------------------------------------- 1. seasonal curve
season = pd.read_csv(SEASON_CSV)
season = season[season['converged'].astype(bool)].copy()
season['leaf_habit'] = season['site_id'].map(habit)
site_curve = season.groupby(['leaf_habit', 'site_id', 'doy_center'])['CUE'].median().reset_index()
curve = site_curve.groupby(['leaf_habit', 'doy_center'])['CUE'].agg(
    median='median', q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75), n_sites='size').reset_index()
curve = curve[curve['n_sites'] >= MIN_SITES_PER_DOY]
curve.to_csv(OUT_CURVE, index=False)
print(f"Seasonal CUE curve (median across sites) -> '{OUT_CURVE}'")
show = curve[curve['doy_center'].isin([75, 105, 135, 165, 195, 225, 255, 285])]
print(show.pivot(index='leaf_habit', columns='doy_center', values='median').round(2).to_string())

fig, ax = plt.subplots(figsize=(7.5, 4.2))
for h, c in (('deciduous', '#2f855a'), ('evergreen', '#2b6cb0'), ('grass/shrub', '#b7791f')):
    d = curve[curve['leaf_habit'] == h]
    if len(d):
        ax.plot(d['doy_center'], d['median'], color=c, label=f"{h} ({int(d['n_sites'].max())} sites)")
        ax.fill_between(d['doy_center'], d['q25'], d['q75'], color=c, alpha=0.15)
ax.axvline(172, color='0.5', lw=0.8, ls='--')
ax.set_xlabel('day of year')
ax.set_ylabel('CUE (30-day windows)')
ax.set_title('Seasonal course of CUE (median and interquartile range across sites)')
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG_DIR / "seasonal_curve.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- 2. annual vs seasonal means
d = pd.read_csv(ec.NPP_CSV, parse_dates=['date'])
d['year'], d['doy'] = d['date'].dt.year, d['date'].dt.dayofyear
sol = d['year'].map(ec.solstice_doy)
d['part'] = np.where((d['doy'] >= PRE_START) & (d['doy'] < sol), 'pre',
                     np.where((d['doy'] >= sol) & (d['doy'] <= POST_END), 'post', None))
g = d.dropna(subset=['part']).groupby(['site_id', 'year', 'part'])
agg = g.agg(CUE_mean=('CUE_daily', 'mean'), n_cue=('CUE_daily', 'count'),
            NPP_sum=('NPP', 'sum'), n_npp=('NPP', 'count'),
            NPPd_sum=('NPPd', 'sum'), n_nppd=('NPPd', 'count')).reset_index()
agg.loc[agg['n_cue'] < MIN_DAYS, 'CUE_mean'] = np.nan
agg.loc[agg['n_npp'] < MIN_DAYS, 'NPP_sum'] = np.nan
agg.loc[agg['n_nppd'] < MIN_DAYS, 'NPPd_sum'] = np.nan
wide = agg.pivot(index=['site_id', 'year'], columns='part', values=['CUE_mean', 'NPP_sum', 'NPPd_sum'])
wide.columns = [f'{a.replace("_mean", "").replace("_sum", "")}_{b}' for a, b in wide.columns]
means = pd.read_csv(CUE_CSV)[['site_id', 'year', 'CUE', 'CUE_sd']].rename(
    columns={'CUE': 'CUE_annual', 'CUE_sd': 'CUE_annual_sd'}).merge(wide.reset_index(), on=['site_id', 'year'], how='left')
means['leaf_habit'] = means['site_id'].map(habit)
means.to_csv(OUT_MEANS, index=False)
print(f"\nSeasonal means per site-year -> '{OUT_MEANS}' "
      f"({int(means['CUE_pre'].notna().sum())} site-years with CUE_pre, {int(means['CUE_post'].notna().sum())} with CUE_post)")
print(means[['CUE_annual', 'CUE_pre', 'CUE_post']].describe().loc[['count', 'mean', 'std', '25%', '50%', '75%']].round(3).to_string())


def within_r(df, a, b):
    x = df[['site_id', a, b]].dropna()
    x = x[x.groupby('site_id')['site_id'].transform('size') >= 3]
    if len(x) < 15 or x['site_id'].nunique() < 3:
        return None
    dm = x[[a, b]] - x.groupby('site_id')[[a, b]].transform('mean')
    return {'a': a, 'b': b, 'n_site_years': len(x), 'n_sites': x['site_id'].nunique(),
            'r_within': float(np.corrcoef(dm[a], dm[b])[0, 1]), 'r_pooled': float(np.corrcoef(x[a], x[b])[0, 1])}


pairs = [('CUE_annual', 'CUE_pre'), ('CUE_annual', 'CUE_post'), ('CUE_pre', 'CUE_post'),
         ('NPP_pre', 'NPPd_pre'), ('NPP_post', 'NPPd_post')]
rows = []
for hab in ['all', 'deciduous', 'evergreen', 'grass/shrub']:
    sub = means if hab == 'all' else means[means['leaf_habit'] == hab]
    for a, b in pairs:
        r = within_r(sub, a, b)
        if r:
            rows.append({'leaf_habit': hab, **r})
agree = pd.DataFrame(rows)
agree.to_csv(OUT_AGREE, index=False)
print(f"\nAgreement between time scales (within-site r) -> '{OUT_AGREE}'")
print(agree.round(2).to_string(index=False))

fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
for ax, (a, b) in zip(axes, pairs[:3]):
    x = means[['site_id', a, b]].dropna()
    dm = x[[a, b]] - x.groupby('site_id')[[a, b]].transform('mean')
    ax.scatter(dm[a], dm[b], s=6, alpha=0.4)
    ax.axhline(0, color='0.8', lw=0.8)
    ax.axvline(0, color='0.8', lw=0.8)
    r = np.corrcoef(dm[a], dm[b])[0, 1] if len(dm) > 2 else np.nan
    ax.set_title(f"within-site r = {r:.2f}  (n = {len(dm)})", fontsize=9)
    ax.set_xlabel(f"{a} anomaly")
    ax.set_ylabel(f"{b} anomaly")
fig.tight_layout()
fig.savefig(FIG_DIR / "annual_vs_seasonal.png", dpi=150)
plt.close(fig)
print(f"\nFigures -> '{FIG_DIR}'")
