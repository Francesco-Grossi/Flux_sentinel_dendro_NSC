"""
PIPELINE STEP 45 - Measured stem growth (the sink) against GPP (the source)
and against the onset of senescence.

The hypothesis behind "early GPP advances senescence" is a sink one: when
the trees have grown enough, the canopy is no longer needed. GPP is only a
proxy for that. Step 29 gives growth itself, from dendrometers, at the few
flux sites that have them. Two questions, both within sites (year minus site
mean, standard errors clustered by site; slopes per +1 within-site SD):

A. Is growth coupled to GPP? Annual growth against annual GPP, and growth
   before the solstice against GPP before the solstice (the CAL_PRE days
   before it, and leaf-out is not needed).
B. Does growth predict the onset of senescence (EOS90)? Predictors:
       growth_pre     growth done by the solstice
       rate_pre       mean growth rate, DOY 100 to the solstice
       growth_annual  the year's total growth
       frac_pre       share of the year's growth done by the solstice
       doy_g90        date growth stops (only years with readings close
                      enough in time; see step 29)
   and, for comparison on the same site-years, pre-solstice GPP and spring
   temperature.
   All EOS sources are also stacked (anomalies relative to each site x source
   mean), because one source alone has few of these site-years.

C. Carbon budget. For the sites where the stand's biomass growth is known in
   carbon units (US-Ha1, AT-Zoe), the dendrometer curve is scaled to the
   stand's mean annual growth and combined with the flux-derived respiration:
       respiration  Ra = (1 - CUE) x GPP        (daily CUE of step 26)
       total sink   = growth + Ra
       NPP          = GPP - Ra
       residual     = GPP - total sink = NPP - growth
   The residual is carbon that was fixed but neither respired nor put into
   the measured growth: roots, (at US-Ha1) leaves, and reserves. Each term is
   summed before and after the solstice and tested against EOS90.
   Caveats: growth is aboveground only, and Ra rests on the CUE estimate, so
   the residual carries the errors of both.

This is a case study: a handful of sites, and US-Ha1 supplies most of the
years. Sites need at least MIN_YEARS years of both growth and EOS90. With
so few sites the standard errors cannot be clustered by site; they are
heteroskedasticity-robust only, and with all sources stacked the same year
enters once per source, so those p-values are too small. Read the signs and
sizes, not the significance.

Input : data/dendro_growth_by_site_year.csv (step 29), phenology tables, daily flux
Output: data/growth_vs_gpp.csv
        data/growth_vs_senescence.csv
        data/carbon_budget_by_site_year.csv
        data/carbon_budget_vs_senescence.csv
        figure/dendro_growth/growth_vs_gpp.png, growth_vs_eos90.png, carbon_budget.png
        Output/growth_vs_senescence_summary.md
"""
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

GROWTH_CSV = ec.DATA_DIR / "dendro_growth_by_site_year.csv"
OUT_A = ec.DATA_DIR / "growth_vs_gpp.csv"
OUT_B = ec.DATA_DIR / "growth_vs_senescence.csv"
OUT_MD = ec.ROOT / "Output" / "growth_vs_senescence_summary.md"
FIG_DIR = ec.FIGURE_DIR / "dendro_growth"
FIG_DIR.mkdir(parents=True, exist_ok=True)
ec.require(GROWTH_CSV, hint="Run 14_dendrometer_download.py and 29_dendrometer_growth.py first.")

TARGET = 'EOS90'
CAL_PRE = 60
MIN_YEARS, MIN_OBS = 3, 12
MIN_SITES_CLUSTER = 5      # with fewer sites, clustered standard errors are unreliable: robust (HC1) ones are used
PREDICTORS = {'growth_pre': 'growth by the solstice', 'rate_pre': 'growth rate before the solstice',
              'growth_annual': 'annual growth', 'frac_pre': 'share of growth done by the solstice',
              'doy_g90': 'date growth stops', 'doy_g50': 'date half of the growth is done',
              'doy_rate_max': 'date of fastest growth', 'rate_max': 'fastest growth rate',
              'gpp_pre': f'GPP, {CAL_PRE} days before the solstice',
              't_spring': f'air temperature, {CAL_PRE} days before the solstice'}
SRC_ORDER = ['GCC', 'NDVI_tower', 'NDVI', 'NIRv']
NAME = {'GCC': 'PhenoCam', 'NDVI_tower': 'tower NDVI', 'NDVI': 'satellite NDVI', 'NIRv': 'satellite NIRv',
        'stacked': 'all sources stacked'}

growth = pd.read_csv(GROWTH_CSV)
flux = ec.load_flux_daily()
lk = ec.FluxLookup(flux, ['GPP', 'TA'])
sol = growth['year'].map(ec.solstice_doy)
growth['gpp_pre'] = [lk.window_mean(s, y, 'GPP', so - CAL_PRE, so - 1) * CAL_PRE
                     for s, y, so in zip(growth['site_id'], growth['year'], sol)]
growth['t_spring'] = [lk.window_mean(s, y, 'TA', so - CAL_PRE, so - 1) for s, y, so in zip(growth['site_id'], growth['year'], sol)]
growth['gpp_annual'] = [lk.window_mean(s, y, 'GPP', 1, 365) * 365 for s, y in zip(growth['site_id'], growth['year'])]


TIMING = {'doy_g90': 'date growth stops', 'doy_g50': 'date half of the growth is done',
          'doy_rate_max': 'date of fastest growth', 'rate_max': 'fastest growth rate'}
MIN_OBS_TIMING, MIN_YEARS_TIMING = 8, 2    # timing needs daily or monthly readings, which exist for very few
                                           # site-years: sites with two years are admitted, and the result is
                                           # shown for what it is


def wfit(d, y, x, unit='site_id', min_obs=None, min_years=None):
    """Within-`unit` slope of y on x (per +1 within-unit SD of x), SE clustered by site."""
    min_obs = MIN_OBS if min_obs is None else min_obs
    min_years = MIN_YEARS if min_years is None else min_years
    d = d[list(dict.fromkeys(['site_id', unit, y, x]))].dropna()
    d = d[d.groupby(unit)[unit].transform('size') >= min_years]
    if len(d) < min_obs:
        return None
    dm = d[[y, x]] - d.groupby(unit)[[y, x]].transform('mean')
    if not dm[x].std() > 0:
        return None
    X = sm.add_constant((dm[x] / dm[x].std()).to_numpy(float))
    n_sites = d['site_id'].nunique()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        f = sm.OLS(dm[y].to_numpy(float), X).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(d['site_id'])[0]}) \
            if n_sites >= MIN_SITES_CLUSTER else sm.OLS(dm[y].to_numpy(float), X).fit(cov_type='HC1')
    return {'beta': float(f.params[1]), 'se': float(f.bse[1]), 'p': float(f.pvalues[1]),
            'r_within': float(np.corrcoef(dm[x], dm[y])[0, 1]), 'n_obs': len(d), 'n_sites': n_sites,
            'sites': ', '.join(f"{s} ({n})" for s, n in d['site_id'].value_counts().items())}


# ---------------------------------------------------------------- A. growth vs GPP
a_rows = []
for y, x, lab in (('growth_annual', 'gpp_annual', 'annual growth ~ annual GPP'),
                  ('growth_pre', 'gpp_pre', 'growth by the solstice ~ pre-solstice GPP'),
                  ('rate_pre', 'gpp_pre', 'growth rate before the solstice ~ pre-solstice GPP'),
                  ('growth_pre', 't_spring', 'growth by the solstice ~ spring temperature')):
    r = wfit(growth, y, x)
    if r:
        a_rows.append({'relation': lab, **r})
A = pd.DataFrame(a_rows)
A.to_csv(OUT_A, index=False)
pd.set_option('display.width', 230)
pd.set_option('display.max_colwidth', 70)
print(f"Growth site-years with flux data: {int(growth['gpp_pre'].notna().sum())} of {len(growth)} "
      f"({growth.loc[growth['gpp_pre'].notna(), 'site_id'].nunique()} sites)")
print(f"\nA. Growth against GPP (within sites; r and slope in growth units per +1 SD) -> '{OUT_A}'")
print(A.round(3).to_string(index=False))

# ---------------------------------------------------------------- B. growth vs EOS90
pheno = ec.load_phenology()[['site_id', 'year', 'vi_index', TARGET]]
m = pheno.merge(growth, on=['site_id', 'year'])
m['unit'] = m['site_id'] + '|' + m['vi_index']
sources = [s for s in SRC_ORDER if s in set(m['vi_index'])]
b_rows = []
for src in sources + ['stacked']:
    d, unit = (m, 'unit') if src == 'stacked' else (m[m['vi_index'] == src], 'site_id')
    for x, lab in PREDICTORS.items():
        r = wfit(d, TARGET, x, unit, MIN_OBS_TIMING, MIN_YEARS_TIMING) if x in TIMING else wfit(d, TARGET, x, unit)
        if r:
            b_rows.append({'eos_source': src, 'predictor': x, 'label': lab, **r})
B = pd.DataFrame(b_rows)
B.to_csv(OUT_B, index=False)
print(f"\nB. {TARGET} against growth (days per +1 within-site SD) -> '{OUT_B}'")
if len(B):
    print(B.drop(columns=['label']).round(3).to_string(index=False))
else:
    print("   too few site-years with both growth and EOS90.")

# ---------------------------------------------------------------- figures
colors = dict(zip(sorted(growth['site_id'].unique()), plt.get_cmap('tab10').colors))
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
for ax, (y, x, yl, xl) in zip(axes, (('growth_annual', 'gpp_annual', 'annual growth anomaly', 'annual GPP anomaly (gC m-2)'),
                                     ('growth_pre', 'gpp_pre', 'growth by the solstice, anomaly',
                                      f'GPP anomaly, {CAL_PRE} days before the solstice (gC m-2)'))):
    d = growth[['site_id', y, x]].dropna()
    d = d[d.groupby('site_id')['site_id'].transform('size') >= MIN_YEARS]
    dm = d[[y, x]] - d.groupby('site_id')[[y, x]].transform('mean')
    for s in d['site_id'].unique():
        k = (d['site_id'] == s).to_numpy()
        ax.scatter(dm[x][k], dm[y][k], s=18, alpha=0.75, color=colors[s], label=f"{s} ({k.sum()})")
    if len(dm) > 2:
        b, a0 = np.polyfit(dm[x], dm[y], 1)
        xs = np.linspace(dm[x].min(), dm[x].max(), 20)
        ax.plot(xs, a0 + b * xs, color='#c53030', lw=1.3)
        ax.set_title(f"r = {np.corrcoef(dm[x], dm[y])[0, 1]:.2f}, n = {len(dm)}", fontsize=9)
    ax.axhline(0, color='0.85', lw=0.8)
    ax.axvline(0, color='0.85', lw=0.8)
    ax.set_xlabel(xl)
    ax.set_ylabel(f"{yl} (1 = a normal year's growth)")
    ax.legend(fontsize=7)
fig.suptitle("Stem growth against GPP, within sites (year minus site mean)", fontsize=10)
fig.tight_layout()
fig.savefig(FIG_DIR / "growth_vs_gpp.png", dpi=150)
plt.close(fig)

show = [p for p in ('growth_pre', 'growth_annual', 'gpp_pre') if p in set(B.get('predictor', []))]
if show:
    fig, axes = plt.subplots(1, len(show), figsize=(5 * len(show), 4.4), sharey=True, squeeze=False)
    for ax, x in zip(axes[0], show):
        d = m[['site_id', 'unit', 'vi_index', TARGET, x]].dropna()
        d = d[d.groupby('unit')['unit'].transform('size') >= MIN_YEARS]
        dm = d[[TARGET, x]] - d.groupby('unit')[[TARGET, x]].transform('mean')
        for s in d['site_id'].unique():
            k = (d['site_id'] == s).to_numpy()
            ax.scatter(dm[x][k], dm[TARGET][k], s=14, alpha=0.6, color=colors[s], label=f"{s} ({k.sum()})")
        b, a0 = np.polyfit(dm[x], dm[TARGET], 1)
        xs = np.linspace(dm[x].min(), dm[x].max(), 20)
        ax.plot(xs, a0 + b * xs, color='#c53030', lw=1.3)
        ax.axhline(0, color='0.85', lw=0.8)
        ax.axvline(0, color='0.85', lw=0.8)
        ax.set_title(f"{PREDICTORS[x]}\nr = {np.corrcoef(dm[x], dm[TARGET])[0, 1]:.2f}, n = {len(dm)}", fontsize=9)
        ax.set_xlabel('anomaly (year minus site mean)')
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel(f'{TARGET} anomaly (days), all EOS sources')
    fig.suptitle(f"Onset of senescence against stem growth, within sites (all EOS sources stacked)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "growth_vs_eos90.png", dpi=150)
    plt.close(fig)

# ---------------------------------------------------------------- C. carbon budget: total sink and what is left
def site_growth_scale():
    """Mean annual growth in gC m-2 for the sites where the stand's biomass growth is known."""
    out = {}
    hf = RAW_DENDRO / "harvard_forest" / "hf069-16-ems-annual.csv"
    if hf.exists():          # aboveground woody increment of surviving trees, MgC ha-1 yr-1 (2023-24 are provisional)
        a = pd.read_csv(hf)
        out['US-Ha1'] = (float(a.loc[a['year'].between(1999, 2022), 'agwi'].mean()) * 100.0,
                         'aboveground woody increment (HF069)')
    zo = list((RAW_DENDRO / "zoebelboden_lter").glob("LTER_EU_AT_003_ZOEBELBODEN_TREEGROWTH*.txt"))
    if zo:                   # stem + branch + foliage growth, kg dry mass m-2 per month, per plot and species
        z = pd.read_csv(zo[0], sep='\t', encoding='latin1', decimal=',')
        z['year'] = z['TIME'].str[:4].astype(int)
        per_plot_year = z.groupby(['STATION_CODE', 'year'])['VALUE'].sum()            # all species and organs of a plot
        out['AT-Zoe'] = (float(per_plot_year.groupby('year').mean().mean()) * 1000.0 * CARBON_FRACTION,
                         f'stem + branch + foliage growth, dry mass x {CARBON_FRACTION:g} (LTER Zoebelboden)')
    return out


RAW_DENDRO = ec.ROOT / "data_raw" / "dendro"
DAILY_CSV = ec.DATA_DIR / "dendro_growth_daily.csv"
OUT_C = ec.DATA_DIR / "carbon_budget_by_site_year.csv"
OUT_C_TEST = ec.DATA_DIR / "carbon_budget_vs_senescence.csv"
CARBON_FRACTION = 0.5
BUDGET_FROM, BUDGET_TO = 90, 305         # the season over which growth is counted in step 29
TERMS = {'gpp': 'GPP', 'ra': 'autotrophic respiration', 'growth': 'growth', 'sink_total': 'total sink (growth + respiration)',
         'npp': 'NPP (GPP - respiration)', 'residual': 'GPP - total sink'}
# the same with the respiration of the fitted model (growth + maintenance respiration, step 26)
MODEL_TERMS = {'rg_model': 'growth respiration (model)', 'rm_model': 'maintenance respiration (model)',
               'ra_model': 'autotrophic respiration (model)', 'sink_total_model': 'total sink (growth + model respiration)',
               'residual_model': 'GPP - total sink (model respiration)'}

scale = site_growth_scale()
flux['Ra'] = flux['GPP'] - flux['NPPd'].fillna(flux['NPP']) if 'NPPd' in flux.columns else np.nan   # (1 - CUE) x GPP
HAS_MODEL = 'Rm_model' in flux.columns and flux['Rm_model'].notna().any()
lkb = ec.FluxLookup(flux, ['GPP', 'Ra', 'Rg_model', 'Rm_model'])
gd = pd.read_csv(DAILY_CSV) if DAILY_CSV.exists() else pd.DataFrame()
c_rows = []
for (site, year), g in (gd[gd['site_id'].isin(scale)].groupby(['site_id', 'year']) if len(gd) else []):
    cum = g.sort_values('doy')['growth'].to_numpy(float) * scale[site][0]        # gC m-2, cumulative
    so = ec.solstice_doy(year)
    r = {'site_id': site, 'year': int(year)}
    for w, (a, b) in {'pre': (BUDGET_FROM, so - 1), 'post': (so, BUDGET_TO), 'season': (BUDGET_FROM, BUDGET_TO)}.items():
        n = b - a + 1
        gpp, ra = lkb.window_mean(site, year, 'GPP', a, b) * n, lkb.window_mean(site, year, 'Ra', a, b) * n
        gr = float(cum[b - 1] - cum[a - 2]) if a > 1 else float(cum[b - 1])
        r.update({f'gpp_{w}': gpp, f'ra_{w}': ra, f'growth_{w}': gr, f'sink_total_{w}': gr + ra, f'npp_{w}': gpp - ra,
                  f'residual_{w}': gpp - ra - gr})
        if HAS_MODEL:
            rg, rm = lkb.window_mean(site, year, 'Rg_model', a, b) * n, lkb.window_mean(site, year, 'Rm_model', a, b) * n
            r.update({f'rg_model_{w}': rg, f'rm_model_{w}': rm, f'ra_model_{w}': rg + rm,
                      f'sink_total_model_{w}': gr + rg + rm, f'residual_model_{w}': gpp - gr - rg - rm})
    c_rows.append(r)
C = pd.DataFrame(c_rows)
if len(C):
    C = C.dropna(subset=['gpp_season', 'ra_season'])
    C['growth_share_of_npp'] = C['growth_season'] / C['npp_season']
    if HAS_MODEL:
        TERMS = {**TERMS, **MODEL_TERMS}
C.to_csv(OUT_C, index=False)
ct_rows = []
if len(C):
    print(f"\nC. Carbon budget, DOY {BUDGET_FROM}-{BUDGET_TO} (gC m-2) -> '{OUT_C}'")
    print(C.groupby('site_id')[[f'{t}_season' for t in TERMS] + ['growth_share_of_npp']].mean().round(1).to_string())
    print("   mean before / after the solstice:")
    print(C.groupby('site_id')[[f'{t}_{w}' for w in ('pre', 'post') for t in TERMS]].mean().round(0).to_string())
    mc = pheno.merge(C, on=['site_id', 'year'])
    mc['unit'] = mc['site_id'] + '|' + mc['vi_index']
    for src in sources + ['stacked']:
        d, unit = (mc, 'unit') if src == 'stacked' else (mc[mc['vi_index'] == src], 'site_id')
        for w in ('pre', 'post'):
            for t, lab in TERMS.items():
                r = wfit(d, TARGET, f'{t}_{w}', unit)
                if r:
                    ct_rows.append({'eos_source': src, 'window': w, 'term': t, 'label': lab, **r})
CT = pd.DataFrame(ct_rows)
CT.to_csv(OUT_C_TEST, index=False)
if len(CT):
    print(f"\n   {TARGET} against the budget terms (days per +1 within-site SD) -> '{OUT_C_TEST}'")
    print(CT[CT['eos_source'].isin(['stacked', 'GCC'])].drop(columns=['label', 'sites']).round(2).to_string(index=False))

if len(C) and len(gd):                    # mean seasonal course of the budget per site
    sites_c = sorted(C['site_id'].unique())
    fig, axes = plt.subplots(1, len(sites_c), figsize=(6 * len(sites_c), 4.3), squeeze=False)
    for ax, site in zip(axes[0], sites_c):
        yrs = C.loc[C['site_id'] == site, 'year']
        f = flux[(flux['site_id'] == site) & flux['year'].isin(yrs)].groupby('doy')[['GPP', 'Ra']].mean()
        gg = gd[(gd['site_id'] == site) & gd['year'].isin(yrs)].sort_values(['year', 'doy'])
        gg = gg.assign(rate=gg.groupby('year')['growth'].diff().fillna(0) * scale[site][0]).groupby('doy')['rate'].mean()
        sm_ = lambda s: s.reindex(range(1, 366)).rolling(21, center=True, min_periods=7).mean()
        gpp, ra, gr = sm_(f['GPP']), sm_(f['Ra']), sm_(gg)
        ax.plot(gpp.index, gpp, color='#2f855a', label='GPP')
        ax.plot(ra.index, ra, color='#c05621', label='autotrophic respiration')
        ax.plot(gr.index, gr, color='#2b6cb0', label='growth')
        ax.plot(gr.index, ra + gr, color='k', lw=1.4, label='total sink (growth + respiration)')
        ax.fill_between(gpp.index, ra + gr, gpp, color='0.85', label='GPP - total sink')
        ax.axvline(172, color='gray', ls='--', lw=0.8)
        ax.set_title(f"{site}: mean of {len(yrs)} years", fontsize=9)
        ax.set_xlabel('day of year (dashed: summer solstice)')
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel('gC m-2 d-1 (21-day mean)')
    fig.suptitle("Carbon budget through the season: source, total sink and what is left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "carbon_budget.png", dpi=150)
    plt.close(fig)

# ---------------------------------------------------------------- written summary
cov =growth.groupby('site_id').agg(dataset=('dataset', 'first'), years=('year', 'size'), first=('year', 'min'),
                                    last=('year', 'max'), units=('n_units', 'median'), timing=('timing_ok', 'sum'),
                                    with_flux=('gpp_pre', 'count'), frac_pre=('frac_pre', 'mean'))
lines = ["# Stem growth against GPP and against the onset of senescence", "",
         "Generated by `Script/45_growth_vs_senescence.py`. Growth from dendrometers (step 29); 1 = a normal year's "
         "growth of the site. Within-site models. A case study on few sites: standard errors are robust but not "
         "clustered, and the stacked rows count each year once per EOS source, so their p-values are too small.", "",
         "## Growth data", "",
         "| site | dataset | years | units (median) | years with flux data | years with timing | share of growth done by the solstice |",
         "|---|---|---|---|---|---|---|"]
for s, r in cov.iterrows():
    lines.append(f"| {s} | {r['dataset']} | {int(r['years'])} ({int(r['first'])}-{int(r['last'])}) | {r['units']:.0f} | "
                 f"{int(r['with_flux'])} | {int(r['timing'])} | {100 * r['frac_pre']:.0f}% |")
lines += ["", "## A. Is growth coupled to GPP?", "", "| relation | within-site r | p | site-years | sites (years) |", "|---|---|---|---|---|"]
for _, r in A.iterrows():
    lines.append(f"| {r['relation']} | {r['r_within']:+.2f} | {r['p']:.3f} | {int(r['n_obs'])} | {r['sites']} |")
lines += ["", f"## B. Does growth predict {TARGET}? (days per +1 within-site SD)", "",
          "| EOS source | predictor | slope [95% CI] | p | site-years | sites (years) |", "|---|---|---|---|---|---|"]
for _, r in B.iterrows():
    lines.append(f"| {NAME[r['eos_source']]} | {r['label']} | {r['beta']:+.1f} [{r['beta'] - 1.96 * r['se']:+.1f}, "
                 f"{r['beta'] + 1.96 * r['se']:+.1f}] | {r['p']:.3f} | {int(r['n_obs'])} | {r['sites']} |")
if len(C):
    lines += ["", f"## C. Carbon budget: total sink, NPP and what is left (DOY {BUDGET_FROM}-{BUDGET_TO}, gC m-2)", "",
              "- respiration = autotrophic respiration = (1 - CUE) x GPP, with the daily CUE of step 26",
              "- growth = measured biomass growth, the dendrometer curve scaled to the stand's mean annual growth: "
              + "; ".join(f"{s}: {v:.0f} gC m-2 per year, {src}" for s, (v, src) in scale.items()),
              "- **total sink = growth + respiration**",
              "- **NPP = GPP - respiration**",
              "- **GPP - total sink = NPP - growth**: carbon fixed but used neither for respiration nor for the measured "
              "growth. It goes to what the dendrometers do not see (roots; at US-Ha1 also leaves) and to reserves.", "",
              "| site | years | window | " + " | ".join(TERMS.values()) + " |", "|---|---|---|" + "---|" * len(TERMS)]
    for site, g in C.groupby('site_id'):
        for w, wl in (('season', 'whole season'), ('pre', 'before the solstice'), ('post', 'after the solstice')):
            v = g[[f'{t}_{w}' for t in TERMS]].mean()
            lines.append(f"| {site} | {len(g)} | {wl} | " + " | ".join(f"{x:.0f}" for x in v) + " |")
    if HAS_MODEL:
        lines += ["", "Model respiration (step 26): growth respiration = gR x CUE x GPP; maintenance respiration = "
                  "temperature response x biomass built since 1 January. The maintenance term does not contain the "
                  "day's GPP, so a total sink built with it is not GPP times a factor."]
    lines += ["", "Growth as a share of NPP over the season: "
              + ", ".join(f"{s} {100 * v:.0f}%" for s, v in C.groupby('site_id')['growth_share_of_npp'].mean().items()) + ".", ""]
if len(CT):
    lines += [f"**{TARGET} against the budget terms** (days per +1 within-site SD)", "",
              "| EOS source | window | term | slope [95% CI] | p | site-years | sites (years) |", "|---|---|---|---|---|---|---|"]
    for _, r in CT.iterrows():
        lines.append(f"| {NAME[r['eos_source']]} | {r['window']}-solstice | {r['label']} | {r['beta']:+.1f} "
                     f"[{r['beta'] - 1.96 * r['se']:+.1f}, {r['beta'] + 1.96 * r['se']:+.1f}] | {r['p']:.3f} | {int(r['n_obs'])} | {r['sites']} |")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print(f"\nSummary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")
