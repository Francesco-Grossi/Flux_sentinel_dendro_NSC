"""
PIPELINE STEP 42 - Do PhenoCam and satellite EOS tell the same story about
the split-GPP hypothesis, and if not, why?

Step 41 runs each EOS source on its own site-years, so PhenoCam (few sites,
from 2008) and satellite NDVI / NIRv (many sites, from 2013) are not
directly comparable. This step puts them on equal terms.

1. Same sample. The within-site model of step 41
       dEOS ~ dpre + dpost      and      dEOS ~ dtotal
   is fitted for each source on
       own       every site-year the source has                (as step 41)
       sites     only the sites that have PhenoCam, all years
       shared    only the site-years where BOTH sources have an EOS
   If PhenoCam and the satellite converge on 'shared', their different
   answers came from which sites and years each one covers; if not, from
   what each one measures.
   To make the predictors identical for both sources they are CALENDAR
   windows, not windows tied to a source's own leaf-out or EOS dates:
       pre    the PRE_DAYS days before the summer solstice
       post   the POST_DAYS days from the solstice (ends before senescence
              normally starts)
       total  pre + post

2. Formal test of the difference. On the shared site-years the difference
       EOS(satellite) - EOS(PhenoCam)
   is regressed on the same predictors (within sites). Everything the two
   sources share cancels; a non-zero slope means GPP shifts the two EOS
   measures differently.

3. Same event? GCC follows leaf colour (chlorophyll loss, top of canopy,
   oblique view); NDVI is seen from above over 1 km, includes the understory
   and stays high until leaves fall. So one source's EOS90 may correspond to
   the other's EOS50. Every PhenoCam level is correlated with every
   satellite level (within sites), with the mean lag in days.

Input : phenology tables (steps 23, 25), daily flux
Output: data/phenocam_vs_satellite_models.csv
        data/phenocam_vs_satellite_difference.csv
        data/phenocam_vs_satellite_levels.csv
        figure/phenocam_vs_satellite/models_<target>.png, levels_<source>.png, scatter_<source>.png
        Output/phenocam_vs_satellite_summary.md
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

OUT_MODELS = ec.DATA_DIR / "phenocam_vs_satellite_models.csv"
OUT_DIFF = ec.DATA_DIR / "phenocam_vs_satellite_difference.csv"
OUT_LEVELS = ec.DATA_DIR / "phenocam_vs_satellite_levels.csv"
OUT_MD = ec.ROOT / "Output" / "phenocam_vs_satellite_summary.md"
FIG_DIR = ec.FIGURE_DIR / "phenocam_vs_satellite"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_MD.parent.mkdir(parents=True, exist_ok=True)

REFERENCE = 'GCC'
SATELLITE = ['NDVI', 'NIRv']
LEVELS = ['EOS90', 'EOS50', 'EOS10']
TARGETS = ['EOS90', 'EOS10']
PRE_DAYS, POST_DAYS = 60, 45
MIN_YEARS, MIN_OBS, MIN_SITES = 3, 25, 5
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous', 'ENF': 'evergreen', 'EBF': 'evergreen'}
NAME = {'GCC': 'PhenoCam', 'NDVI': 'satellite NDVI', 'NIRv': 'satellite NIRv'}

pheno = ec.load_phenology()
if REFERENCE not in set(pheno['vi_index']):
    raise SystemExit("No PhenoCam phenology - run steps 13 and 25 first.")
flux = ec.load_flux_daily()
lk = ec.FluxLookup(flux, ['GPP'])
habit = flux.drop_duplicates('site_id').set_index('site_id')['igbp'].map(LEAF_HABIT).fillna('grass/shrub')

# identical calendar-window predictors for every source (gC m-2)
sy = pheno[['site_id', 'year']].drop_duplicates().reset_index(drop=True)
sol = sy['year'].map(ec.solstice_doy)
sy['pre'] = [lk.window_mean(s, y, 'GPP', so - PRE_DAYS, so - 1) * PRE_DAYS for s, y, so in zip(sy['site_id'], sy['year'], sol)]
sy['post'] = [lk.window_mean(s, y, 'GPP', so, so + POST_DAYS - 1) * POST_DAYS for s, y, so in zip(sy['site_id'], sy['year'], sol)]
sy['total'] = sy['pre'] + sy['post']

wide = {v: pheno[pheno['vi_index'] == v][['site_id', 'year'] + LEVELS].merge(sy, on=['site_id', 'year'])
        for v in [REFERENCE] + [s for s in SATELLITE if s in set(pheno['vi_index'])]}
satellites = [s for s in SATELLITE if s in wide]
ref_sites = set(wide[REFERENCE]['site_id'])


def fit_rows(d, target, tag):
    """Within-site joint and total model -> one result row (or None)."""
    j = ec.within_fit(d, target, ['pre', 'post'], MIN_YEARS, MIN_OBS, MIN_SITES)
    t = ec.within_fit(d, target, ['total'], MIN_YEARS, MIN_OBS, MIN_SITES)
    if j is None or t is None:
        return None
    (bp, sp, pp), (bq, sq, pq), (bt, st, pt) = j[0]['pre'], j[0]['post'], t[0]['total']
    return {**tag, 'n_obs': j[1], 'n_sites': j[2], 'b_pre': bp, 'se_pre': sp, 'p_pre': pp,
            'b_post': bq, 'se_post': sq, 'p_post': pq, 'b_total': bt, 'se_total': st, 'p_total': pt}


# ---------------------------------------------------------------- 1. same model, three samples
rows = []
for target in TARGETS:
    for sat in satellites:
        shared = wide[REFERENCE][['site_id', 'year', target]].merge(
            wide[sat][['site_id', 'year', target]], on=['site_id', 'year'], suffixes=('_ref', '_sat')).dropna()
        keys = shared[['site_id', 'year']]
        samples = {
            'own': {REFERENCE: wide[REFERENCE], sat: wide[sat]},
            'sites': {REFERENCE: wide[REFERENCE], sat: wide[sat][wide[sat]['site_id'].isin(ref_sites)]},
            'shared': {REFERENCE: wide[REFERENCE].merge(keys, on=['site_id', 'year']),
                       sat: wide[sat].merge(keys, on=['site_id', 'year'])},
        }
        for sample, by_src in samples.items():
            for src, d in by_src.items():
                if src == REFERENCE and sample == 'sites':
                    continue                                  # identical to 'own' for PhenoCam
                for group in ('all', 'deciduous'):
                    dg = d if group == 'all' else d[d['site_id'].map(habit) == 'deciduous']
                    r = fit_rows(dg, target, {'target': target, 'pair': f'{REFERENCE} vs {sat}', 'sample': sample,
                                              'group': group, 'eos_source': src})
                    if r:
                        rows.append(r)
models = pd.DataFrame(rows).drop_duplicates(['target', 'sample', 'group', 'eos_source', 'n_obs'])
models.to_csv(OUT_MODELS, index=False)
pd.set_option('display.width', 240)
print(f"Same model on three samples -> '{OUT_MODELS}'")
print(models[models['group'] == 'all'].drop(columns=['group', 'se_pre', 'se_post', 'se_total']).round(2).to_string(index=False))

# ---------------------------------------------------------------- 2. difference model on shared site-years
drow = []
for sat in satellites:
    for lvl in LEVELS:
        m = wide[REFERENCE][['site_id', 'year', lvl]].merge(wide[sat][['site_id', 'year', lvl, 'pre', 'post', 'total']],
                                                            on=['site_id', 'year'], suffixes=('_ref', '_sat')).dropna()
        m['diff'] = m[f'{lvl}_sat'] - m[f'{lvl}_ref']
        r = fit_rows(m, 'diff', {'level': lvl, 'difference': f'{sat} - {REFERENCE}'})
        if r:
            r['mean_difference_days'] = float(m['diff'].mean())
            drow.append(r)
diff = pd.DataFrame(drow)
diff.to_csv(OUT_DIFF, index=False)
print(f"\nDifference model (satellite EOS - PhenoCam EOS ~ GPP, within sites) -> '{OUT_DIFF}'")
if len(diff):
    print(diff.drop(columns=['se_pre', 'se_post', 'se_total']).round(2).to_string(index=False))

# ---------------------------------------------------------------- 3. which levels correspond?
lrow = []
for sat in satellites:
    for a in LEVELS:
        for b in LEVELS:
            m = wide[REFERENCE][['site_id', 'year', a]].rename(columns={a: 'ref'}).merge(
                wide[sat][['site_id', 'year', b]].rename(columns={b: 'sat'}), on=['site_id', 'year']).dropna()
            m = m[m.groupby('site_id')['year'].transform('size') >= MIN_YEARS]
            if len(m) < MIN_OBS:
                continue
            dm = m[['ref', 'sat']] - m.groupby('site_id')[['ref', 'sat']].transform('mean')
            lrow.append({'satellite': sat, 'phenocam_level': a, 'satellite_level': b, 'n_obs': len(m),
                         'n_sites': m['site_id'].nunique(), 'r_within': float(np.corrcoef(dm['ref'], dm['sat'])[0, 1]),
                         'lag_days_satellite_minus_phenocam': float((m['sat'] - m['ref']).mean()),
                         'rmse_within_days': float(np.sqrt(((dm['sat'] - dm['ref']) ** 2).mean()))})
levels = pd.DataFrame(lrow)
levels.to_csv(OUT_LEVELS, index=False)
print(f"\nWhich levels correspond? -> '{OUT_LEVELS}'")
for sat in satellites:
    lv = levels[levels['satellite'] == sat]
    if len(lv):
        print(f"\n  within-site r, PhenoCam level (rows) x {sat} level (columns):")
        print(lv.pivot(index='phenocam_level', columns='satellite_level', values='r_within').reindex(index=LEVELS, columns=LEVELS).round(2).to_string())
        print("  lag (days, satellite - PhenoCam):")
        print(lv.pivot(index='phenocam_level', columns='satellite_level', values='lag_days_satellite_minus_phenocam').reindex(index=LEVELS, columns=LEVELS).round(0).to_string())

# ---------------------------------------------------------------- figures
colors = {'GCC': '#2f855a', 'NDVI': '#2b6cb0', 'NIRv': '#805ad5'}
for target in TARGETS:
    sub = models[(models['target'] == target) & (models['group'] == 'all')]
    if sub.empty:
        continue
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, (sample, title) in zip(axes, (('own', "each source's own site-years"), ('sites', 'PhenoCam sites only'),
                                           ('shared', 'shared site-years only'))):
        ss = sub[sub['sample'] == sample]
        if sample == 'sites':                       # PhenoCam is its own reference here
            ss = pd.concat([sub[(sub['sample'] == 'own') & (sub['eos_source'] == REFERENCE)], ss])
        k = 0
        for (src, pair), r in ss.drop_duplicates(['eos_source', 'pair']).set_index(['eos_source', 'pair']).iterrows():
            x = np.arange(3) + (k - 1.5) * 0.17
            ax.errorbar(x, [r['b_pre'], r['b_post'], r['b_total']],
                        yerr=[1.96 * r['se_pre'], 1.96 * r['se_post'], 1.96 * r['se_total']], fmt='o', ms=5, capsize=2,
                        color=colors.get(src, 'gray'), alpha=1.0 if src != REFERENCE else 0.9,
                        label=f"{NAME[src]}{'' if src != REFERENCE else ' (' + pair.split(' vs ')[1] + ' rows)' if sample == 'shared' else ''} n={int(r['n_obs'])}")
            k += 1
        ax.axhline(0, color='k', lw=0.7)
        ax.set_xticks(range(3), [f'pre\n({PRE_DAYS} d before SOL)', f'post\n({POST_DAYS} d after SOL)', 'total'])
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=7)
    axes[0].set_ylabel(f'{target} shift (days per +1 within-site SD of GPP)')
    fig.suptitle(f"{target}: the same within-site model on three samples (bars: 95% CI)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"models_{target}.png", dpi=150)
    plt.close(fig)

for sat in satellites:
    lv = levels[levels['satellite'] == sat]
    if lv.empty:
        continue
    mat = lv.pivot(index='phenocam_level', columns='satellite_level', values='r_within').reindex(index=LEVELS, columns=LEVELS)
    lag = lv.pivot(index='phenocam_level', columns='satellite_level', values='lag_days_satellite_minus_phenocam').reindex(index=LEVELS, columns=LEVELS)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    im = ax.imshow(mat.to_numpy(float), vmin=0, vmax=1, cmap='viridis')
    for i in range(3):
        for j in range(3):
            if np.isfinite(mat.iloc[i, j]):
                ax.text(j, i, f"r={mat.iloc[i, j]:.2f}\n{lag.iloc[i, j]:+.0f} d", ha='center', va='center', color='w', fontsize=9)
    ax.set_xticks(range(3), LEVELS)
    ax.set_yticks(range(3), LEVELS)
    ax.set_xlabel(f'{NAME[sat]} level')
    ax.set_ylabel('PhenoCam level')
    fig.colorbar(im, ax=ax, label='within-site r')
    ax.set_title(f"Which stages correspond? (lag = {sat} - PhenoCam)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"levels_{sat}.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, lvl in zip(axes, LEVELS):
        m = wide[REFERENCE][['site_id', 'year', lvl]].merge(wide[sat][['site_id', 'year', lvl]], on=['site_id', 'year'],
                                                            suffixes=('_ref', '_sat')).dropna()
        m = m[m.groupby('site_id')['year'].transform('size') >= MIN_YEARS]
        ax.axhline(0, color='0.85', lw=0.8)
        ax.axvline(0, color='0.85', lw=0.8)
        if len(m) > 3:
            dm = m[[f'{lvl}_ref', f'{lvl}_sat']] - m.groupby('site_id')[[f'{lvl}_ref', f'{lvl}_sat']].transform('mean')
            dec = (m['site_id'].map(habit) == 'deciduous').to_numpy()
            ax.scatter(dm.iloc[:, 0][~dec], dm.iloc[:, 1][~dec], s=12, alpha=0.6, color='0.5', label='other')
            ax.scatter(dm.iloc[:, 0][dec], dm.iloc[:, 1][dec], s=14, alpha=0.7, color=colors[sat], label='deciduous')
            lim = float(np.nanpercentile(np.abs(dm.to_numpy()), 98)) + 2
            ax.plot([-lim, lim], [-lim, lim], 'k--', lw=0.8)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_title(f"{lvl}: r = {np.corrcoef(dm.iloc[:, 0], dm.iloc[:, 1])[0, 1]:.2f}, n = {len(m)}", fontsize=9)
        ax.set_xlabel('PhenoCam anomaly (days)')
    axes[0].set_ylabel(f'{NAME[sat]} anomaly (days)')
    axes[0].legend(fontsize=7)
    fig.suptitle(f"PhenoCam vs {NAME[sat]}: EOS anomalies on shared site-years (1:1 line dashed)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"scatter_{sat}.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- written summary
def cell(b, se, p):
    return f"{b:+.1f}{ec.stars(p)} [{b - 1.96 * se:+.1f}, {b + 1.96 * se:+.1f}]"


lines = ["# PhenoCam vs satellite EOS - the split-GPP result on equal terms", "",
         "Generated by `Script/42_compare_phenocam_satellite.py`. Within-site models; days of EOS shift per +1 "
         f"within-site SD of GPP [95% CI]. Predictors are calendar windows: pre = {PRE_DAYS} days before the solstice, "
         f"post = {POST_DAYS} days from the solstice. Stars: two-sided * p<0.05, ** p<0.01, *** p<0.001.", ""]
for target in TARGETS:
    lines += [f"## {target}: same model, three samples", "",
              "| sample | sites | EOS source | n (sites) | pre | post | total |", "|---|---|---|---|---|---|---|"]
    for _, r in models[models['target'] == target].sort_values(['pair', 'group', 'sample', 'eos_source']).iterrows():
        lines.append(f"| {r['sample']} ({r['pair']}) | {r['group']} | {NAME[r['eos_source']]} | {int(r['n_obs'])} ({int(r['n_sites'])}) | "
                     f"{cell(r['b_pre'], r['se_pre'], r['p_pre'])} | {cell(r['b_post'], r['se_post'], r['p_post'])} | "
                     f"{cell(r['b_total'], r['se_total'], r['p_total'])} |")
    lines.append("")
lines += ["## Does GPP shift the two EOS measures differently? (difference model, shared site-years)", "",
          "| difference | level | n (sites) | mean difference (days) | pre | post | total |", "|---|---|---|---|---|---|---|"]
for _, r in diff.iterrows():
    lines.append(f"| {r['difference']} | {r['level']} | {int(r['n_obs'])} ({int(r['n_sites'])}) | {r['mean_difference_days']:+.1f} | "
                 f"{cell(r['b_pre'], r['se_pre'], r['p_pre'])} | {cell(r['b_post'], r['se_post'], r['p_post'])} | "
                 f"{cell(r['b_total'], r['se_total'], r['p_total'])} |")
lines += ["", "## Which stages correspond? Within-site r (lag in days, satellite minus PhenoCam)", ""]
for sat in satellites:
    lv = levels[levels['satellite'] == sat]
    if lv.empty:
        continue
    lines += [f"**PhenoCam (rows) x {NAME[sat]} (columns)**", "", "| | " + " | ".join(LEVELS) + " |", "|---|---|---|---|"]
    for a in LEVELS:
        cells = []
        for b in LEVELS:
            x = lv[(lv['phenocam_level'] == a) & (lv['satellite_level'] == b)]
            cells.append(f"{x['r_within'].iloc[0]:.2f} ({x['lag_days_satellite_minus_phenocam'].iloc[0]:+.0f} d)" if len(x) else "-")
        lines.append(f"| {a} | " + " | ".join(cells) + " |")
    lines.append("")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print(f"\nSummary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")
