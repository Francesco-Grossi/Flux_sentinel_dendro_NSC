"""
PIPELINE STEP 21 - How well do the different EOS sources agree, and which
one should be the outcome of the EOS <-> carbon analyses (steps 12-17)?

Sources (every vi_index present after QC, plus the flux-derived dates):
    NDVI / NIRv         satellite, HLS            (step 5)
    NDVI_tower          tower broadband NDVI      (step 19)
    GCC                 PhenoCam greenness        (step 20)
    GPP                 flux-derived, EOS10 only  (step 18) - NOT independent of
                        the GPP predictors; shown for reference only

The EOS <-> carbon analyses use year-to-year differences WITHIN a site, so
agreement is judged on within-site anomalies (each source's site mean
removed), over sites with >= MIN_YEARS_PER_SITE shared years:
    r_within     correlation of the anomalies       <- the number that matters
    rmse_within  RMSE of the anomalies (days)
    bias         mean (source_b - source_a) in days: a constant offset between
                 sources (e.g. leaf colour lagging GPP) is harmless
Also reported per source: site-years, sites, and the typical within-site
standard deviation of EOS (very large values indicate noise, not biology).

Input : phenology tables of steps 5 / 19 / 20, data/gpp_derived_eos_by_site_year.csv
Output: data/eos_source_summary.csv
        data/eos_source_agreement.csv           (all pairs x EOS90/50/10 x leaf habit)
        figure/eos_source_agreement/r_within_<level>.png
        figure/eos_source_agreement/anomaly_scatter_vs_<reference>.png
"""
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

LEVELS = ['EOS90', 'EOS50', 'EOS10']
MIN_YEARS_PER_SITE = 3
MIN_PAIR_OBS, MIN_PAIR_SITES = 15, 3
GPP_EOS_CSV = ec.DATA_DIR / "gpp_derived_eos_by_site_year.csv"
OUT_SUMMARY = ec.DATA_DIR / "eos_source_summary.csv"
OUT_AGREE = ec.DATA_DIR / "eos_source_agreement.csv"
FIG_DIR = ec.FIGURE_DIR / "eos_source_agreement"
FIG_DIR.mkdir(parents=True, exist_ok=True)
# independent ground references, best first; the first one present is used for the scatter figure
REFERENCES = ['GCC', 'NDVI_tower']
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous',
              'ENF': 'evergreen', 'EBF': 'evergreen'}   # everything else: grass/shrub

pheno = ec.load_phenology()[['site_id', 'year', 'vi_index'] + LEVELS].rename(columns={'vi_index': 'source'})
if os.path.exists(GPP_EOS_CSV):
    g = pd.read_csv(GPP_EOS_CSV).rename(columns={f'GPP_{l}': l for l in LEVELS})
    pheno = pd.concat([pheno, g.reindex(columns=['site_id', 'year'] + LEVELS).assign(source='GPP')],
                      ignore_index=True)
else:
    print(f"NOTE: '{GPP_EOS_CSV}' not found - run 18_vi_eos_vs_gpp_eos.py to include flux-derived EOS.")

igbp = pd.read_csv(ec.DATA_DIR / "fluxnet_daily_all_vars.csv", usecols=['site_id', 'igbp']).drop_duplicates('site_id')
pheno = pheno.merge(igbp, on='site_id', how='left')
pheno['leaf_habit'] = pheno['igbp'].map(LEAF_HABIT).fillna('grass/shrub')
sources = [s for s in ['NDVI', 'NIRv', 'NDVI_tower', 'GCC', 'GPP'] if s in set(pheno['source'])]
print(f"EOS sources: {sources}")

# ------------------------------------------------------------ per-source summary
summary = []
for src in sources:
    d = pheno[pheno['source'] == src]
    for lvl in LEVELS:
        x = d.dropna(subset=[lvl])
        if x.empty:
            continue   # e.g. GPP-derived EOS90 / EOS50, which are not computed
        sd = x.groupby('site_id')[lvl].std()
        summary.append({'source': src, 'level': lvl, 'n_site_years': len(x), 'n_sites': x['site_id'].nunique(),
                        'first_year': int(x['year'].min()) if len(x) else np.nan,
                        'last_year': int(x['year'].max()) if len(x) else np.nan,
                        'median_doy': x[lvl].median(), 'median_within_site_sd_days': sd.median()})
summary = pd.DataFrame(summary)
summary.to_csv(OUT_SUMMARY, index=False)
print(f"\nPer-source summary -> '{OUT_SUMMARY}'")
print(summary.round(1).to_string(index=False))


# ------------------------------------------------------------ pairwise agreement
def paired_anomalies(a, b, lvl, habit=None):
    """Site-years where both sources have `lvl`, with within-site anomalies."""
    da = pheno[pheno['source'] == a]
    if habit:
        da = da[da['leaf_habit'] == habit]
    m = da[['site_id', 'year', lvl]].merge(pheno.loc[pheno['source'] == b, ['site_id', 'year', lvl]],
                                           on=['site_id', 'year'], suffixes=('_a', '_b')).dropna()
    m = m[m.groupby('site_id')['year'].transform('size') >= MIN_YEARS_PER_SITE].copy()
    for s in ('a', 'b'):
        m[f'anom_{s}'] = m[f'{lvl}_{s}'] - m.groupby('site_id')[f'{lvl}_{s}'].transform('mean')
    return m


rows = []
for a, b in itertools.combinations(sources, 2):
    for lvl in LEVELS:
        for habit in [None, 'deciduous', 'evergreen', 'grass/shrub']:
            m = paired_anomalies(a, b, lvl, habit)
            if len(m) < MIN_PAIR_OBS or m['site_id'].nunique() < MIN_PAIR_SITES:
                continue
            ok = m['anom_a'].std() > 0 and m['anom_b'].std() > 0
            rows.append({'source_a': a, 'source_b': b, 'level': lvl, 'leaf_habit': habit or 'all',
                         'n_site_years': len(m), 'n_sites': m['site_id'].nunique(),
                         'r_within': float(np.corrcoef(m['anom_a'], m['anom_b'])[0, 1]) if ok else np.nan,
                         'rmse_within_days': float(np.sqrt(((m['anom_b'] - m['anom_a']) ** 2).mean())),
                         'bias_days': float((m[f'{lvl}_b'] - m[f'{lvl}_a']).mean()),
                         'r_pooled': float(np.corrcoef(m[f'{lvl}_a'], m[f'{lvl}_b'])[0, 1])})
agree = pd.DataFrame(rows)
agree.to_csv(OUT_AGREE, index=False)
if agree.empty:
    raise SystemExit("No pair of sources shares enough site-years.")
print(f"\nPairwise agreement -> '{OUT_AGREE}'")
print(agree[agree['leaf_habit'] == 'all'].drop(columns='leaf_habit').round(2).to_string(index=False))

# ------------------------------------------------------------ figures
for lvl in LEVELS:
    d = agree[(agree['level'] == lvl) & (agree['leaf_habit'] == 'all')]
    mat = pd.DataFrame(np.nan, index=sources, columns=sources)
    n = pd.DataFrame(0, index=sources, columns=sources)
    for _, r in d.iterrows():
        mat.loc[r['source_a'], r['source_b']] = mat.loc[r['source_b'], r['source_a']] = r['r_within']
        n.loc[r['source_a'], r['source_b']] = n.loc[r['source_b'], r['source_a']] = r['n_site_years']
    fig, ax = plt.subplots(figsize=(1.1 * len(sources) + 2, 1.0 * len(sources) + 1.5))
    im = ax.imshow(mat.to_numpy(float), vmin=-1, vmax=1, cmap='RdBu')
    ax.set_xticks(range(len(sources)), sources, rotation=45, ha='right')
    ax.set_yticks(range(len(sources)), sources)
    for i, j in itertools.product(range(len(sources)), repeat=2):
        if np.isfinite(mat.iloc[i, j]):
            ax.text(j, i, f"{mat.iloc[i, j]:.2f}\nn={n.iloc[i, j]}", ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, label='within-site r of EOS anomalies')
    ax.set_title(f"{lvl}: year-to-year agreement between EOS sources")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"r_within_{lvl}.png", dpi=150)
    plt.close(fig)

ref = next((r for r in REFERENCES if r in sources), None)
if ref:
    others = [s for s in sources if s != ref]
    fig, axes = plt.subplots(len(LEVELS), len(others), figsize=(3 * len(others), 3 * len(LEVELS)), squeeze=False)
    for i, lvl in enumerate(LEVELS):
        for j, src in enumerate(others):
            ax, m = axes[i, j], paired_anomalies(ref, src, lvl)
            ax.axhline(0, color='0.8', lw=0.8)
            ax.axvline(0, color='0.8', lw=0.8)
            if len(m) >= MIN_PAIR_OBS:
                ax.scatter(m['anom_a'], m['anom_b'], s=10, alpha=0.6)
                lim = float(np.nanpercentile(np.abs(m[['anom_a', 'anom_b']].to_numpy()), 98)) + 2
                ax.plot([-lim, lim], [-lim, lim], color='k', lw=0.8, ls='--')
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
                ax.set_title(f"{src} {lvl}  r={np.corrcoef(m['anom_a'], m['anom_b'])[0, 1]:.2f}  n={len(m)}",
                             fontsize=9)
            else:
                ax.set_title(f"{src} {lvl}  (too few shared years)", fontsize=9)
            if i == len(LEVELS) - 1:
                ax.set_xlabel(f"{ref} anomaly (days)")
            if j == 0:
                ax.set_ylabel("source anomaly (days)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"anomaly_scatter_vs_{ref}.png", dpi=150)
    plt.close(fig)

# ------------------------------------------------------------ which source to use
print("\nAgreement of each source with the independent references (within-site r, all sites):")
for r in [x for x in REFERENCES if x in sources]:
    d = agree[(agree['leaf_habit'] == 'all') & ((agree['source_a'] == r) | (agree['source_b'] == r))].copy()
    d['other'] = np.where(d['source_a'] == r, d['source_b'], d['source_a'])
    tab = d.pivot(index='other', columns='level', values='r_within').reindex(columns=LEVELS)
    print(f"\n  reference = {r}")
    print(tab.round(2).to_string())
print(f"\nFigures -> '{FIG_DIR}'")
