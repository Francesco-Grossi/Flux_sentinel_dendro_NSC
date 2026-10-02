"""
PIPELINE STEP 38 - One results table for the central question, split by leaf
habit: does carbon uptake (source: GPP; sink: NPP, CUE) before / after the
summer solstice shift the end of season (EOS)?

For every combination of
    leaf habit   all | deciduous (DBF, DNF, MF) | evergreen (ENF, EBF) | grass/shrub
    EOS source   every vi_index after QC (NDVI, NIRv, NDVI_tower, GCC)
    target       EOS10, EOS50
    predictor    see PREDICTORS below
two WITHIN-SITE models are fitted (year minus site mean, sites with >= 3 years,
SE clustered by site):
    raw     EOS ~ predictor
    adj_T   EOS ~ predictor + mean air temperature before and after the solstice
            (the predictor's effect net of the main climatic driver)
beta = days of EOS shift per +1 within-site SD of the predictor (negative =
earlier EOS).

Predictor windows:
    pre    leaf-out (SOS) -> solstice                      (from step 27, fixed anchors)
    post   solstice -> the site's mean EOS10               (from step 27, fixed anchors)
    cal    the CAL_DAYS days before the solstice           (computed here; does not depend
           on the fitted leaf-out date, so it is immune to noise in SOS)

Reading the result:
  * p_fdr is the Benjamini-Hochberg corrected p-value over ALL cells of a
    model type - with this many tests, uncorrected p < 0.05 is expected by
    chance in ~5% of cells.
  * The EOS sources are measured on the same site-years, so they are not
    independent replicates; they are compared by SIGN AGREEMENT
    (n_sources_same_sign), not pooled. A robust effect has the same sign in
    all sources and survives FDR in the independent ones (GCC, NDVI_tower).

Input : data/eos_window_predictors_fixed_anchor.csv (step 27), daily flux + NPP file
Output: data/eos_results_by_leaf_habit.csv          every cell
        data/eos_results_consistency.csv            per habit x target x predictor: sign agreement
        figure/eos_results_by_leaf_habit/<target>_<model>.png
        Output/eos_results_summary.md               plain-text summary of the above
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

OUT_CELLS = ec.DATA_DIR / "eos_results_by_leaf_habit.csv"
OUT_CONS = ec.DATA_DIR / "eos_results_consistency.csv"
OUT_MD = ec.ROOT / "Output" / "eos_results_summary.md"
FIG_DIR = ec.FIGURE_DIR / "eos_results_by_leaf_habit"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_MD.parent.mkdir(parents=True, exist_ok=True)
ec.require(ec.WINDOW_FIXED_CSV, hint="Run 27_window_predictors.py first.")

TARGETS = ['EOS10', 'EOS50']
PRE, POST, CAL_DAYS = 'SOS_to_SOL', 'SOL_to_EOS10', 60
T_COVARS = [f'TA_mean__{PRE}', f'TA_mean__{POST}']
HABITS = ['all', 'deciduous', 'evergreen', 'grass/shrub']
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous', 'ENF': 'evergreen', 'EBF': 'evergreen'}
INDEPENDENT_SOURCES = ['GCC', 'NDVI_tower']        # ground references, independent of the satellite
MIN_OBS, MIN_SITES = 30, 5                         # per cell (subsets are small; below this, skip)
# label -> column; '{w}' is filled with the window name
PREDICTORS = {
    'GPP cumulative (source)':      'GPP_cum__{w}',
    'NPP cumulative, annual CUE':   'NPP_cum__{w}',
    'NPP cumulative, daily CUE':    'NPPd_cum__{w}',
    'GPP rate (mean)':              'GPP_mean__{w}',
    'CUE seasonal mean':            'CUEd_mean__{w}',
}
CAL_VARS = {'GPP rate (mean)': 'GPP', 'NPP rate, daily CUE': 'NPPd', 'CUE seasonal mean': 'CUEd'}

w = pd.read_csv(ec.WINDOW_FIXED_CSV)
w['leaf_habit'] = w['igbp'].map(LEAF_HABIT).fillna('grass/shrub')

# calendar window before the solstice, independent of the fitted SOS
flux = ec.load_flux_daily()
cal_vars = {lab: v for lab, v in CAL_VARS.items() if v in flux.columns}
lk = ec.FluxLookup(flux, list(cal_vars.values()))
for lab, v in cal_vars.items():
    w[f'{v}_mean__cal'] = [lk.window_mean(s, y, v, ec.solstice_doy(y) - CAL_DAYS, ec.solstice_doy(y) - 1)
                           for s, y in zip(w['site_id'], w['year'])]

specs = []                                          # (window label, predictor label, column)
for win_lab, win in (('pre', PRE), ('post', POST)):
    for lab, col in PREDICTORS.items():
        if col.format(w=win) in w.columns:
            specs.append((win_lab, lab, col.format(w=win)))
for lab, v in cal_vars.items():
    specs.append((f'cal{CAL_DAYS}', lab, f'{v}_mean__cal'))
t_covars = [c for c in T_COVARS if c in w.columns]
sources = sorted(w['vi_index'].unique())
print(f"EOS sources: {sources}\nPredictors: {len(specs)} | temperature covariates: {t_covars}")

# ------------------------------------------------------------------ all cells
rows = []
for habit in HABITS:
    wh = w if habit == 'all' else w[w['leaf_habit'] == habit]
    for src in sources:
        ws = wh[wh['vi_index'] == src]
        for target in TARGETS:
            for win_lab, lab, col in specs:
                for model, covars in (('raw', []), ('adj_T', t_covars)):
                    d = ws[[target, 'site_id', col] + covars].dropna()
                    if len(d) < MIN_OBS or d['site_id'].nunique() < MIN_SITES:
                        continue
                    wf = ec.within_fit(d, target, [col] + covars, ec.MIN_YEARS_WITHIN, MIN_OBS, MIN_SITES)
                    if wf is None:
                        continue
                    b, se, p = wf[0][col]
                    rows.append({'leaf_habit': habit, 'eos_source': src, 'target': target, 'window': win_lab,
                                 'predictor': lab, 'model': model, 'beta_days_per_sd': b,
                                 'std_err': se, 'p_value': p, 'n_obs': wf[1], 'n_sites': wf[2]})
cells = pd.DataFrame(rows)
if cells.empty:
    raise SystemExit("No cell had enough data.")


def bh_fdr(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    out = np.empty(len(p))
    out[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(out, 0, 1)


cells['p_fdr'] = cells.groupby('model')['p_value'].transform(bh_fdr)
cells.to_csv(OUT_CELLS, index=False)
n_raw = int((cells['p_value'] < 0.05).sum())
n_fdr = int((cells['p_fdr'] < 0.05).sum())
print(f"\n{len(cells)} cells -> '{OUT_CELLS}'")
print(f"  p < 0.05 uncorrected: {n_raw} ({100 * n_raw / len(cells):.1f}% of cells; ~5% expected by chance)")
print(f"  p < 0.05 after FDR  : {n_fdr}")

# ------------------------------------------------------------------ sign agreement across sources
key = ['leaf_habit', 'target', 'window', 'predictor', 'model']
cons = []
for k, g in cells.groupby(key):
    sign = np.sign(g['beta_days_per_sd'])
    major = 1.0 if (sign > 0).sum() >= (sign < 0).sum() else -1.0
    ind = g[g['eos_source'].isin(INDEPENDENT_SOURCES)]
    cons.append({**dict(zip(key, k)), 'n_sources': len(g), 'n_sources_same_sign': int((sign == major).sum()),
                 'direction': 'later EOS' if major > 0 else 'earlier EOS',
                 'median_beta': g['beta_days_per_sd'].median(),
                 'min_beta': g['beta_days_per_sd'].min(), 'max_beta': g['beta_days_per_sd'].max(),
                 'n_sources_p05': int((g['p_value'] < 0.05).sum()), 'n_sources_fdr05': int((g['p_fdr'] < 0.05).sum()),
                 'independent_sources_agree': bool(len(ind) == len(INDEPENDENT_SOURCES)
                                                   and (np.sign(ind['beta_days_per_sd']) == major).all()),
                 'sources': ', '.join(f"{s}:{b:+.1f}" for s, b in zip(g['eos_source'], g['beta_days_per_sd']))})
cons = pd.DataFrame(cons)
cons['robust'] = (cons['n_sources'] >= 3) & (cons['n_sources_same_sign'] == cons['n_sources']) \
    & (cons['n_sources_fdr05'] >= 1) & cons['independent_sources_agree']
cons.to_csv(OUT_CONS, index=False)
print(f"Sign agreement across EOS sources -> '{OUT_CONS}'")

# ------------------------------------------------------------------ figures (forest plots)
colors = {'NDVI': '#2b6cb0', 'NIRv': '#805ad5', 'NDVI_tower': '#dd6b20', 'GCC': '#2f855a'}
for target in TARGETS:
    for model in ('raw', 'adj_T'):
        sub = cells[(cells['target'] == target) & (cells['model'] == model)]
        if sub.empty:
            continue
        labels = [f"{wl}: {lab}" for wl, lab, _ in specs]
        fig, axes = plt.subplots(1, len(HABITS), figsize=(4.2 * len(HABITS), 0.42 * len(labels) + 1.8), sharey=True)
        for ax, habit in zip(axes, HABITS):
            sh = sub[sub['leaf_habit'] == habit]
            for si, src in enumerate(sources):
                ss = sh[sh['eos_source'] == src].assign(label=lambda d: d['window'] + ': ' + d['predictor'])
                ss = ss.set_index('label').reindex(labels)
                y = np.arange(len(labels)) + (si - (len(sources) - 1) / 2) * 0.18
                ax.errorbar(ss['beta_days_per_sd'], y, xerr=1.96 * ss['std_err'], fmt='o', ms=3.5, lw=1,
                            color=colors.get(src, 'gray'), label=src)
                sig = ss['p_fdr'] < 0.05
                ax.scatter(ss['beta_days_per_sd'][sig], y[sig.to_numpy()], s=60, facecolors='none',
                           edgecolors=colors.get(src, 'gray'))
            ax.axvline(0, color='k', lw=0.6)
            n = sh.groupby('eos_source')['n_obs'].max()
            ax.set_title(f"{habit}\n(n = {int(n.min()) if len(n) else 0}-{int(n.max()) if len(n) else 0} site-years)",
                         fontsize=9)
            ax.set_xlabel('EOS shift (days per +1 within-site SD)')
        axes[0].set_yticks(np.arange(len(labels)))
        axes[0].set_yticklabels(labels, fontsize=8)
        axes[0].invert_yaxis()
        axes[-1].legend(fontsize=7, loc='lower right')
        fig.suptitle(f"{target}, model '{model}' - points: beta with 95% CI; circled: p < 0.05 after FDR", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"{target}_{model}.png", dpi=150)
        plt.close(fig)

# ------------------------------------------------------------------ text summary
lines = ["# EOS vs carbon source / sink - results by leaf habit", "",
         f"Generated by `Script/38_eos_results_by_leaf_habit.py` from `{ec.WINDOW_FIXED_CSV.name}`.",
         "Within-site models: beta = days of EOS shift per +1 within-site SD of the predictor (negative = earlier EOS).", "",
         f"- EOS sources: {', '.join(sources)}",
         f"- Cells tested: {len(cells)} ({n_raw} with uncorrected p < 0.05 = {100 * n_raw / len(cells):.1f}%; "
         f"about 5% expected by chance; {n_fdr} with p < 0.05 after FDR)", ""]
robust = cons[cons['robust']]
lines += ["## Robust effects",
          "Same sign in every EOS source (at least 3), both independent ground sources agree, "
          "and at least one source significant after FDR.", ""]
if robust.empty:
    lines += ["**None.**", ""]
else:
    for _, r in robust.sort_values(['leaf_habit', 'target', 'model']).iterrows():
        lines.append(f"- {r['leaf_habit']} / {r['target']} / {r['model']} / {r['window']}: {r['predictor']} -> "
                     f"{r['direction']} (median {r['median_beta']:+.1f} d/SD; {r['sources']})")
    lines.append("")
lines += ["## Core hypothesis: pre-solstice carbon uptake -> EOS10", "",
          "| leaf habit | model | predictor | window | sources with same sign | median beta | per source |",
          "|---|---|---|---|---|---|---|"]
core = cons[(cons['target'] == 'EOS10') & cons['window'].isin(['pre', f'cal{CAL_DAYS}'])]
for _, r in core.sort_values(['leaf_habit', 'model', 'predictor', 'window']).iterrows():
    lines.append(f"| {r['leaf_habit']} | {r['model']} | {r['predictor']} | {r['window']} | "
                 f"{r['n_sources_same_sign']}/{r['n_sources']} {r['direction']} | {r['median_beta']:+.1f} | {r['sources']} |")
lines += ["", "## Significant after FDR (any cell)", ""]
sig = cells[cells['p_fdr'] < 0.05].sort_values('p_fdr')
if sig.empty:
    lines.append("None.")
else:
    lines += ["| leaf habit | EOS source | target | model | window | predictor | beta | p_fdr | n |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in sig.iterrows():
        lines.append(f"| {r['leaf_habit']} | {r['eos_source']} | {r['target']} | {r['model']} | {r['window']} | "
                     f"{r['predictor']} | {r['beta_days_per_sd']:+.1f} | {r['p_fdr']:.3f} | {r['n_obs']} |")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print(f"Summary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")

# ------------------------------------------------------------------ console: the core table
pd.set_option('display.width', 250)
pd.set_option('display.max_colwidth', 80)
for habit in HABITS:
    c = core[(core['leaf_habit'] == habit) & (core['model'] == 'raw')]
    if len(c):
        print(f"\n[{habit}] pre-solstice predictors -> EOS10 (raw model), beta per EOS source:")
        print(c[['window', 'predictor', 'n_sources_same_sign', 'n_sources', 'median_beta', 'sources']]
              .round(1).to_string(index=False))
print(f"\nRobust effects: {len(robust)}")
