"""
PIPELINE STEP 43 - How solid is the pre-solstice GPP effect on EOS90?

Step 41 finds that a year with high GPP before the summer solstice has an
earlier EOS90 at the same site. This step checks three ways in which that
result could be fragile. The model is the within-site model of step 41 with
one predictor (year minus site mean, sites with >= 3 years, SE clustered by
site; slope in days per +1 within-site SD), for two versions of the predictor:
    sos   cumulative GPP from the year's leaf-out (SOS10) to the solstice
    cal   cumulative GPP over the CAL_DAYS days before the solstice
          (does not depend on leaf-out)

1. Leave one site out. The model is refitted with each site removed in turn.
   If one site drives the result, the slope changes sign or loses most of
   its size when that site is left out.

2. QC thresholds. The phenology QC of pheno_fit.py is re-applied to the stored
   curve fits with other thresholds: minimum fit R2 (R2_GRID) and the maximum
   distance of a year's EOS90 from the site median (DEV_GRID; None = rule
   off). Transition dates are recomputed from the stored fit parameters, so
   no curve is refitted. The effect should keep its sign and rough size.

3. Power. From the standard error of each fit: the smallest effect detectable
   with 80% power at the present sample size, and the number of site-years
   needed to detect TARGET_EFFECT days per SD (assuming the standard error
   shrinks with 1/sqrt(n) and the same number of years per site). Also for
   the site-years PhenoCam shares with each satellite index (step 42).

Input : phenology tables (steps 23, 24, 25), daily flux
Output: data/presolstice_loso.csv                 slope with each site left out
        data/presolstice_loso_summary.csv
        data/presolstice_qc_sensitivity.csv
        data/presolstice_power.csv
        figure/presolstice_robustness/loso.png, qc_sensitivity.png
        Output/presolstice_robustness_summary.md
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
import eos_common as ec
import pheno_fit as pf

OUT_LOSO = ec.DATA_DIR / "presolstice_loso.csv"
OUT_LOSO_SUM = ec.DATA_DIR / "presolstice_loso_summary.csv"
OUT_QC = ec.DATA_DIR / "presolstice_qc_sensitivity.csv"
OUT_POWER = ec.DATA_DIR / "presolstice_power.csv"
OUT_MD = ec.ROOT / "Output" / "presolstice_robustness_summary.md"
FIG_DIR = ec.FIGURE_DIR / "presolstice_robustness"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_MD.parent.mkdir(parents=True, exist_ok=True)

TARGET = 'EOS90'
CAL_DAYS = 60
PREDICTORS = {'sos': 'GPP from leaf-out to the solstice', 'cal': f'GPP in the {CAL_DAYS} days before the solstice'}
R2_GRID = [0.7, 0.8, 0.9]
DEV_GRID = [30, 60, 90, None]
MIN_YEARS, MIN_OBS, MIN_SITES = 3, 30, 5
TARGET_EFFECT, POWER, ALPHA = 1.5, 0.80, 0.05
Z_SUM = norm.ppf(1 - ALPHA / 2) + norm.ppf(POWER)          # 1.96 + 0.84
SRC_ORDER = ['GCC', 'NDVI_tower', 'NDVI', 'NIRv']
NAME = {'GCC': 'PhenoCam', 'NDVI_tower': 'tower NDVI', 'NDVI': 'satellite NDVI', 'NIRv': 'satellite NIRv'}
COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#805ad5', 'NDVI_tower': '#dd6b20', 'GCC': '#2f855a'}

flux = ec.load_flux_daily()
lk = ec.FluxLookup(flux, ['GPP'])


def add_predictors(p):
    """Adds the two pre-solstice GPP sums (gC m-2) to a phenology table."""
    p = p.copy()
    sol = p['year'].map(ec.solstice_doy).astype(float)
    p['cal'] = [lk.window_mean(s, y, 'GPP', so - CAL_DAYS, so - 1) * CAL_DAYS
                for s, y, so in zip(p['site_id'], p['year'], sol)]
    p['sos'] = [lk.window_mean(s, y, 'GPP', a, so) * ec.FluxLookup.n_days(a, so)
                for s, y, a, so in zip(p['site_id'], p['year'], p['leaf_out_10'], sol)]
    return p


def fit(d, x):
    wf = ec.within_fit(d, TARGET, [x], MIN_YEARS, MIN_OBS, MIN_SITES)
    if wf is None:
        return None
    b, se, p = wf[0][x]
    return {'beta': b, 'se': se, 'p': p, 'n_obs': wf[1], 'n_sites': wf[2]}


pheno = add_predictors(ec.load_phenology())
sources = [s for s in SRC_ORDER if s in set(pheno['vi_index'])]

# ---------------------------------------------------------------- 1. leave one site out
loso_rows, sum_rows = [], []
for src in sources:
    d = pheno[pheno['vi_index'] == src]
    for x in PREDICTORS:
        full = fit(d, x)
        if full is None:
            continue
        dd = d[['site_id', TARGET, x]].dropna()
        dd = dd[dd.groupby('site_id')['site_id'].transform('size') >= MIN_YEARS]
        part = []
        for site in sorted(dd['site_id'].unique()):
            r = fit(dd[dd['site_id'] != site], x)
            if r:
                part.append({'eos_source': src, 'predictor': x, 'site_left_out': site,
                             'n_years_of_site': int((dd['site_id'] == site).sum()), **r,
                             'change_vs_all_sites': r['beta'] - full['beta']})
        if not part:
            continue
        loso_rows += part
        b = np.array([r['beta'] for r in part])
        top = max(part, key=lambda r: abs(r['change_vs_all_sites']))
        sum_rows.append({'eos_source': src, 'predictor': x, **{f'all_{k}': v for k, v in full.items()},
                         'loso_min_beta': b.min(), 'loso_max_beta': b.max(),
                         'share_negative': float((b < 0).mean()),
                         'share_p05': float(np.mean([r['p'] < 0.05 for r in part])),
                         'most_influential_site': top['site_left_out'],
                         'beta_without_it': top['beta'], 'p_without_it': top['p']})
loso, loso_sum = pd.DataFrame(loso_rows), pd.DataFrame(sum_rows)
loso.to_csv(OUT_LOSO, index=False)
loso_sum.to_csv(OUT_LOSO_SUM, index=False)
pd.set_option('display.width', 250)
print(f"1. Leave one site out -> '{OUT_LOSO_SUM}'")
print(loso_sum.round(3).to_string(index=False))

# ---------------------------------------------------------------- 2. QC thresholds
PARAMS = ['vmin', 'vmax', 'S', 'greenup_kinetic_i', 'A', 'senescence_kinetic_i']
DATES = ['leaf_out_10', 'leaf_out_50', 'leaf_out_90', 'peak_doy', 'EOS90', 'EOS50', 'EOS10']
raw = pd.concat([pd.read_csv(p) for p in [ec.PHENOLOGY_CSV] + ec.EXTRA_PHENOLOGY_CSVS if os.path.exists(p)],
                ignore_index=True)
raw = raw[(raw['method'] == 'double_logistic') & raw[PARAMS].notna().all(axis=1)].reset_index(drop=True)
stored_pass = raw['qc_pass'].astype(bool).to_numpy()
dates = pd.DataFrame([pf.transition_dates(r) for r in raw[PARAMS].to_numpy(float)]).reindex(columns=DATES + ['amplitude'])
raw[DATES + ['amplitude']] = dates.to_numpy()
seq = raw[DATES].to_numpy(float)
raw['base_ok'] = (raw['n_obs_after_peak'] >= pf.MIN_OBS_AFTER_PEAK) \
    & (raw['amplitude'] >= pf.MIN_AMP_TO_RMSE * raw['rmse']) \
    & ~np.isnan(seq).any(axis=1) & (np.diff(seq, axis=1) > 0).all(axis=1) & (raw['corr'] >= ec.MIN_FIT_CORR)
raw = add_predictors(raw)


def apply_qc(min_r2, max_dev):
    ok = raw['base_ok'] & (raw['r2'] >= min_r2)
    if max_dev is not None:
        med = raw[TARGET].where(ok).groupby([raw['site_id'], raw['vi_index']]).transform('median')
        cnt = ok.groupby([raw['site_id'], raw['vi_index']]).transform('sum')
        ok = ok & ~((cnt >= pf.MIN_YEARS_FOR_SITE_CHECK) & ((raw[TARGET] - med).abs() > max_dev))
    return ok


check = apply_qc(pf.MIN_R2, pf.MAX_EOS90_DEV_DAYS).to_numpy()
print(f"\n2. QC thresholds. Re-applying the default QC reproduces the stored flag for "
      f"{100 * (check == stored_pass).mean():.1f}% of {len(raw)} fits.")
qc_rows = []
for min_r2 in R2_GRID:
    for max_dev in DEV_GRID:
        ok = apply_qc(min_r2, max_dev)
        for src in sources:
            d = raw[ok & (raw['vi_index'] == src)]
            for x in PREDICTORS:
                r = fit(d, x)
                if r:
                    qc_rows.append({'min_r2': min_r2, 'max_eos90_dev_days': max_dev if max_dev is not None else np.inf,
                                    'default': min_r2 == pf.MIN_R2 and max_dev == pf.MAX_EOS90_DEV_DAYS,
                                    'eos_source': src, 'predictor': x, 'n_fits_passing': int(len(d)), **r})
qc = pd.DataFrame(qc_rows)
qc.to_csv(OUT_QC, index=False)
print(f"   -> '{OUT_QC}'")
print(qc.pivot_table(index=['predictor', 'min_r2', 'max_eos90_dev_days'], columns='eos_source', values='beta')
      .reindex(columns=sources).round(2).to_string())

# ---------------------------------------------------------------- 3. power
pw_rows = []


def power_row(sample, src, x, r):
    need = r['n_obs'] * (Z_SUM * r['se'] / TARGET_EFFECT) ** 2
    pw_rows.append({'sample': sample, 'eos_source': src, 'predictor': x, **r,
                    'years_per_site': r['n_obs'] / r['n_sites'],
                    'detectable_effect_days_per_sd': Z_SUM * r['se'],
                    'power_for_target_effect': float(norm.cdf(TARGET_EFFECT / r['se'] - norm.ppf(1 - ALPHA / 2))),
                    'site_years_needed': need, 'sites_needed': need / (r['n_obs'] / r['n_sites']),
                    'enough': r['n_obs'] >= need})


for src in sources:
    for x in PREDICTORS:
        r = fit(pheno[pheno['vi_index'] == src], x)
        if r:
            power_row('own site-years', src, x, r)
if 'GCC' in sources:
    ref = pheno[pheno['vi_index'] == 'GCC'][['site_id', 'year', TARGET]].dropna()
    for sat in ('NDVI', 'NIRv'):
        if sat not in sources:
            continue
        keys = ref.merge(pheno[pheno['vi_index'] == sat][['site_id', 'year', TARGET]].dropna(),
                         on=['site_id', 'year'])[['site_id', 'year']]
        for src in ('GCC', sat):
            for x in PREDICTORS:
                r = ec.within_fit(pheno[pheno['vi_index'] == src].merge(keys, on=['site_id', 'year']),
                                  TARGET, [x], MIN_YEARS, 20, MIN_SITES)
                if r:
                    b, se, p = r[0][x]
                    power_row(f'shared PhenoCam + {sat}', src, x, {'beta': b, 'se': se, 'p': p, 'n_obs': r[1], 'n_sites': r[2]})
power = pd.DataFrame(pw_rows)
power.to_csv(OUT_POWER, index=False)
print(f"\n3. Power to detect {TARGET_EFFECT} days per SD ({100 * POWER:.0f}% power, two-sided alpha {ALPHA}) -> '{OUT_POWER}'")
print(power[['sample', 'eos_source', 'predictor', 'n_obs', 'n_sites', 'beta', 'se', 'detectable_effect_days_per_sd',
             'power_for_target_effect', 'site_years_needed', 'sites_needed']].round(2).to_string(index=False))

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, len(PREDICTORS), figsize=(6 * len(PREDICTORS), 4.2), sharey=True, squeeze=False)
for ax, (x, lab) in zip(axes[0], PREDICTORS.items()):
    for i, src in enumerate(sources):
        s = loso[(loso['eos_source'] == src) & (loso['predictor'] == x)]
        a = loso_sum[(loso_sum['eos_source'] == src) & (loso_sum['predictor'] == x)]
        if s.empty or a.empty:
            continue
        jitter = np.linspace(-0.25, 0.25, len(s))
        ax.scatter(i + jitter, np.sort(s['beta']), s=9, alpha=0.6, color=COLORS.get(src, 'gray'))
        ax.errorbar(i + 0.38, a['all_beta'].iloc[0], yerr=1.96 * a['all_se'].iloc[0], fmt='D', ms=6, capsize=3, color='k')
    ax.axhline(0, color='k', lw=0.7)
    ax.set_xticks(range(len(sources)), [NAME[s] for s in sources], fontsize=8)
    ax.set_title(lab, fontsize=9)
axes[0][0].set_ylabel(f'{TARGET} shift (days per +1 within-site SD)')
fig.suptitle("Leave one site out: each dot is the slope with one site removed (diamond: all sites, 95% CI)", fontsize=10)
fig.tight_layout()
fig.savefig(FIG_DIR / "loso.png", dpi=150)
plt.close(fig)

if len(qc):
    labels = [f"R2>={a:g}\n{'no rule' if not np.isfinite(b) else f'{b:.0f} d'}" for a in R2_GRID
              for b in [np.inf if v is None else v for v in DEV_GRID]]
    keys = [(a, np.inf if v is None else v) for a in R2_GRID for v in DEV_GRID]
    fig, axes = plt.subplots(1, len(PREDICTORS), figsize=(7.5 * len(PREDICTORS), 4.4), sharey=True, squeeze=False)
    for ax, (x, lab) in zip(axes[0], PREDICTORS.items()):
        for i, src in enumerate(sources):
            s = qc[(qc['eos_source'] == src) & (qc['predictor'] == x)].set_index(['min_r2', 'max_eos90_dev_days']).reindex(keys)
            xx = np.arange(len(keys)) + (i - (len(sources) - 1) / 2) * 0.17
            ax.errorbar(xx, s['beta'], yerr=1.96 * s['se'], fmt='o', ms=4, lw=1, capsize=2, color=COLORS.get(src, 'gray'),
                        label=NAME[src])
        ax.axhline(0, color='k', lw=0.7)
        ax.axvspan(keys.index((pf.MIN_R2, pf.MAX_EOS90_DEV_DAYS)) - 0.5, keys.index((pf.MIN_R2, pf.MAX_EOS90_DEV_DAYS)) + 0.5,
                   color='0.9', zorder=0)
        ax.set_xticks(range(len(keys)), labels, fontsize=7)
        ax.set_title(lab, fontsize=9)
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel(f'{TARGET} shift (days per +1 within-site SD)')
    fig.suptitle("QC thresholds: minimum fit R2 and maximum distance of EOS90 from the site median (grey: default)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "qc_sensitivity.png", dpi=150)
    plt.close(fig)

# ---------------------------------------------------------------- written summary
lines = ["# Pre-solstice GPP and EOS90 - robustness", "",
         "Generated by `Script/43_presolstice_robustness.py`. Within-site model, one predictor; days of EOS90 shift "
         "per +1 within-site SD of pre-solstice GPP.", "",
         "## 1. Leave one site out", "",
         "| EOS source | predictor | all sites (p) | range with one site removed | negative in | p < 0.05 in | most influential site | slope without it (p) |",
         "|---|---|---|---|---|---|---|---|"]
for _, r in loso_sum.iterrows():
    lines.append(f"| {NAME[r['eos_source']]} | {r['predictor']} | {r['all_beta']:+.1f} ({r['all_p']:.3f}) | "
                 f"{r['loso_min_beta']:+.1f} to {r['loso_max_beta']:+.1f} | {100 * r['share_negative']:.0f}% | "
                 f"{100 * r['share_p05']:.0f}% | {r['most_influential_site']} | {r['beta_without_it']:+.1f} ({r['p_without_it']:.3f}) |")
lines += ["", "## 2. QC thresholds", ""]
for x, lab in PREDICTORS.items():
    q = qc[qc['predictor'] == x]
    if q.empty:
        continue
    lines += [f"**{lab}** - slope (site-years)", "", "| min R2 | max EOS90 distance | " + " | ".join(NAME[s] for s in sources) + " |",
              "|---|---|" + "---|" * len(sources)]
    for (a, b), g in q.groupby(['min_r2', 'max_eos90_dev_days']):
        g = g.set_index('eos_source')
        cells = [f"{g.at[s, 'beta']:+.1f}{ec.stars(g.at[s, 'p'])} ({int(g.at[s, 'n_obs'])})" if s in g.index else '-' for s in sources]
        mark = ' (default)' if bool(g['default'].iloc[0]) else ''
        lines.append(f"| {a:g}{mark} | {'no rule' if not np.isfinite(b) else f'{b:.0f} d'} | " + " | ".join(cells) + " |")
    lines.append("")
lines += [f"## 3. Power to detect {TARGET_EFFECT} days per SD ({100 * POWER:.0f}% power, two-sided alpha {ALPHA})", "",
          "| sample | EOS source | predictor | site-years (sites) | slope | smallest detectable effect | power now | site-years needed | sites needed |",
          "|---|---|---|---|---|---|---|---|---|"]
for _, r in power.iterrows():
    lines.append(f"| {r['sample']} | {NAME[r['eos_source']]} | {r['predictor']} | {int(r['n_obs'])} ({int(r['n_sites'])}) | "
                 f"{r['beta']:+.1f} | {r['detectable_effect_days_per_sd']:.1f} | {100 * r['power_for_target_effect']:.0f}% | "
                 f"{r['site_years_needed']:.0f} | {r['sites_needed']:.0f} |")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print(f"\nSummary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")
