"""
PIPELINE STEP 36 - Pathway (piecewise structural) analysis: how do
environment -> carbon uptake -> EOS chains work, and is carbon uptake a
mediator of environmental effects on EOS? (Notion D.1-3b: "construct
different concept structures ... like Fig. 1 in Zani et al. 2020 or Fig. 4
in Lu et al. 2022".)

Method: within-site (site-demeaned, then standardized) OLS for every
equation in EQUATIONS below = fixed-effects equivalent of the random-
intercept models used elsewhere, fast enough to bootstrap. Path
coefficients are assembled into a matrix B; total effects are
(I - B)^-1 - I, indirect = total - direct. 95% CIs come from a CLUSTER
(site) bootstrap.

*** EQUATIONS is the ONLY thing to edit to mirror the paper diagrams. ***
The default is a generic structure, NOT a transcription of Zani Fig.1 /
Lu Fig.4 (those figures are not available to the script):

    SOS, env(pre)  -> C_pre            carbon uptake, [SOS, SOL]
    C_pre, env(post) -> C_post         carbon uptake, [SOL, EOS site mean]
    C_pre, C_post, SOS, env -> EOS     direct effects on the target

Run once per carbon variable C in {GPP, NPP (if available)}, per vi_index,
per target in PRIMARY_TARGETS (fixed-anchor windows from step 27).

Input : data/eos_window_predictors_fixed_anchor.csv   (step 27)
Output: data/eos_path_coefficients.csv   (each equation's coefficients)
        data/eos_path_effects.csv        (direct / indirect / total, with CI)
        figure/eos_path_analysis/<vi>_<target>_<carbon>.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

N_BOOT = int(os.environ.get("PHENO_NBOOT", 500))
SEED = 42
PRE, POST = 'SOS_to_SOL', 'SOL_to_EOS10'

EQUATIONS = {                                   # dependent -> predictors
    'C_pre':  ['SOS', 'TA_pre', 'SW_pre', 'P_pre'],
    'C_post': ['C_pre', 'TA_post', 'SW_post', 'P_post'],
    'EOS':    ['C_pre', 'C_post', 'SOS', 'TA_pre', 'SW_pre', 'P_pre', 'TA_post', 'SW_post', 'P_post'],
}
OUT_COEF = ec.DATA_DIR / "eos_path_coefficients.csv"
OUT_EFF = ec.DATA_DIR / "eos_path_effects.csv"
FIG_DIR = ec.FIGURE_DIR / "eos_path_analysis"
FIG_DIR.mkdir(parents=True, exist_ok=True)
ec.require(ec.WINDOW_FIXED_CSV, hint="Run 27_window_predictors.py first.")

w = pd.read_csv(ec.WINDOW_FIXED_CSV)


def column_map(cv, target):
    m = {'SOS': 'SOS', 'EOS': target,
         'C_pre': f'{cv}_cum__{PRE}', 'C_post': f'{cv}_cum__{POST}',
         'TA_pre': f'TA_mean__{PRE}', 'TA_post': f'TA_mean__{POST}',
         'SW_pre': f'SW_mean__{PRE}', 'SW_post': f'SW_mean__{POST}',
         'P_pre': f'P_sum__{PRE}', 'P_post': f'P_sum__{POST}'}
    return {k: v for k, v in m.items() if v in w.columns}


def within_z(d, cols):
    dm = d[cols] - d.groupby('site_id')[cols].transform('mean')
    return dm / dm.std()


def solve(Z, eqs, order):
    """Fit each equation; return (B matrix, {eq: coef dict})."""
    idx = {v: i for i, v in enumerate(order)}
    B, coefs = np.zeros((len(order), len(order))), {}
    for dep, preds in eqs.items():
        X = Z[:, [idx[p] for p in preds]]
        beta, *_ = np.linalg.lstsq(X, Z[:, idx[dep]], rcond=None)
        for p, b in zip(preds, beta):
            B[idx[dep], idx[p]] = b
        coefs[dep] = dict(zip(preds, beta))
    return B, coefs


def effects_on_eos(B, order):
    n, e = len(order), order.index('EOS')
    total = np.linalg.inv(np.eye(n) - B) - np.eye(n)
    return {v: (B[e, i], total[e, i] - B[e, i], total[e, i]) for i, v in enumerate(order) if v != 'EOS'}


rng = np.random.default_rng(SEED)
carbon_vars = sorted({c.split('_cum__')[0] for c in w.columns if '_cum__' in c} & {'GPP', 'NPP', 'NPPd'})
coef_rows, eff_rows = [], []

for cv in carbon_vars:
    for vi in sorted(w['vi_index'].unique()):
        for target in ec.PRIMARY_TARGETS:
            cmap = column_map(cv, target)
            eqs = {dep: [p for p in preds if p in cmap] for dep, preds in EQUATIONS.items() if dep in cmap}
            eqs = {dep: p for dep, p in eqs.items() if p}
            order = [v for v in cmap if v not in eqs] + list(eqs)          # exogenous first, then equation order
            if 'EOS' not in eqs or len(order) < 3:
                continue
            d = w[w['vi_index'] == vi][['site_id'] + list(cmap.values())].dropna()
            d = d.groupby('site_id').filter(lambda g: len(g) >= 3)
            n_obs, n_sites = len(d), d['site_id'].nunique()
            if n_obs < ec.MIN_OBS or n_sites < ec.MIN_SITES:
                continue
            inv = {v: k for k, v in cmap.items()}
            zdf = within_z(d, list(cmap.values())).rename(columns=inv)[order]
            Z = zdf.to_numpy()
            B, coefs = solve(Z, eqs, order)
            point = effects_on_eos(B, order)

            site_rows = [np.where(d['site_id'].to_numpy() == s)[0] for s in d['site_id'].unique()]
            boot = {v: [] for v in point}
            for _ in range(N_BOOT):
                pick = rng.integers(0, len(site_rows), len(site_rows))
                Zb = Z[np.concatenate([site_rows[i] for i in pick])]
                try:
                    Bb, _ = solve(Zb, eqs, order)
                    for v, vals in effects_on_eos(Bb, order).items():
                        boot[v].append(vals)
                except np.linalg.LinAlgError:
                    continue

            for dep, cd in coefs.items():
                for p, b in cd.items():
                    coef_rows.append({'vi_index': vi, 'target': target, 'carbon': cv, 'equation': dep,
                                      'predictor': p, 'std_coef': float(b), 'n_obs': n_obs, 'n_sites': n_sites})
            for v, (dr, ind, tot) in point.items():
                arr = np.array(boot[v]) if boot[v] else np.full((1, 3), np.nan)
                for j, kind, est in ((0, 'direct', dr), (1, 'indirect', ind), (2, 'total', tot)):
                    lo, hi = np.nanpercentile(arr[:, j], [2.5, 97.5])
                    eff_rows.append({'vi_index': vi, 'target': target, 'carbon': cv, 'source': v,
                                     'effect': kind, 'estimate': float(est), 'ci_lo': lo, 'ci_hi': hi,
                                     'n_obs': n_obs, 'n_sites': n_sites, 'n_boot': N_BOOT})

            # ---- figure: direct vs total effect on EOS per source variable
            e = pd.DataFrame([r for r in eff_rows if (r['vi_index'], r['target'], r['carbon']) == (vi, target, cv)
                              and r['effect'] in ('direct', 'total')])
            srcs = list(dict.fromkeys(e['source']))
            fig, ax = plt.subplots(figsize=(6.5, 0.45 * len(srcs) + 1.8))
            for k, (kind, color) in enumerate((('direct', '#a0aec0'), ('total', '#2b6cb0'))):
                s = e[e['effect'] == kind].set_index('source').reindex(srcs)
                y = np.arange(len(srcs)) + (k - 0.5) * 0.38
                ax.barh(y, s['estimate'], 0.38, color=color, label=kind,
                        xerr=[s['estimate'] - s['ci_lo'], s['ci_hi'] - s['estimate']], capsize=2)
            ax.set_yticks(np.arange(len(srcs)))
            ax.set_yticklabels([cmap[s_] if s_ in cmap else s_ for s_ in srcs], fontsize=7)
            ax.axvline(0, color='k', lw=0.6)
            ax.set_xlabel('standardized effect on EOS (within-site SD units), 95% site-bootstrap CI')
            ax.set_title(f'{vi} / {target} / {cv}  (n={n_obs}, sites={n_sites})', fontsize=9)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(FIG_DIR / f'{vi}_{target}_{cv}.png', dpi=140)
            plt.close(fig)
            print(f"[{vi}/{target}/{cv}] n={n_obs}, sites={n_sites} - "
                  + ", ".join(f"{v}: tot={t:+.2f}" for v, (_, _, t) in point.items() if v in ('C_pre', 'C_post')))

pd.DataFrame(coef_rows).to_csv(OUT_COEF, index=False)
pd.DataFrame(eff_rows).to_csv(OUT_EFF, index=False)
print(f"\nWritten: {OUT_COEF.name}, {OUT_EFF.name}; figures -> '{FIG_DIR}'")
