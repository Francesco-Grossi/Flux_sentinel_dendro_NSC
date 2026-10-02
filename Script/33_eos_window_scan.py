"""
PIPELINE STEP 33 - Which carbon window links to which EOS, and with what sign?

Covers Notion sections A (previous analysis) and D.1 (GPP negative
relationship with EOS, sink activities via NPP):

  * carbon variable  : GPP, NEP (NEE-based) and - if data/npp_luo2025_daily.csv
                       exists - NPP (sink activity, Luo et al. 2025)
  * metric           : cumulative (cGPP) vs mean/rate (GPPrate)
  * windows          : [SOS,SOL] [SOL,EOS90] [SOL,EOS10] [SOS,EOS90] [SOS,EOS50]
                       [SOS,EOS10] [EOS90,EOS50]
  * targets          : EOS10, EOS50, EOS90
  * anchor mode      : 'fixed' (site-mean EOS as window end; the correct test)
                       and 'year' (this year's own EOS; the earlier analysis).
                       Comparing the two shows how much of a "positive cGPP -
                       EOS" relationship is just season length.

Each cell = WITHIN-SITE slope of EOS on the predictor (year minus site mean,
SE clustered by site; beta = days of EOS shift per +1 within-site SD), plus
the site-random-intercept mixed-model slope for comparison (lme_ columns;
inflated by differences between sites) and pooled / within-site Pearson r.

Input : data/eos_window_predictors_{fixed,year}_anchor.csv   (step 27)
Output: data/eos_window_scan.csv
        figure/eos_window_scan/<anchor>_<vi>_<carbon>.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

FIG_DIR = ec.FIGURE_DIR / "eos_window_scan"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV = ec.DATA_DIR / "eos_window_scan.csv"
ec.require(ec.WINDOW_FIXED_CSV, ec.WINDOW_YEAR_CSV, hint="Run 27_window_predictors.py first.")

rows = []
for mode, path in (('fixed', ec.WINDOW_FIXED_CSV), ('year', ec.WINDOW_YEAR_CSV)):
    w = pd.read_csv(path)
    carbon_vars = sorted({c.split('_cum__')[0] for c in w.columns if '_cum__' in c})
    ratio_vars = [v for v in ec.RATIO_VARS if any(c.startswith(f'{v}_mean__') for c in w.columns)]
    for vi in sorted(w['vi_index'].unique()):
        wv = w[w['vi_index'] == vi]
        for target in ec.TARGETS:
            for cv in carbon_vars + ratio_vars:
                windows = [c.split('__')[1] for c in w.columns if c.startswith(f'{cv}_mean__')]
                for metric in (('mean',) if cv in ratio_vars else ('cum', 'mean')):
                    for win in windows:
                        res = ec.site_slope(wv, target, f'{cv}_{metric}__{win}')
                        if res:
                            rows.append({'anchor': mode, 'vi_index': vi, 'target': target,
                                         'carbon': cv, 'metric': metric, 'window': win, **res})

scan = pd.DataFrame(rows)
scan.to_csv(OUTPUT_CSV, index=False)
print(f"{len(scan)} window x metric x target cells -> '{OUTPUT_CSV}'")
if scan.empty:
    raise SystemExit("No cell had enough data (need >= "
                     f"{ec.MIN_OBS} obs / {ec.MIN_SITES} sites). Check step 27 output.")

both = scan.dropna(subset=['lme_beta_days_per_sd']) if 'lme_beta_days_per_sd' in scan else scan.iloc[:0]
if len(both):
    ratio = (both['lme_beta_days_per_sd'].abs() / both['beta_days_per_sd'].abs().clip(lower=1e-9)).median()
    print(f"Mixed model vs within-site: median |beta| {both['lme_beta_days_per_sd'].abs().median():.2f} vs "
          f"{both['beta_days_per_sd'].abs().median():.2f} d/SD (median ratio {ratio:.1f}); "
          f"p < 0.05 in {int((both['lme_p_value'] < 0.05).sum())} vs {int((both['p_value'] < 0.05).sum())} of {len(both)} cells; "
          f"same sign in {100 * (np.sign(both['lme_beta_days_per_sd']) == np.sign(both['beta_days_per_sd'])).mean():.0f}%.")

# ---- console summary: the "sign shift" of the Notion page, EOS10 & EOS50, GPP
pd.set_option('display.width', 200)
for mode in ('fixed', 'year'):
    for target in ('EOS10', 'EOS50'):
        sub = scan[(scan['anchor'] == mode) & (scan['target'] == target) & (scan['carbon'] == 'GPP')]
        if sub.empty:
            continue
        pv = sub.pivot_table(index='window', columns=['vi_index', 'metric'], values='beta_days_per_sd')
        print(f"\n[{mode} anchors] GPP -> {target}: within-site beta (days per +1 SD), rows = window")
        print(pv.round(2).to_string())

# ---- figures: grouped bars (cum vs mean) with 95% CI, one panel per target
for (mode, vi, cv), sub in scan.groupby(['anchor', 'vi_index', 'carbon']):
    targets = [t for t in ec.TARGETS if t in sub['target'].unique()]
    fig, axes = plt.subplots(1, len(targets), figsize=(5.2 * len(targets), 4.2), sharey=True, squeeze=False)
    for ax, target in zip(axes[0], targets):
        s = sub[sub['target'] == target]
        wins = list(dict.fromkeys(s['window']))
        x = np.arange(len(wins))
        for k, (metric, color) in enumerate((('cum', '#2b6cb0'), ('mean', '#dd6b20'))):
            m = s[s['metric'] == metric].set_index('window').reindex(wins)
            ax.bar(x + (k - 0.5) * 0.38, m['beta_days_per_sd'], 0.38, color=color,
                   yerr=1.96 * m['std_err'], capsize=2, label=f'{cv} {metric}')
            for xi, (b, p) in enumerate(zip(m['beta_days_per_sd'], m['p_value'])):
                if np.isfinite(b) and ec.stars(p):
                    ax.text(xi + (k - 0.5) * 0.38, b, ec.stars(p), ha='center',
                            va='bottom' if b >= 0 else 'top', fontsize=8)
        ax.axhline(0, color='k', lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([w_.replace('_to_', '→') for w_ in wins], rotation=40, ha='right', fontsize=8)
        ax.set_title(f'{target}')
    axes[0][0].set_ylabel('EOS shift (days per +1 within-site SD)')
    axes[0][0].legend(fontsize=8)
    fig.suptitle(f'{cv} vs EOS - {vi} - {mode} anchors', y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f'{mode}_{vi}_{cv}.png', dpi=140, bbox_inches='tight')
    plt.close(fig)
print(f"Figures -> '{FIG_DIR}'")
