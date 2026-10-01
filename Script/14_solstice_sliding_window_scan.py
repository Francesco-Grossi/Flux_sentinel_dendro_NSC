"""
PIPELINE STEP 14 - Sliding-window scan of carbon flux around the summer
solstice (Notion D.1 test (b), after Zohner et al. 2023 Fig. S12).

Question: which period of the season has the MOST NEGATIVE relationship
between carbon uptake and EOS, and is it close to the summer solstice?

For each window length L (15 / 30 d) and start offset k = -75 ... +75 d
relative to the solstice, take the mean daily flux over
[SOL + k, SOL + k + L - 1] and regress EOS10 / EOS50 / EOS90 on it
(site-random-intercept LME, predictor z-scored; beta = days per +1 SD).
Windows do not depend on any EOS value, so there is no circularity.

Input : data/phenology_double_logistic_by_site_year_index.csv (step 5)
        data/fluxnet_landsat_merged.csv                        (step 4)
        data/npp_luo2025_daily.csv                             (optional)
Output: data/eos_solstice_sliding_scan.csv
        data/eos_solstice_sliding_scan_minima.csv   (most-negative window per curve)
        figure/eos_solstice_sliding_scan/<carbon>_L<len>.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

OFFSETS = list(range(-75, 76, 5))
LENGTHS = (15, 30)
OUTPUT_CSV = ec.DATA_DIR / "eos_solstice_sliding_scan.csv"
OUTPUT_MIN_CSV = ec.DATA_DIR / "eos_solstice_sliding_scan_minima.csv"
FIG_DIR = ec.FIGURE_DIR / "eos_solstice_sliding_scan"
FIG_DIR.mkdir(parents=True, exist_ok=True)

pheno = ec.load_phenology()
flux = ec.load_flux_daily()
carbon_vars = [v for v in ('GPP', 'NPP', 'NPPd', 'CUEd') if v in flux.columns]
lk = ec.FluxLookup(flux, carbon_vars)

sy = pheno[['site_id', 'year']].drop_duplicates().reset_index(drop=True)
sol = sy['year'].map(ec.solstice_doy).astype(float)

rows = []
for cv in carbon_vars:
    for L in LENGTHS:
        for k in OFFSETS:
            x = [lk.window_mean(s, y, cv, so + k, so + k + L - 1) for s, y, so in zip(sy['site_id'], sy['year'], sol)]
            xs = sy.assign(x=x)
            for vi in sorted(pheno['vi_index'].unique()):
                d = pheno[pheno['vi_index'] == vi].merge(xs, on=['site_id', 'year'], how='left')
                for target in ec.TARGETS:
                    res = ec.lme_slope(d, target, 'x')
                    if res:
                        rows.append({'carbon': cv, 'vi_index': vi, 'target': target, 'length_days': L,
                                     'offset_from_solstice_days': k, **res})

scan = pd.DataFrame(rows)
scan.to_csv(OUTPUT_CSV, index=False)
print(f"{len(scan)} cells -> '{OUTPUT_CSV}'")
if scan.empty:
    raise SystemExit("No window had enough data.")

minima = (scan.sort_values('beta_days_per_sd')
          .groupby(['carbon', 'vi_index', 'target', 'length_days'], as_index=False).first()
          [['carbon', 'vi_index', 'target', 'length_days', 'offset_from_solstice_days',
            'beta_days_per_sd', 'p_value', 'n_obs', 'n_sites']])
minima['near_solstice_pm15d'] = minima['offset_from_solstice_days'].abs() <= 15
minima.to_csv(OUTPUT_MIN_CSV, index=False)
print("\nMost negative window per curve (offset in days from solstice; window START):")
print(minima.round(3).to_string(index=False))

for (cv, L), sub in scan.groupby(['carbon', 'length_days']):
    vis = sorted(sub['vi_index'].unique())
    fig, axes = plt.subplots(len(vis), 1, figsize=(7.5, 3.2 * len(vis)), sharex=True, squeeze=False)
    for ax, vi in zip(axes[:, 0], vis):
        for target, color in zip(ec.TARGETS, ('#c53030', '#dd6b20', '#2b6cb0')):
            s = sub[(sub['vi_index'] == vi) & (sub['target'] == target)].sort_values('offset_from_solstice_days')
            if s.empty:
                continue
            ax.plot(s['offset_from_solstice_days'], s['beta_days_per_sd'], color=color, label=target, lw=1.6)
            ax.fill_between(s['offset_from_solstice_days'],
                            s['beta_days_per_sd'] - 1.96 * s['std_err'],
                            s['beta_days_per_sd'] + 1.96 * s['std_err'], color=color, alpha=0.12)
        ax.axhline(0, color='k', lw=0.6)
        ax.axvline(0, color='gray', ls='--', lw=0.8)
        ax.set_ylabel(f'{vi}\nEOS shift (d / SD)')
        ax.legend(fontsize=8, ncol=3)
    axes[-1, 0].set_xlabel(f'window start relative to summer solstice (days); window = {L} d')
    fig.suptitle(f'{cv} mean flux vs EOS: sliding window', y=1.0)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f'{cv}_L{L}.png', dpi=140, bbox_inches='tight')
    plt.close(fig)
print(f"Figures -> '{FIG_DIR}'")
