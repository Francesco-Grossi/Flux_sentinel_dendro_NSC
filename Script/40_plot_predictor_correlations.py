"""
PIPELINE STEP 40 - Scatter + best-fit line for every (autumn_parameter x
predictor) pair from step 28's correlation table, NDVI and NIRv side by
side, so the correlation strength can actually be looked at rather than
just read off a table of r/p values.

For the split-GPP effect, compare the pairs of figures
    <EOS>_vs_gpp_solstice_to_eos90.png   and   <EOS>_vs_gpp_solstice_to_eos90_fixed.png
    <EOS>_vs_gpp_solstice_to_eos10.png   and   <EOS>_vs_gpp_solstice_to_eos10_fixed.png
(window ending at the same year's EOS vs at the site's mean EOS; see step 28).

Inputs: data/phenology_flux_predictors_by_site_year_index.csv  (step 28)
        data/autumn_phenology_correlations.csv                 (step 28)

Every pair is drawn twice:
    pooled        all site-years as they are. Mixes differences BETWEEN sites
                  (productive vs unproductive, dry-summer vs summer-green) with
                  year-to-year changes.
    within-site   each value minus its site's mean (sites with >= MIN_YEARS
                  years). Only year-to-year changes at a site remain - this is
                  the version that corresponds to the hypothesis tests (step 41).

Output: figure/predictor_correlations/<autumn_param>_vs_<predictor>.png
        figure/predictor_correlations_within_site/<autumn_param>_vs_<predictor>.png
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "predictor_correlations"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
WITHIN_DIR = FIGURE_DIR.parent / "predictor_correlations_within_site"
WITHIN_DIR.mkdir(parents=True, exist_ok=True)
MIN_YEARS = 3

SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # step 28
CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"                    # step 28

VI_INDICES = ['NDVI', 'NIRv']
VI_COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#38a169'}

for p in (SITEYEAR_CSV, CORR_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run 28_autumn_phenology_predictors.py first.")

analysis_df = pd.read_csv(SITEYEAR_CSV)
corr_df = pd.read_csv(CORR_CSV)


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


def plot_panel(ax, vi, autumn_p, pred, within=False):
    sub = analysis_df[analysis_df['vi_index'] == vi]
    if pred not in sub.columns or autumn_p not in sub.columns:
        ax.text(0.5, 0.5, "column not found\nre-run step 28", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return

    pair = sub[['site_id', pred, autumn_p]].dropna()
    if within:      # year minus site mean
        pair = pair[pair.groupby('site_id')['site_id'].transform('size') >= MIN_YEARS]
        pair = pair[[pred, autumn_p]] - pair.groupby('site_id')[[pred, autumn_p]].transform('mean')
    if len(pair) < 2 or pair[pred].std() == 0:
        ax.text(0.5, 0.5, "not enough data", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return

    x = pair[pred].to_numpy(dtype=float)
    y = pair[autumn_p].to_numpy(dtype=float)
    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    color = VI_COLORS.get(vi, '#2b6cb0')
    ax.scatter(x, y, alpha=0.6, s=30, color=color, edgecolor='white', linewidth=0.5)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, color='#c53030', linewidth=2,
             label=f"y = {slope:.3g}x + {intercept:.3g}")

    ax.set_xlabel(f"{pred} anomaly (year - site mean)" if within else pred)
    ax.set_ylabel(f"{autumn_p} anomaly (days)" if within else autumn_p)
    if within:
        ax.axhline(0, color='0.8', lw=0.8, zorder=0)
        ax.axvline(0, color='0.8', lw=0.8, zorder=0)
    ax.set_title(f"{vi}\nr = {r_value:.2f}, p = {p_value:.3g}, n = {len(pair)}")
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.4)


pairs = corr_df[['autumn_parameter', 'predictor']].drop_duplicates()
print(f"Plotting {len(pairs)} autumn-parameter x predictor pairs "
      f"(NDVI and NIRv side by side) to '{FIGURE_DIR}'...")

n_plotted = 0
for _, prow in pairs.iterrows():
    autumn_p, pred = prow['autumn_parameter'], prow['predictor']
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, vi in zip(axes, VI_INDICES):
        plot_panel(ax, vi, autumn_p, pred)
    fig.suptitle(f"{autumn_p} vs {pred}", fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out_path = FIGURE_DIR / f"{safe_filename(autumn_p)}_vs_{safe_filename(pred)}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    n_plotted += 1

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, vi in zip(axes, VI_INDICES):
        plot_panel(ax, vi, autumn_p, pred, within=True)
    fig.suptitle(f"{autumn_p} vs {pred} - within-site anomalies", fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(WITHIN_DIR / f"{safe_filename(autumn_p)}_vs_{safe_filename(pred)}.png", dpi=150)
    plt.close(fig)

print(f"\nDone. {n_plotted} pooled plots under '{FIGURE_DIR}', and the same {n_plotted} as within-site "
      f"anomalies under '{WITHIN_DIR}'.")
