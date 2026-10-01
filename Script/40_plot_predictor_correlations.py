"""
PIPELINE STEP 10 - Scatter + best-fit line for every (autumn_parameter x
predictor) pair from step 6's correlation table, NDVI and NIRv side by
side, so the correlation strength can actually be looked at rather than
just read off a table of r/p values.

Inputs: data/phenology_flux_predictors_by_site_year_index.csv  (step 6)
        data/autumn_phenology_correlations.csv                 (step 6)

Output: figure/predictor_correlations/<autumn_param>_vs_<predictor>.png
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

SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # step 6
CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"                    # step 6

VI_INDICES = ['NDVI', 'NIRv']
VI_COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#38a169'}

for p in (SITEYEAR_CSV, CORR_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run 06_autumn_phenology_predictors.py first.")

analysis_df = pd.read_csv(SITEYEAR_CSV)
corr_df = pd.read_csv(CORR_CSV)


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


def plot_panel(ax, vi, autumn_p, pred):
    sub = analysis_df[analysis_df['vi_index'] == vi]
    if pred not in sub.columns or autumn_p not in sub.columns:
        ax.text(0.5, 0.5, "column not found\nre-run step 06", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return

    pair = sub[[pred, autumn_p]].dropna()
    if len(pair) < 2:
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

    ax.set_xlabel(pred)
    ax.set_ylabel(autumn_p)
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

print(f"\nDone. {n_plotted} plots saved under '{FIGURE_DIR}'.")
