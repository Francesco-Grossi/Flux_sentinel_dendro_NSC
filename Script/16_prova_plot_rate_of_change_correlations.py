import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "rate_of_change_correlations"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

SITEYEAR_CSV = DATA_DIR / "rate_of_change_predictors_by_site_year_index.csv"  # script 15 output
CORR_CSV = DATA_DIR / "rate_of_change_correlations.csv"                      # script 15 output

VI_INDICES = ['NDVI', 'NIRv']
VI_COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#38a169'}

if not os.path.exists(SITEYEAR_CSV):
    raise FileNotFoundError(f"Missing '{SITEYEAR_CSV}'. Run the rate-of-change correlation script first.")
if not os.path.exists(CORR_CSV):
    raise FileNotFoundError(f"Missing '{CORR_CSV}'. Run the rate-of-change correlation script first.")

analysis_df = pd.read_csv(SITEYEAR_CSV)
corr_df = pd.read_csv(CORR_CSV)


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


def plot_panel(ax, vi, autumn_p, pred):
    sub = analysis_df[analysis_df['vi_index'] == vi]
    missing = [c for c in (pred, autumn_p) if c not in sub.columns]
    if missing:
        ax.text(0.5, 0.5, f"column(s) {missing} not found\nre-run script 15",
                ha='center', va='center', fontsize=9, color='gray', transform=ax.transAxes)
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
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, color='#c53030', linewidth=2,
             label=f"y = {slope:.3g}x + {intercept:.3g}")

    ax.set_xlabel(pred)
    ax.set_ylabel(autumn_p)
    ax.set_title(f"{vi}\nr = {r_value:.2f}, p = {p_value:.3g}, n = {len(pair)}")
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.4)


# ---------------------------------------------------------------------------
# 1. One figure per (autumn_parameter, predictor) pair, NDVI and NIRv side
#    by side - reuses exactly the pairs script 15 correlated.
# ---------------------------------------------------------------------------
pairs = corr_df[['autumn_parameter', 'predictor']].drop_duplicates()
print(f"Plotting {len(pairs)} rate-of-change correlation pairs "
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
