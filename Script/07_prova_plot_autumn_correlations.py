import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/, data/ and figure/ are sibling folders under the repo root, so
# this works regardless of the working directory the script is launched
# from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "autumn_phenology_correlations"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # script 6 output
CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"                    # script 6 output

# growing_season_length is excluded from these plots per request - it's
# EOS50 minus leaf_out_50, i.e. partly a function of the response variable
# itself when the autumn parameter is EOS50, which makes it a less useful
# thing to visually compare side by side with the other predictors here.
EXCLUDED_PREDICTORS = {'growing_season_length'}

VI_INDICES = ['NDVI', 'NIRv']
VI_COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#38a169'}

if not os.path.exists(SITEYEAR_CSV):
    raise FileNotFoundError(f"Missing '{SITEYEAR_CSV}'. Run the autumn-phenology-correlation script first.")
if not os.path.exists(CORR_CSV):
    raise FileNotFoundError(f"Missing '{CORR_CSV}'. Run the autumn-phenology-correlation script first.")

analysis_df = pd.read_csv(SITEYEAR_CSV)
corr_df = pd.read_csv(CORR_CSV)
corr_df = corr_df[~corr_df['predictor'].isin(EXCLUDED_PREDICTORS)].copy()


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


def plot_panel(ax, vi, autumn_p, pred):
    """Draw one index's scatter + best-fit line into the given axis.
    Returns True if something was plotted, False if the data/column wasn't
    available (in which case the axis is left with a 'no data' message)."""
    sub = analysis_df[analysis_df['vi_index'] == vi]
    missing = [c for c in (pred, autumn_p) if c not in sub.columns]
    if missing:
        ax.text(0.5, 0.5, f"column(s) {missing} not found\nre-run script 6",
                ha='center', va='center', fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return False

    pair = sub[[pred, autumn_p]].dropna()
    if len(pair) < 2:
        ax.text(0.5, 0.5, "not enough data", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return False

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
    return True


# ---------------------------------------------------------------------------
# 1. One figure per (autumn_parameter, predictor) pair, with NDVI and NIRv
#    side by side as two subplots - lets the two indices be compared
#    directly on the same axes scale/layout rather than as separate files.
# ---------------------------------------------------------------------------
pairs = corr_df[['autumn_parameter', 'predictor']].drop_duplicates()
print(f"Plotting {len(pairs)} autumn-parameter x predictor pairs "
      f"(NDVI and NIRv side by side, growing_season_length excluded) to '{FIGURE_DIR}'...")

n_plotted = 0
for _, prow in pairs.iterrows():
    autumn_p, pred = prow['autumn_parameter'], prow['predictor']

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    any_plotted = False
    for ax, vi in zip(axes, VI_INDICES):
        any_plotted |= plot_panel(ax, vi, autumn_p, pred)

    fig.suptitle(f"{autumn_p} vs {pred}", fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_path = FIGURE_DIR / f"{safe_filename(autumn_p)}_vs_{safe_filename(pred)}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    if any_plotted:
        n_plotted += 1

print(f"\nDone. {n_plotted} side-by-side (NDVI | NIRv) plots saved under '{FIGURE_DIR}'.")