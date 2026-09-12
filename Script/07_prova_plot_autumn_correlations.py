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

SITEYEAR_CSV = DATA_DIR / "phenology_gpp_predictors_by_site_year_index.csv"  # script 6 output
CORR_CSV = DATA_DIR / "autumn_phenology_correlations.csv"                    # script 6 output

if not os.path.exists(SITEYEAR_CSV):
    raise FileNotFoundError(f"Missing '{SITEYEAR_CSV}'. Run the autumn-phenology-correlation script first.")
if not os.path.exists(CORR_CSV):
    raise FileNotFoundError(f"Missing '{CORR_CSV}'. Run the autumn-phenology-correlation script first.")

analysis_df = pd.read_csv(SITEYEAR_CSV)
corr_df = pd.read_csv(CORR_CSV)

print(f"Plotting {len(corr_df)} autumn-parameter x predictor pairs "
      f"(one scatter + best-fit-line plot each) to '{FIGURE_DIR}'...")


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


# ---------------------------------------------------------------------------
# 1. One scatter + best-linear-fit plot per (vi_index, autumn_parameter,
#    predictor) row already computed in CORR_CSV - reuses exactly the same
#    pairs script 6 correlated (same MIN_PAIRS_FOR_CORR / non-zero-std
#    gating), so a plot exists if and only if a correlation was reported.
# ---------------------------------------------------------------------------
n_plotted = 0
for _, row in corr_df.iterrows():
    vi, autumn_p, pred = row['vi_index'], row['autumn_parameter'], row['predictor']

    sub = analysis_df[analysis_df['vi_index'] == vi]
    pair = sub[[pred, autumn_p]].dropna()
    if len(pair) < 2:
        continue

    x = pair[pred].to_numpy(dtype=float)
    y = pair[autumn_p].to_numpy(dtype=float)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x, y, alpha=0.6, s=30, color='#2b6cb0', edgecolor='white', linewidth=0.5)

    x_line = np.linspace(x.min(), x.max(), 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, color='#c53030', linewidth=2,
             label=f"y = {slope:.3g}x + {intercept:.3g}")

    ax.set_xlabel(pred)
    ax.set_ylabel(autumn_p)
    ax.set_title(f"{vi}: {autumn_p} vs {pred}\n"
                 f"r = {r_value:.2f}, p = {p_value:.3g}, n = {len(pair)}")
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.4)

    fig.tight_layout()

    out_dir = FIGURE_DIR / vi
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{safe_filename(autumn_p)}_vs_{safe_filename(pred)}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    n_plotted += 1

print(f"\nDone. {n_plotted} plots saved under '{FIGURE_DIR}' (one subfolder per VI index).")
