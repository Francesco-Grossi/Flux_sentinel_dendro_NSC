import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "drought_correlations"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

COMPARISON_CSV = DATA_DIR / "drought_vs_nondrought_comparison.csv"  # script 13 output

VI_INDICES = ['NDVI', 'NIRv']
GROUP_COLORS = {'drought': '#c53030', 'non_drought': '#2b6cb0'}

if not os.path.exists(COMPARISON_CSV):
    raise FileNotFoundError(f"Missing '{COMPARISON_CSV}'. Run the drought-correlation script first.")

comparison_df = pd.read_csv(COMPARISON_CSV)


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


def plot_panel(ax, vi, autumn_p):
    sub = comparison_df[(comparison_df['vi_index'] == vi) & (comparison_df['autumn_parameter'] == autumn_p)]
    if sub.empty:
        ax.text(0.5, 0.5, "no data", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return

    sub = sub.sort_values('predictor')
    predictors = sub['predictor'].tolist()
    x = np.arange(len(predictors))
    width = 0.38

    r_drought = sub['pearson_r_drought'].to_numpy(dtype=float) if 'pearson_r_drought' in sub.columns else np.full(len(sub), np.nan)
    r_non_drought = sub['pearson_r_non_drought'].to_numpy(dtype=float) if 'pearson_r_non_drought' in sub.columns else np.full(len(sub), np.nan)
    n_drought = sub['n_drought'].to_numpy() if 'n_drought' in sub.columns else None
    n_non_drought = sub['n_non_drought'].to_numpy() if 'n_non_drought' in sub.columns else None

    ax.bar(x - width / 2, r_drought, width,
           label='Drought' + (f' (n~{int(np.nanmedian(n_drought))})' if n_drought is not None else ''),
           color=GROUP_COLORS['drought'])
    ax.bar(x + width / 2, r_non_drought, width,
           label='Non-drought' + (f' (n~{int(np.nanmedian(n_non_drought))})' if n_non_drought is not None else ''),
           color=GROUP_COLORS['non_drought'])

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(predictors, rotation=60, ha='right', fontsize=8)
    ax.set_ylabel('Pearson r')
    ax.set_ylim(-1, 1)
    ax.set_title(vi)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, axis='y', linestyle='--', alpha=0.4)


# ---------------------------------------------------------------------------
# 1. One figure per autumn parameter, NDVI and NIRv side by side - each
#    panel shows every predictor's drought-vs-non-drought r as paired bars,
#    so a flipped sign or a big gap between the two bars is immediately
#    visible.
# ---------------------------------------------------------------------------
autumn_params = sorted(comparison_df['autumn_parameter'].unique())
print(f"Plotting {len(autumn_params)} autumn-parameter figures (NDVI and NIRv side by side) to '{FIGURE_DIR}'...")

for autumn_p in autumn_params:
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, vi in zip(axes, VI_INDICES):
        plot_panel(ax, vi, autumn_p)

    fig.suptitle(f"{autumn_p}: drought vs non-drought correlation by predictor", fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_path = FIGURE_DIR / f"{safe_filename(autumn_p)}_drought_comparison.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved '{out_path}'.")

print("\nDone.")
