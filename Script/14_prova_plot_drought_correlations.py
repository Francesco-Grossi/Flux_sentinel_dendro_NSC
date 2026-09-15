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
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "drought_correlations"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

COMPARISON_CSV = DATA_DIR / "drought_vs_nondrought_comparison.csv"  # script 13 output
SITEYEAR_CSV = DATA_DIR / "phenology_flux_predictors_with_drought_flag.csv"  # script 13 output

VI_INDICES = ['NDVI', 'NIRv']
GROUP_COLORS = {'drought': '#c53030', 'non_drought': '#2b6cb0'}
GROUPS = {'Drought': True, 'Non-drought': False}

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

# ---------------------------------------------------------------------------
# 2. 2x2 scatter+fit grid per (autumn_parameter, predictor) pair - columns
#    are Drought / Non-drought, rows are NDVI / NIRv, so the same pair can
#    be compared across both dimensions in one figure.
# ---------------------------------------------------------------------------
if not os.path.exists(SITEYEAR_CSV):
    print(f"\nNote: '{SITEYEAR_CSV}' not found - skipping the 2x2 scatter grids. "
          "Re-run the drought-correlation script to generate it.")
else:
    siteyear_df = pd.read_csv(SITEYEAR_CSV)
    scatter_dir = FIGURE_DIR / "scatter_2x2"
    scatter_dir.mkdir(parents=True, exist_ok=True)

    pairs = comparison_df[['autumn_parameter', 'predictor']].drop_duplicates()
    print(f"\nPlotting {len(pairs)} 2x2 drought x VI-index scatter grids to '{scatter_dir}'...")

    n_plotted = 0
    for _, prow in pairs.iterrows():
        autumn_p, pred = prow['autumn_parameter'], prow['predictor']
        missing = [c for c in (pred, autumn_p) if c not in siteyear_df.columns]
        if missing:
            continue

        fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex='col', sharey='row')
        for row_idx, vi in enumerate(VI_INDICES):
            for col_idx, (group_label, flag_value) in enumerate(GROUPS.items()):
                ax = axes[row_idx, col_idx]
                sub = siteyear_df[(siteyear_df['vi_index'] == vi) & (siteyear_df['is_drought'] == flag_value)]
                pair = sub[[pred, autumn_p]].dropna()

                if len(pair) < 2:
                    ax.text(0.5, 0.5, "not enough data", ha='center', va='center',
                            fontsize=9, color='gray', transform=ax.transAxes)
                else:
                    x = pair[pred].to_numpy(dtype=float)
                    y = pair[autumn_p].to_numpy(dtype=float)
                    color = GROUP_COLORS['drought' if flag_value else 'non_drought']
                    ax.scatter(x, y, alpha=0.6, s=30, color=color, edgecolor='white', linewidth=0.5)
                    if np.std(x) > 0 and np.std(y) > 0:
                        slope, intercept, r_value, p_value, std_err = linregress(x, y)
                        x_line = np.linspace(x.min(), x.max(), 100)
                        ax.plot(x_line, slope * x_line + intercept, color='#333333', linewidth=2)
                        ax.text(0.03, 0.96, f"r = {r_value:.2f}, p = {p_value:.3g}, n = {len(pair)}",
                                transform=ax.transAxes, ha='left', va='top', fontsize=8)

                if row_idx == 0:
                    ax.set_title(group_label, fontsize=11, fontweight='bold')
                if col_idx == 0:
                    ax.set_ylabel(f"{vi}\n{autumn_p}")
                if row_idx == 1:
                    ax.set_xlabel(pred)
                ax.grid(True, linestyle='--', alpha=0.4)

        fig.suptitle(f"{autumn_p} vs {pred}: drought x VI index", fontsize=13, fontweight='bold')
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        out_path = scatter_dir / f"{safe_filename(autumn_p)}_vs_{safe_filename(pred)}_2x2.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        n_plotted += 1

    print(f"Done. {n_plotted} 2x2 scatter grids saved under '{scatter_dir}'.")