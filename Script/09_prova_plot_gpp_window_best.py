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
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "gpp_window_best_correlations"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # script 5 output
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                            # script 4 output
BEST_CSV = DATA_DIR / "gpp_window_best_correlations.csv"                     # script 8 output

VI_INDICES = ['NDVI', 'NIRv']
VI_COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#38a169'}

GPP_COLUMN = 'GPP_NT_VUT_REF'
GPP_PLAUSIBLE_RANGE = (-5, 50)
MIN_COVERAGE_FRAC = 0.7   # must match script 8's setting to reproduce the same window sums
MIN_FIT_CORR = 0.8        # must match script 8's setting to reproduce the same site-year set

for path in (PHENOLOGY_CSV, FLUX_CSV, BEST_CSV):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing '{path}'. Run the phenology, merge, and GPP-window-scan "
                                 "scripts first.")


# ---------------------------------------------------------------------------
# 1. Rebuild the same (site_id, year) x day-of-year GPP prefix-sum matrix
#    used in the window-scan script, so the exact window sums behind each
#    "best window" row can be recomputed for plotting.
# ---------------------------------------------------------------------------
pheno = pd.read_csv(PHENOLOGY_CSV)
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')
              & (pheno['corr'] >= MIN_FIT_CORR)].copy()
pheno['year'] = pheno['year'].astype(int)

flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year
flux['doy'] = flux['date'].dt.dayofyear

gpp_lo, gpp_hi = GPP_PLAUSIBLE_RANGE
implausible = (flux[GPP_COLUMN] < gpp_lo) | (flux[GPP_COLUMN] > gpp_hi)
flux.loc[implausible, GPP_COLUMN] = np.nan

flux_pivot = flux.pivot_table(index=['site_id', 'year'], columns='doy', values=GPP_COLUMN, aggfunc='mean')
flux_pivot = flux_pivot.reindex(columns=range(1, 367))

matrix_keys = list(flux_pivot.index)
matrix_pos = {key: i for i, key in enumerate(matrix_keys)}
matrix_values = flux_pivot.to_numpy(dtype=float)

valid_mask = ~np.isnan(matrix_values)
filled = np.nan_to_num(matrix_values, nan=0.0)
cumsum_vals = np.concatenate([np.zeros((filled.shape[0], 1)), np.cumsum(filled, axis=1)], axis=1)
cumsum_counts = np.concatenate([np.zeros((valid_mask.shape[0], 1)),
                                 np.cumsum(valid_mask.astype(float), axis=1)], axis=1)


def window_sums_for_positions(positions, start_doy, end_doy):
    s = cumsum_vals[positions, end_doy] - cumsum_vals[positions, start_doy - 1]
    counts = cumsum_counts[positions, end_doy] - cumsum_counts[positions, start_doy - 1]
    coverage = counts / (end_doy - start_doy + 1)
    return np.where(coverage >= MIN_COVERAGE_FRAC, s, np.nan)


def get_window_pair(vi, autumn_p, start_doy, end_doy):
    """Recompute the (window-sum, autumn-parameter) pairs for one VI index
    and one specific window, for plotting."""
    sub = pheno[pheno['vi_index'] == vi].copy()
    if autumn_p not in sub.columns:
        return None, None
    sub['_pos'] = sub.apply(lambda r: matrix_pos.get((r['site_id'], r['year'])), axis=1)
    sub = sub.dropna(subset=['_pos', autumn_p]).copy()
    if sub.empty:
        return None, None
    sub['_pos'] = sub['_pos'].astype(int)
    positions = sub['_pos'].to_numpy()

    w_sum = window_sums_for_positions(positions, start_doy, end_doy)
    y = sub[autumn_p].to_numpy(dtype=float)
    mask = ~np.isnan(w_sum) & ~np.isnan(y)
    if mask.sum() < 2:
        return None, None
    return w_sum[mask], y[mask]


def safe_filename(text):
    return text.replace('/', '_').replace(' ', '_')


def plot_panel(ax, vi, autumn_p, best_row):
    if best_row is None:
        ax.text(0.5, 0.5, f"no window found\nfor {vi}", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return

    start_doy, end_doy = int(best_row['window_start_doy']), int(best_row['window_end_doy'])
    x, y = get_window_pair(vi, autumn_p, start_doy, end_doy)
    if x is None:
        ax.text(0.5, 0.5, "not enough data", ha='center', va='center',
                fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(vi)
        return

    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    color = VI_COLORS.get(vi, '#2b6cb0')
    ax.scatter(x, y, alpha=0.6, s=30, color=color, edgecolor='white', linewidth=0.5)

    x_line = np.linspace(x.min(), x.max(), 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, color='#c53030', linewidth=2,
             label=f"y = {slope:.3g}x + {intercept:.3g}")

    ax.set_xlabel(f"cumulative GPP, DOY {start_doy}-{end_doy}")
    ax.set_ylabel(autumn_p)
    ax.set_title(f"{vi}: DOY {start_doy}-{end_doy}\nr = {r_value:.2f}, p = {p_value:.3g}, n = {len(x)}")
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.4)


# ---------------------------------------------------------------------------
# 2. One figure per (autumn_parameter, correlation_type), NDVI and NIRv
#    side by side - each panel uses ITS OWN best window (start/end may
#    differ between the two indices for the same autumn parameter).
# ---------------------------------------------------------------------------
best_df = pd.read_csv(BEST_CSV)

combos = best_df[['autumn_parameter', 'correlation_type']].drop_duplicates()
print(f"Plotting {len(combos)} (autumn-parameter x correlation-type) figures "
      f"(NDVI and NIRv side by side, each using its own best window) to '{FIGURE_DIR}'...")

n_plotted = 0
for _, crow in combos.iterrows():
    autumn_p, corr_type = crow['autumn_parameter'], crow['correlation_type']

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, vi in zip(axes, VI_INDICES):
        match = best_df[(best_df['vi_index'] == vi) &
                         (best_df['autumn_parameter'] == autumn_p) &
                         (best_df['correlation_type'] == corr_type)]
        best_row = match.iloc[0] if not match.empty else None
        plot_panel(ax, vi, autumn_p, best_row)

    label = "Strongest positive correlation" if corr_type == 'most_positive' else "Strongest negative correlation"
    fig.suptitle(f"{autumn_p} - {label}", fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_path = FIGURE_DIR / f"{safe_filename(autumn_p)}_{corr_type}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    n_plotted += 1

print(f"\nDone. {n_plotted} plots saved under '{FIGURE_DIR}'.")