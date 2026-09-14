import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "climate_clusters"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

CLUSTERS_CSV = DATA_DIR / "site_year_climate_clusters.csv"           # script 10 output
CLUSTER_MEANS_CSV = DATA_DIR / "climate_cluster_mean_timeseries.csv"  # script 10 output
PROJECTION_CSV = DATA_DIR / "site_year_climate_cluster_projection_2d.csv"  # script 10 output

CLIMATE_VARS = ['TA_F', 'SW_IN_F', 'VPD_F', 'P_F']
VAR_LABELS = {
    'TA_F': 'Air temperature (\N{DEGREE SIGN}C)',
    'SW_IN_F': 'Incoming shortwave radiation (W m$^{-2}$)',
    'VPD_F': 'VPD (hPa)',
    'P_F': 'Precipitation (mm)',
}

if not os.path.exists(CLUSTERS_CSV):
    raise FileNotFoundError(f"Missing '{CLUSTERS_CSV}'. Run the climate-clustering script first.")
if not os.path.exists(CLUSTER_MEANS_CSV):
    raise FileNotFoundError(f"Missing '{CLUSTER_MEANS_CSV}'. Run the climate-clustering script first.")

clusters_df = pd.read_csv(CLUSTERS_CSV)
means_df = pd.read_csv(CLUSTER_MEANS_CSV)

n_clusters = clusters_df['climate_cluster'].nunique()
colors = plt.cm.tab10(np.linspace(0, 1, max(n_clusters, 1)))

print(f"Plotting {n_clusters} climate clusters to '{FIGURE_DIR}'...")

# ---------------------------------------------------------------------------
# 1. Cluster sizes
# ---------------------------------------------------------------------------
cluster_sizes = clusters_df['climate_cluster'].value_counts().sort_index()

fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(cluster_sizes.index.astype(str), cluster_sizes.values,
       color=[colors[i] for i in cluster_sizes.index])
ax.set_xlabel('Climate cluster')
ax.set_ylabel('Number of site-years')
ax.set_title(f'Site-years per climate cluster (k={n_clusters})')
ax.grid(True, axis='y', linestyle='--', alpha=0.4)
fig.tight_layout()

out_path = FIGURE_DIR / 'climate_cluster_sizes.png'
fig.savefig(out_path, dpi=150)
plt.close(fig)
print(f"Saved '{out_path}'.")

# ---------------------------------------------------------------------------
# 2. Mean seasonal curve per variable, one line per cluster - lets the 8
#    groups actually be interpreted (e.g. "cluster 3 = warm, high-
#    radiation, low-precip") rather than just being arbitrary labels.
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
axes = axes.flatten()

for ax, var in zip(axes, CLIMATE_VARS):
    var_df = means_df[means_df['variable'] == var]
    if var_df.empty:
        ax.text(0.5, 0.5, f"'{var}' not found\nin {CLUSTER_MEANS_CSV.name}",
                ha='center', va='center', fontsize=9, color='gray', transform=ax.transAxes)
        ax.set_title(var)
        continue

    for cluster_id in sorted(var_df['climate_cluster'].unique()):
        curve = var_df[var_df['climate_cluster'] == cluster_id].sort_values('doy')
        n_members = int((clusters_df['climate_cluster'] == cluster_id).sum())
        ax.plot(curve['doy'], curve['mean_value'], color=colors[int(cluster_id)],
                 label=f'Cluster {cluster_id} (n={n_members})', linewidth=2)

    ax.set_title(var)
    ax.set_ylabel(VAR_LABELS.get(var, var))
    ax.grid(True, linestyle='--', alpha=0.4)

axes[-2].set_xlabel('Day of year')
axes[-1].set_xlabel('Day of year')
# Single shared legend below the whole figure rather than repeating it in
# every subplot - all 4 panels use the same cluster -> color mapping.
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=8, bbox_to_anchor=(0.5, -0.03))
fig.suptitle('Mean growing-season climate curve by cluster', fontsize=14, fontweight='bold')
fig.tight_layout(rect=[0, 0.05, 1, 0.96])

out_path = FIGURE_DIR / 'climate_cluster_profiles.png'
fig.savefig(out_path, dpi=150)
plt.close(fig)
print(f"Saved '{out_path}'.")

# ---------------------------------------------------------------------------
# 3. Separation plot - 2D projection of the clustering feature space,
#    colored by cluster, so you can see how distinct the 8 groups actually
#    are (tight, well-separated blobs vs. overlapping clusters). Drought
#    site-years (from the clustering script's within-site precip/VPD
#    anomaly flag) are additionally ringed in black, since they're an
#    explicit feature in the clustering now, not just incidentally captured.
# ---------------------------------------------------------------------------
if not os.path.exists(PROJECTION_CSV):
    print(f"\nNote: '{PROJECTION_CSV}' not found - skipping the separation plot. "
          "Re-run the climate-clustering script to generate it.")
else:
    proj_df = pd.read_csv(PROJECTION_CSV)
    if 'is_drought' in clusters_df.columns:
        proj_df = proj_df.merge(clusters_df[['site_id', 'year', 'is_drought']], on=['site_id', 'year'], how='left')
    else:
        proj_df['is_drought'] = False

    fig, ax = plt.subplots(figsize=(7, 6))
    for cluster_id in sorted(proj_df['climate_cluster'].unique()):
        member = proj_df[proj_df['climate_cluster'] == cluster_id]
        ax.scatter(member['pc1'], member['pc2'], color=colors[int(cluster_id)],
                   alpha=0.7, s=35, edgecolor='white', linewidth=0.4,
                   label=f'Cluster {cluster_id} (n={len(member)})')

    drought_pts = proj_df[proj_df['is_drought'] == True]
    if not drought_pts.empty:
        ax.scatter(drought_pts['pc1'], drought_pts['pc2'], facecolors='none',
                   edgecolors='black', linewidth=1.2, s=90,
                   label=f'Drought site-year (n={len(drought_pts)})')

    ax.set_xlabel('PC1 (climate feature space)')
    ax.set_ylabel('PC2 (climate feature space)')
    ax.set_title('Climate cluster separation (2D projection)')
    ax.legend(loc='best', fontsize=8, ncol=2)
    ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()

    out_path = FIGURE_DIR / 'climate_cluster_separation.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved '{out_path}'.")