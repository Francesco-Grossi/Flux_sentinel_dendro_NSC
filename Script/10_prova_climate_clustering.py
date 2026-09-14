import os
import sys
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
except ImportError:
    print("Installing scikit-learn (required for PCA/KMeans clustering)...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-learn"])
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

FLUX_CSV = DATA_DIR / "fluxnet_daily_all_vars_qc_filtered.csv"  # script 3 output

OUTPUT_CLUSTERS_CSV = DATA_DIR / "site_year_climate_clusters.csv"
OUTPUT_CLUSTER_MEANS_CSV = DATA_DIR / "climate_cluster_mean_timeseries.csv"
OUTPUT_CLUSTER_PROJECTION_CSV = DATA_DIR / "site_year_climate_cluster_projection_2d.csv"

# The 4 climate variables, used as their FULL daily time series (not summed
# or otherwise reduced to a single number) - this is what "not the sum, all
# 4 time series" means in practice: each variable contributes its whole
# seasonal shape to the clustering, not just its total/mean.
CLIMATE_VARS = ['TA_F', 'SW_IN_F', 'VPD_F', 'P_F']

# Growing season window (Mar 1 - Oct 31), matching the rest of this
# pipeline - a fixed 245-day template regardless of leap year, so every
# site-year contributes a same-length vector per variable.
WINDOW_START_DOY = 60
WINDOW_END_DOY = 304

# A site-year is dropped from clustering if any of the 4 variables has less
# than this fraction of its growing-season days actually present - small
# gaps are linearly interpolated, but a site-year that's mostly missing one
# variable shouldn't be clustered on a mostly-invented curve.
MIN_COMPLETENESS = 0.7

# Cumulative variance to retain per variable's PCA (captures the SHAPE of
# that variable's seasonal curve in a handful of components instead of
# clustering on ~245 raw, highly autocorrelated daily values directly).
PCA_VARIANCE_RETAINED = 0.90

N_CLUSTERS = 4
RANDOM_STATE = 42

if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run the growing-season QC filter script first.")


# ---------------------------------------------------------------------------
# 1. Load daily FLUXNET data, restrict to the growing-season window
# ---------------------------------------------------------------------------
flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year
flux['doy'] = flux['date'].dt.dayofyear
flux = flux[(flux['doy'] >= WINDOW_START_DOY) & (flux['doy'] <= WINDOW_END_DOY)].copy()

missing_vars = [v for v in CLIMATE_VARS if v not in flux.columns]
if missing_vars:
    raise ValueError(f"Missing climate variable column(s) in FLUX_CSV: {missing_vars}")

site_meta = flux.groupby('site_id')[['lat', 'lon', 'igbp']].first()

# ---------------------------------------------------------------------------
# 2. Build one (site_id, year) x day-of-year matrix PER VARIABLE, fill small
#    gaps by linear interpolation, and drop site-years with too much
#    missing data in any variable.
# ---------------------------------------------------------------------------
day_range = list(range(WINDOW_START_DOY, WINDOW_END_DOY + 1))
n_days = len(day_range)

var_matrices = {}
for var in CLIMATE_VARS:
    pivot = flux.pivot_table(index=['site_id', 'year'], columns='doy', values=var, aggfunc='mean')
    pivot = pivot.reindex(columns=day_range)
    var_matrices[var] = pivot

# Site-years present in every variable's pivot (should normally be all of
# them, since they all come from the same table, but this stays robust to
# a variable that happens to be entirely absent for some site-year).
common_keys = set(var_matrices[CLIMATE_VARS[0]].index)
for var in CLIMATE_VARS[1:]:
    common_keys &= set(var_matrices[var].index)
common_keys = sorted(common_keys)

kept_keys = []
for key in common_keys:
    ok = True
    for var in CLIMATE_VARS:
        row = var_matrices[var].loc[key]
        completeness = row.notna().mean()
        if completeness < MIN_COMPLETENESS:
            ok = False
    if ok:
        kept_keys.append(key)

n_dropped = len(common_keys) - len(kept_keys)
print(f"Site-years with >= {int(MIN_COMPLETENESS * 100)}% coverage in all 4 climate variables: "
      f"{len(kept_keys)} (dropped {n_dropped} for insufficient coverage in at least one variable).")

# Interpolate remaining small gaps (linear, within the growing-season
# window only - no extrapolation past the window edges) for the kept
# site-years, then stack into a (n_site_years, n_days) array per variable.
clean_matrices = {}
for var in CLIMATE_VARS:
    sub = var_matrices[var].loc[kept_keys]
    sub = sub.interpolate(axis=1, limit_direction='both')
    clean_matrices[var] = sub.to_numpy(dtype=float)

# ---------------------------------------------------------------------------
# 3. Per-variable PCA on the standardized daily curves - this is what
#    encodes "the whole time series shape" (timing, amplitude, curve
#    features) into a manageable number of features per variable, instead
#    of reducing each variable to a single sum/mean.
# ---------------------------------------------------------------------------
pc_blocks = []
pc_variable_labels = []
for var in CLIMATE_VARS:
    X = clean_matrices[var]
    X_scaled = StandardScaler().fit_transform(X)  # standardize each day-of-year column across site-years
    pca = PCA(n_components=PCA_VARIANCE_RETAINED, random_state=RANDOM_STATE)
    scores = pca.fit_transform(X_scaled)
    pc_blocks.append(scores)
    pc_variable_labels += [f"{var}_PC{i+1}" for i in range(scores.shape[1])]
    print(f"{var}: {scores.shape[1]} PCs retained to explain "
          f"{pca.explained_variance_ratio_.sum() * 100:.1f}% of variance.")

combined_features = np.concatenate(pc_blocks, axis=1)
# Re-standardize the combined feature set - PCA components from a variable
# with more retained PCs (or larger raw variance) would otherwise dominate
# the distance metric KMeans uses, even though each variable is meant to
# contribute equally to what "similar climate" means here.
combined_features_scaled = StandardScaler().fit_transform(combined_features)

# ---------------------------------------------------------------------------
# 4. K-means into 8 climate groups
# ---------------------------------------------------------------------------
kmeans = KMeans(n_clusters=N_CLUSTERS, random_state=RANDOM_STATE, n_init=10)
cluster_labels = kmeans.fit_predict(combined_features_scaled)

clusters_df = pd.DataFrame(kept_keys, columns=['site_id', 'year'])
clusters_df['climate_cluster'] = cluster_labels
clusters_df = clusters_df.merge(site_meta, on='site_id', how='left')
clusters_df.to_csv(OUTPUT_CLUSTERS_CSV, index=False)

print(f"\nCluster sizes:\n{clusters_df['climate_cluster'].value_counts().sort_index().to_string()}")
print(f"\nSite-year cluster assignments written to '{OUTPUT_CLUSTERS_CSV}'.")

# ---------------------------------------------------------------------------
# 4b. 2D projection for visualizing cluster separation. This is a SEPARATE,
#     purely-for-plotting PCA on the same combined_features_scaled used for
#     clustering - it compresses the (potentially dozens-of-dimensions)
#     feature space KMeans actually used down to 2 components, so the 8
#     groups can be looked at on a scatter plot. It's an approximation of
#     the true separation (some cluster structure that's clear in the full
#     feature space can look muddier in just 2D), not a re-derivation of
#     the clusters themselves.
# ---------------------------------------------------------------------------
projection_pca = PCA(n_components=2, random_state=RANDOM_STATE)
projection_2d = projection_pca.fit_transform(combined_features_scaled)
explained = projection_pca.explained_variance_ratio_
print(f"\n2D projection (for visualization only) explains "
      f"{explained[0]*100:.1f}% + {explained[1]*100:.1f}% = {explained.sum()*100:.1f}% "
      "of the clustering feature variance.")

projection_df = pd.DataFrame(kept_keys, columns=['site_id', 'year'])
projection_df['climate_cluster'] = cluster_labels
projection_df['pc1'] = projection_2d[:, 0]
projection_df['pc2'] = projection_2d[:, 1]
projection_df.to_csv(OUTPUT_CLUSTER_PROJECTION_CSV, index=False)
print(f"2D cluster projection written to '{OUTPUT_CLUSTER_PROJECTION_CSV}'.")

# ---------------------------------------------------------------------------
# 5. Per-cluster mean daily curve for each variable - lets you actually
#    interpret what each of the 8 groups' climate looks like (e.g. plot
#    these later: warm/dry vs cool/wet vs high-radiation groups, etc.)
# ---------------------------------------------------------------------------
mean_rows = []
for var in CLIMATE_VARS:
    for cluster_id in range(N_CLUSTERS):
        member_mask = cluster_labels == cluster_id
        if not member_mask.any():
            continue
        mean_curve = clean_matrices[var][member_mask].mean(axis=0)
        for doy, val in zip(day_range, mean_curve):
            mean_rows.append({'climate_cluster': cluster_id, 'variable': var, 'doy': doy, 'mean_value': val})

pd.DataFrame(mean_rows).to_csv(OUTPUT_CLUSTER_MEANS_CSV, index=False)
print(f"Per-cluster mean daily time series written to '{OUTPUT_CLUSTER_MEANS_CSV}'.")