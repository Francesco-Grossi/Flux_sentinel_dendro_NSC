import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # script 5 output
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                            # script 4 output

OUTPUT_GRID_CSV = DATA_DIR / "gpp_window_correlation_scan.csv"       # every window tested
OUTPUT_BEST_CSV = DATA_DIR / "gpp_window_best_correlations.csv"      # best +/- window per param x index

VI_INDICES = ['NDVI', 'NIRv']
AUTUMN_PARAMS = ['EOS90', 'EOS50', 'EOS10', 'senescence_kinetic_i']

GPP_COLUMN = 'GPP_NT_VUT_REF'
GPP_PLAUSIBLE_RANGE = (-5, 50)  # gC m-2 d-1 - guards against leftover sentinel values

# Search grid for candidate [start_doy, end_doy] windows. Restricted to the
# Mar-Oct growing season range used elsewhere in this pipeline (DOY 60-305)
# rather than the full calendar year, since a GPP window entirely outside
# the growing season isn't a meaningful "cumulative GPP" period to test.
SEARCH_START_DOY = 60
SEARCH_END_DOY = 305
WINDOW_STEP_DAYS = 5      # granularity of both the start and end grid
MIN_WINDOW_LENGTH_DAYS = 15  # shortest window tested - very short windows are noisy

MIN_COVERAGE_FRAC = 0.7   # a window's GPP sum is only trusted if at least
                           # this fraction of its days have real data
MIN_PAIRS_FOR_CORR = 10   # minimum non-NaN (window-sum, autumn-param) pairs
                           # required before a correlation is computed

if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run the double-logistic phenology script first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run the FLUXNET + Landsat merge script first.")


# ---------------------------------------------------------------------------
# 1. Load phenology results - only successful fits, only NDVI/NIRv
# ---------------------------------------------------------------------------
pheno = pd.read_csv(PHENOLOGY_CSV)
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')].copy()
pheno['year'] = pheno['year'].astype(int)

# ---------------------------------------------------------------------------
# 2. Load daily FLUXNET data and build a (site_id, year) x day-of-year GPP
#    matrix, so any window's cumulative sum can be read off with a prefix
#    sum instead of re-summing per window (there can be 1,000+ candidate
#    windows, so this matters).
# ---------------------------------------------------------------------------
flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year
flux['doy'] = flux['date'].dt.dayofyear

if GPP_COLUMN not in flux.columns:
    raise ValueError(f"'{GPP_COLUMN}' not found in FLUX_CSV.")

gpp_lo, gpp_hi = GPP_PLAUSIBLE_RANGE
implausible = (flux[GPP_COLUMN] < gpp_lo) | (flux[GPP_COLUMN] > gpp_hi)
if implausible.any():
    print(f"Excluding {int(implausible.sum())} implausible '{GPP_COLUMN}' value(s) outside "
          f"[{gpp_lo}, {gpp_hi}] gC m-2 d-1 before building the window-sum matrix.")
    flux.loc[implausible, GPP_COLUMN] = np.nan

# Pivot to one row per site-year, one column per day-of-year (1-366).
# Missing days (gaps in the record, or simply days that fall outside a
# site's actual span) come out as NaN - handled explicitly below via the
# coverage-fraction check rather than silently treated as 0.
flux_pivot = flux.pivot_table(index=['site_id', 'year'], columns='doy', values=GPP_COLUMN, aggfunc='mean')
flux_pivot = flux_pivot.reindex(columns=range(1, 367))

matrix_keys = list(flux_pivot.index)          # list of (site_id, year) tuples
matrix_pos = {key: i for i, key in enumerate(matrix_keys)}
matrix_values = flux_pivot.to_numpy(dtype=float)   # shape (n_site_years, 366)

valid_mask = ~np.isnan(matrix_values)
filled = np.nan_to_num(matrix_values, nan=0.0)

# 1-indexed prefix sums: cumsum_vals[:, d] = sum of days 1..d (d=0 -> 0).
cumsum_vals = np.concatenate([np.zeros((filled.shape[0], 1)), np.cumsum(filled, axis=1)], axis=1)
cumsum_counts = np.concatenate([np.zeros((valid_mask.shape[0], 1)),
                                 np.cumsum(valid_mask.astype(float), axis=1)], axis=1)


def window_sums_for_positions(positions, start_doy, end_doy):
    """Vectorized cumulative-GPP sum over [start_doy, end_doy] (inclusive)
    for every row in `positions` (an array of matrix row indices), with
    NaN wherever that site-year's data coverage in the window falls below
    MIN_COVERAGE_FRAC."""
    s = cumsum_vals[positions, end_doy] - cumsum_vals[positions, start_doy - 1]
    counts = cumsum_counts[positions, end_doy] - cumsum_counts[positions, start_doy - 1]
    coverage = counts / (end_doy - start_doy + 1)
    return np.where(coverage >= MIN_COVERAGE_FRAC, s, np.nan)


# ---------------------------------------------------------------------------
# 3. Build the candidate window grid
# ---------------------------------------------------------------------------
starts = list(range(SEARCH_START_DOY, SEARCH_END_DOY, WINDOW_STEP_DAYS))
windows = [
    (s, e)
    for s in starts
    for e in range(s + MIN_WINDOW_LENGTH_DAYS, SEARCH_END_DOY + 1, WINDOW_STEP_DAYS)
]
print(f"Scanning {len(windows)} candidate GPP windows "
      f"(DOY {SEARCH_START_DOY}-{SEARCH_END_DOY}, step {WINDOW_STEP_DAYS}d, "
      f"min length {MIN_WINDOW_LENGTH_DAYS}d) x {len(AUTUMN_PARAMS)} autumn parameters "
      f"x {len(VI_INDICES)} VI indices...")

# ---------------------------------------------------------------------------
# 4. Scan: for each VI index, align its site-years to the GPP matrix once,
#    then test every window against every autumn parameter.
# ---------------------------------------------------------------------------
results = []
for vi in VI_INDICES:
    sub = pheno[pheno['vi_index'] == vi].copy()
    sub['_pos'] = sub.apply(lambda r: matrix_pos.get((r['site_id'], r['year'])), axis=1)
    sub = sub.dropna(subset=['_pos']).copy()
    sub['_pos'] = sub['_pos'].astype(int)
    positions = sub['_pos'].to_numpy()

    y_by_param = {p: sub[p].to_numpy(dtype=float) for p in AUTUMN_PARAMS if p in sub.columns}

    for start_doy, end_doy in windows:
        w_sum = window_sums_for_positions(positions, start_doy, end_doy)
        for autumn_p, y in y_by_param.items():
            mask = ~np.isnan(w_sum) & ~np.isnan(y)
            n = int(mask.sum())
            if n < MIN_PAIRS_FOR_CORR or np.std(w_sum[mask]) == 0 or np.std(y[mask]) == 0:
                continue
            r, p_value = pearsonr(w_sum[mask], y[mask])
            results.append({
                'vi_index': vi, 'autumn_parameter': autumn_p,
                'window_start_doy': start_doy, 'window_end_doy': end_doy,
                'window_length_days': end_doy - start_doy + 1,
                'n': n, 'pearson_r': r, 'p_value': p_value,
            })

scan_df = pd.DataFrame(results)
scan_df.to_csv(OUTPUT_GRID_CSV, index=False)
print(f"\nFull window scan ({len(scan_df):,} rows) written to '{OUTPUT_GRID_CSV}'.")

# ---------------------------------------------------------------------------
# 5. Best positive and best negative window per (vi_index, autumn_parameter)
# ---------------------------------------------------------------------------
best_rows = []
for (vi, autumn_p), group in scan_df.groupby(['vi_index', 'autumn_parameter']):
    best_pos = group.loc[group['pearson_r'].idxmax()]
    best_neg = group.loc[group['pearson_r'].idxmin()]
    for label, row in [('most_positive', best_pos), ('most_negative', best_neg)]:
        best_rows.append({
            'vi_index': vi, 'autumn_parameter': autumn_p, 'correlation_type': label,
            'window_start_doy': int(row['window_start_doy']), 'window_end_doy': int(row['window_end_doy']),
            'window_length_days': int(row['window_length_days']),
            'n': int(row['n']), 'pearson_r': row['pearson_r'], 'p_value': row['p_value'],
        })

best_df = pd.DataFrame(best_rows)
best_df.to_csv(OUTPUT_BEST_CSV, index=False)
print(f"Best +/- window per parameter x index written to '{OUTPUT_BEST_CSV}'.")

print("\nStrongest windows found:")
print(best_df.to_string(index=False))
