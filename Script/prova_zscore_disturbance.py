import os
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

FLUX_CSV = DATA_DIR / "fluxnet_daily_selected_vars.csv"
# Output stays at daily resolution: the original FLUXNET table with the new
# z-score columns appended (no site-year aggregation/summing here - that
# happens downstream, e.g. per-window in whatever correlates against EOS10).
OUTPUT_DAILY_CSV = DATA_DIR / "fluxnet_daily_with_disturbance_zscores.csv"

# Which FLUXNET columns count as "flux" (GPP) and "stress" for this analysis.
# Add more entries here (e.g. 'GPP_DT': 'GPP_DT_VUT_REF', 'temp_stress':
# 'TA_F') if you want additional disturbance variables run through the same
# moving-window z-score.
VARS_FOR_ZSCORE = {
    'GPP': 'GPP_NT_VUT_REF',
    'stress': 'VPD_F',   # vapor pressure deficit as the water-stress proxy
}

WINDOW_DAYS = 7          # +/- DOY window, pooled across ALL years available at a site
MIN_OBS_PER_WINDOW = 10  # minimum pooled (non-NaN) obs required to trust a window's mean/std

# prova_fluxent_data.py's own comment says GPP_NT/DT_VUT_REF and RECO_NT_VUT_REF
# inherit NEE_VUT_REF's QC flag (they're outputs of the NEE partitioning step
# and carry no QC flag of their own) - but its QC_PAIRS dict only actually
# blanks NEE_VUT_REF itself, not these three. Re-applying it here closes that
# gap so the z-scores aren't built on low-confidence gap-filled GPP values.
QC_THRESHOLD = 2
QC_INHERIT_FROM_NEE = ['GPP_NT_VUT_REF', 'GPP_DT_VUT_REF', 'RECO_NT_VUT_REF']

DOY_MOD = 366  # circular wraparound base for the window (Dec-Jan boundary)


# ---------------------------------------------------------------------------
# 1. Load data, re-apply the missing QC inheritance, derive DOY/year
# ---------------------------------------------------------------------------
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run prova_fluxent_data.py first.")

df = pd.read_csv(FLUX_CSV)
df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
df['year'] = df['date'].dt.year
df['doy'] = df['date'].dt.dayofyear

if 'NEE_VUT_REF_QC' in df.columns:
    qc_numeric = pd.to_numeric(df['NEE_VUT_REF_QC'], errors='coerce')
    bad = qc_numeric > QC_THRESHOLD
    for col in QC_INHERIT_FROM_NEE:
        if col in df.columns:
            df.loc[bad, col] = np.nan
else:
    print("Note: 'NEE_VUT_REF_QC' not found - skipping the GPP/RECO QC-inheritance re-application.")


# ---------------------------------------------------------------------------
# 2. Circular +/- WINDOW_DAYS DOY-window climatology, pooled across years
# ---------------------------------------------------------------------------
def circular_distance(doy_a, doy_b, mod=DOY_MOD):
    """Shortest distance between two DOYs on a circular year (handles the
    Dec 31 / Jan 1 wraparound so a window centered near day 1 or day 365
    still pools the days on the other side of the boundary)."""
    diff = np.abs(doy_a - doy_b) % mod
    return np.minimum(diff, mod - diff)


def window_climatology(doy_values, window=WINDOW_DAYS, min_obs=MIN_OBS_PER_WINDOW):
    """Given a Series of (doy -> value) observations pooled across all years
    at one site, return per-target-DOY (mean, std, n) for the +/-window
    circular DOY neighborhood, pooling every year's observations that fall
    in that neighborhood."""
    doy_arr = doy_values.index.to_numpy()
    val_arr = doy_values.to_numpy()

    stats = {}
    for target_doy in range(1, DOY_MOD + 1):
        mask = circular_distance(doy_arr, target_doy) <= window
        pool = val_arr[mask]
        pool = pool[~np.isnan(pool)]
        if len(pool) >= min_obs:
            stats[target_doy] = (np.mean(pool), np.std(pool, ddof=1), len(pool))
        else:
            stats[target_doy] = (np.nan, np.nan, len(pool))
    return stats


def add_zscore_for_site(df_site, value_col, out_prefix):
    """Compute the DOY-window z-score for one variable at one site, and
    split it into positive-only / negative-only columns (0.0 where the
    z-score doesn't apply, so both sums are non-cancelling and directly
    summable)."""
    doy_values = df_site.set_index('doy')[value_col].dropna()
    doy_values = doy_values[doy_values.index.notna()]

    stats = window_climatology(doy_values)

    means = df_site['doy'].map(lambda d: stats.get(int(d), (np.nan, np.nan, 0))[0])
    stds = df_site['doy'].map(lambda d: stats.get(int(d), (np.nan, np.nan, 0))[1])

    z = (df_site[value_col] - means) / stds
    z = z.replace([np.inf, -np.inf], np.nan)

    df_site[f'{out_prefix}_zscore'] = z
    df_site[f'{out_prefix}_zscore_pos'] = z.where(z > 0, 0.0)
    df_site[f'{out_prefix}_zscore_neg'] = z.where(z < 0, 0.0)
    # keep NaN (undefined) days as NaN in the raw column but 0.0 in the
    # pos/neg split so nansum-style aggregation below isn't affected
    undefined = z.isna()
    df_site.loc[undefined, [f'{out_prefix}_zscore_pos', f'{out_prefix}_zscore_neg']] = np.nan
    return df_site


# ---------------------------------------------------------------------------
# 3. Apply per site, per variable
# ---------------------------------------------------------------------------
site_frames = []
total_sites = df['site_id'].nunique()
for i, (site_id, df_site) in enumerate(df.groupby('site_id'), start=1):
    df_site = df_site.copy()
    for prefix, col in VARS_FOR_ZSCORE.items():
        if col not in df_site.columns:
            print(f"   -> {site_id}: column '{col}' not found, skipping '{prefix}'")
            continue
        df_site = add_zscore_for_site(df_site, col, prefix)
    site_frames.append(df_site)
    if i % 25 == 0 or i == total_sites:
        print(f"Processed {i}/{total_sites} sites")

df_out = pd.concat(site_frames, ignore_index=True)
df_out.to_csv(OUTPUT_DAILY_CSV, index=False)

added_cols = []
for prefix in VARS_FOR_ZSCORE:
    added_cols += [f'{prefix}_zscore', f'{prefix}_zscore_pos', f'{prefix}_zscore_neg']
added_cols = [c for c in added_cols if c in df_out.columns]

print(f"\nDaily data with disturbance z-score columns written to '{OUTPUT_DAILY_CSV}'.")
print(f"New columns added (one row per site-date, same granularity as the input): {added_cols}")