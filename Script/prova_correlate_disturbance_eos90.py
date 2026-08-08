import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Script/ and data/ are sibling folders under the repo root, so this works
# regardless of the working directory the script is launched from.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

PHENOLOGY_CSV = DATA_DIR / "phenology_sos_eos_by_site_year_index.csv"          # prova_sigmoid_fit.py output
DAILY_ZSCORE_CSV = DATA_DIR / "fluxnet_daily_with_disturbance_zscores.csv"     # prova_zscore_disturbance.py output

OUTPUT_SITEYEAR_CSV = DATA_DIR / "disturbance_sum_vs_delta_eos90_site_year.csv"
OUTPUT_CORRELATION_CSV = DATA_DIR / "disturbance_eos90_correlation_by_site.csv"

ANOMALY_VARS = ['GPP', 'stress']  # must match the prefixes used in prova_zscore_disturbance.py

# How to define the summation window for each site-year:
#   'growing_season' - that site-year's own SOS90 -> EOS90 (from PHENOLOGY_CSV,
#                       same vi_index), i.e. the active-canopy period. This is
#                       the scientifically motivated default: it tests whether
#                       stress/GPP anomalies accrued while the canopy was
#                       active feed forward into senescence timing.
#   'fixed'           - a fixed DOY range (SUM_START_DOY to SUM_END_DOY) for
#                       every site-year, e.g. to test a specific pre-senescence
#                       window instead of the full active season.
SUM_WINDOW_MODE = 'growing_season'
SUM_START_DOY, SUM_END_DOY = 1, 365  # only used if SUM_WINDOW_MODE == 'fixed'

MIN_DAYS_IN_WINDOW = 20   # minimum daily obs required to trust a site-year's sum
MIN_YEARS_PER_SITE = 4    # minimum site-years required before correlating at a site


# ---------------------------------------------------------------------------
# 1. Load inputs
# ---------------------------------------------------------------------------
if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run prova_sigmoid_fit.py first.")
if not os.path.exists(DAILY_ZSCORE_CSV):
    raise FileNotFoundError(f"Missing '{DAILY_ZSCORE_CSV}'. Run prova_zscore_disturbance.py first.")

pheno = pd.read_csv(PHENOLOGY_CSV)
daily = pd.read_csv(DAILY_ZSCORE_CSV)
daily['date'] = pd.to_datetime(daily['date'])
if 'doy' not in daily.columns:
    daily['doy'] = daily['date'].dt.dayofyear
if 'year' not in daily.columns:
    daily['year'] = daily['date'].dt.year

pos_neg_cols = []
for v in ANOMALY_VARS:
    pos_neg_cols += [f'{v}_zscore_pos', f'{v}_zscore_neg']
pos_neg_cols = [c for c in pos_neg_cols if c in daily.columns]
if not pos_neg_cols:
    raise ValueError("None of the expected *_zscore_pos/_neg columns were found in DAILY_ZSCORE_CSV.")


# ---------------------------------------------------------------------------
# 2. Delta-EOS90 per site-year-vi_index (site's own multi-year mean as baseline)
# ---------------------------------------------------------------------------
pheno = pheno.dropna(subset=['EOS90']).copy()
pheno['EOS90_site_mean'] = pheno.groupby(['site_id', 'vi_index'])['EOS90'].transform('mean')
pheno['delta_EOS90'] = pheno['EOS90'] - pheno['EOS90_site_mean']


# ---------------------------------------------------------------------------
# 3. Sum anomalies within each site-year's window, then attach delta-EOS90
# ---------------------------------------------------------------------------
def window_bounds(row):
    if SUM_WINDOW_MODE == 'growing_season':
        start, end = row.get('SOS90'), row.get('EOS90')
    else:
        start, end = SUM_START_DOY, SUM_END_DOY
    if pd.isna(start) or pd.isna(end) or start >= end:
        return None, None
    return start, end


records = []
for _, row in pheno.iterrows():
    start, end = window_bounds(row)
    if start is None:
        continue

    site_daily = daily[(daily['site_id'] == row['site_id']) & (daily['year'] == row['year'])]
    in_window = site_daily[(site_daily['doy'] >= start) & (site_daily['doy'] <= end)]
    if len(in_window) < MIN_DAYS_IN_WINDOW:
        continue

    rec = {
        'site_id': row['site_id'], 'year': row['year'], 'vi_index': row['vi_index'],
        'window_start_doy': start, 'window_end_doy': end, 'n_days_in_window': len(in_window),
        'EOS90': row['EOS90'], 'delta_EOS90': row['delta_EOS90'],
    }
    for col in pos_neg_cols:
        rec[f'{col}_sum'] = in_window[col].sum(skipna=True)
    records.append(rec)

site_year = pd.DataFrame(records)
site_year.to_csv(OUTPUT_SITEYEAR_CSV, index=False)
print(f"Site-year table (delta-EOS90 + windowed anomaly sums) written to '{OUTPUT_SITEYEAR_CSV}' "
      f"({len(site_year)} rows).")


# ---------------------------------------------------------------------------
# 4. Per-site correlation: delta-EOS90 vs. each anomaly sum, separately per site
# ---------------------------------------------------------------------------
sum_cols = [f'{c}_sum' for c in pos_neg_cols]
corr_rows = []

for (site_id, vi_index), group in site_year.groupby(['site_id', 'vi_index']):
    group = group.dropna(subset=['delta_EOS90'])
    if len(group) < MIN_YEARS_PER_SITE:
        continue
    for col in sum_cols:
        sub = group.dropna(subset=[col])
        if len(sub) < MIN_YEARS_PER_SITE or sub[col].std() == 0:
            continue
        r, p = pearsonr(sub[col], sub['delta_EOS90'])
        corr_rows.append({
            'site_id': site_id, 'vi_index': vi_index, 'anomaly_variable': col,
            'n_years': len(sub), 'pearson_r': r, 'p_value': p
        })

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUTPUT_CORRELATION_CSV, index=False)
print(f"Per-site delta-EOS90 correlations written to '{OUTPUT_CORRELATION_CSV}' ({len(corr_df)} rows).")

if not corr_df.empty:
    n_sig = (corr_df['p_value'] < 0.05).sum()
    print(f"\n{n_sig}/{len(corr_df)} site x vi_index x anomaly-variable correlations are significant at p<0.05.")
    print("Remember: with typically ~4-8 years per site, these are exploratory, low-power correlations - "
          "treat individual-site significance cautiously and look at the sign/magnitude pattern across sites.")