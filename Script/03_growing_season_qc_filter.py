"""
PIPELINE STEP 3 - Growing-season (Mar 1 - Oct 31) QC filter, per site-year.

A site-year is dropped if >=50% of its growing-season days are low-quality
(QC>1) or NaN for GPP/NEE/RECO, or if any run of >=50 consecutive bad days
occurs. A failing growing season disqualifies that site-year's whole daily
record (all months), not just the growing-season window.

Output: data/fluxnet_daily_all_vars_qc_filtered.csv
        data/site_year_growing_season_qc_summary.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
INPUT_CSV = DATA_DIR / "fluxnet_daily_all_vars.csv"
OUTPUT_CSV = DATA_DIR / "fluxnet_daily_all_vars_qc_filtered.csv"
SUMMARY_CSV = DATA_DIR / "site_year_growing_season_qc_summary.csv"

BIOLOGICAL_VARS = {
    'GPP_NT_VUT_REF': 'GPP_NT_VUT_REF_QC',
    'NEE_VUT_REF': 'NEE_VUT_REF_QC',
    'RECO_NT_VUT_REF': 'RECO_NT_VUT_REF_QC',
}
GS_START_MONTH, GS_START_DAY = 3, 1
GS_END_MONTH, GS_END_DAY = 10, 31
LOW_QUALITY_QC_THRESHOLD = 1
MAX_BAD_FRACTION = 0.5
MAX_CONSECUTIVE_BAD_DAYS = 50

if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run 01_fluxnet_download.py first.")

df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
df = df.dropna(subset=['date']).copy()
df['year'] = df['date'].dt.year

available_bio_vars = {v: qc for v, qc in BIOLOGICAL_VARS.items() if v in df.columns}


def longest_true_run(bool_array):
    max_run = cur = 0
    for b in bool_array:
        if b:
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 0
    return max_run


summary_rows = []
passing_site_years = set()

for site_id, site_df in df.groupby('site_id'):
    site_df = site_df.sort_values('date').set_index('date')
    min_year, max_year = site_df.index.min().year, site_df.index.max().year

    for year in range(min_year, max_year + 1):
        gs_start = pd.Timestamp(year=year, month=GS_START_MONTH, day=GS_START_DAY)
        gs_end = pd.Timestamp(year=year, month=GS_END_MONTH, day=GS_END_DAY)
        full_range = pd.date_range(gs_start, gs_end, freq='D')
        sub = site_df.reindex(full_range)

        fail_vars = []
        var_details = {}
        for var, qc_col in available_bio_vars.items():
            if qc_col in sub.columns:
                qc_numeric = pd.to_numeric(sub[qc_col], errors='coerce')
                bad = qc_numeric.isna() | (qc_numeric > LOW_QUALITY_QC_THRESHOLD)
            else:
                bad = pd.Series(True, index=full_range)
            frac_bad = float(bad.mean())
            max_consecutive_bad = longest_true_run(bad.to_numpy())
            var_fail = (frac_bad >= MAX_BAD_FRACTION) or (max_consecutive_bad >= MAX_CONSECUTIVE_BAD_DAYS)
            var_details[var] = (frac_bad, max_consecutive_bad, var_fail)
            if var_fail:
                fail_vars.append(var)

        passed = len(fail_vars) == 0
        if passed:
            passing_site_years.add((site_id, year))

        row = {'site_id': site_id, 'year': year, 'n_growing_season_days': len(full_range), 'passed': passed,
               'failed_variables': ';'.join(fail_vars) if fail_vars else ''}
        for var, (frac_bad, max_consecutive_bad, var_fail) in var_details.items():
            row[f'{var}_frac_bad'] = round(frac_bad, 3)
            row[f'{var}_max_consecutive_bad_days'] = max_consecutive_bad
        summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(SUMMARY_CSV, index=False)
n_total, n_passed = len(summary_df), int(summary_df['passed'].sum())
print(f"Site-year QC: {n_passed}/{n_total} passed. Summary -> '{SUMMARY_CSV}'.")

passing_df = pd.DataFrame(list(passing_site_years), columns=['site_id', 'year'])
filtered_df = df.merge(passing_df, on=['site_id', 'year'], how='inner').drop(columns=['date', 'year'])
filtered_df.to_csv(OUTPUT_CSV, index=False)
print(f"Filtered daily data -> '{OUTPUT_CSV}' ({len(filtered_df):,} rows, {n_passed} passing site-years).")
