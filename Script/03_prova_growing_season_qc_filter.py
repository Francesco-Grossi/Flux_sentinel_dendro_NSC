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

INPUT_CSV = DATA_DIR / "fluxnet_daily_all_vars.csv"   # output of the FLUXNET extraction script

OUTPUT_CSV = DATA_DIR / "fluxnet_daily_all_vars_qc_filtered.csv"
SUMMARY_CSV = DATA_DIR / "site_year_growing_season_qc_summary.csv"

# "Biological" (carbon-flux) variables, as opposed to the environmental
# ones (TA_F, VPD_F, SW_IN_F, P_F) and the energy-flux one (LE_F_MDS) -
# this matches the project description's own split ("mainly source (GPP,
# NEE, NPP) and environmental conditions (Tair, VPD, precipitation...)").
# Each entry maps the value column to its QC column. Note that GPP and RECO
# don't carry an independent QC flag in the FLUXNET daily product - both
# columns were written out as copies of NEE_VUT_REF_QC upstream (see the
# FLUXNET extraction script) - so in practice all three checks below key
# off the same QC series. They're kept as three separate checks anyway so
# this still does the right thing if a future data source supplies GPP/RECO
# with their own independent QC.
BIOLOGICAL_VARS = {
    'GPP_NT_VUT_REF': 'GPP_NT_VUT_REF_QC',
    'NEE_VUT_REF': 'NEE_VUT_REF_QC',
    'RECO_NT_VUT_REF': 'RECO_NT_VUT_REF_QC',
}

# Growing season window: March 1 - October 31 inclusive, every year.
GS_START_MONTH, GS_START_DAY = 3, 1
GS_END_MONTH, GS_END_DAY = 10, 31

# FLUXNET QC convention: 0=measured, 1=good gap-fill, 2=medium, 3=poor.
# This filter uses a stricter cut than the QC_THRESHOLD=2 used upstream to
# blank individual bad values (which only nulls out QC=3/poor) - here,
# QC=2 (medium) counts as "low quality" too, since this is a coarser
# site-year reliability gate rather than a per-value cleaning step. A day
# with no row at all (missing from the record entirely) or a missing/NaN
# QC value is also treated as low-quality/bad.
LOW_QUALITY_QC_THRESHOLD = 1  # QC > 1 (i.e. 2 or 3) counts as low quality

MAX_BAD_FRACTION = 0.5        # exclude site-year if >= 50% of growing-season days are low quality
MAX_CONSECUTIVE_BAD_DAYS = 50  # exclude site-year if any run of >= 50 consecutive bad days


# ---------------------------------------------------------------------------
# 1. Load + parse dates
# ---------------------------------------------------------------------------
if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing '{INPUT_CSV}'. Run the FLUXNET extraction script first.")

df = pd.read_csv(INPUT_CSV)
df['date'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
df['date'] = df['date'].fillna(pd.to_datetime(df['TIMESTAMP'].astype(str), errors='coerce'))
df = df.dropna(subset=['date']).copy()
df['year'] = df['date'].dt.year

available_bio_vars = {v: qc for v, qc in BIOLOGICAL_VARS.items() if v in df.columns}
missing_bio_vars = set(BIOLOGICAL_VARS) - set(available_bio_vars)
if missing_bio_vars:
    print(f"Note: biological variable(s) not found in input, skipped: {sorted(missing_bio_vars)}")


def longest_true_run(bool_array):
    """Length of the longest run of consecutive True values."""
    max_run = 0
    cur = 0
    for b in bool_array:
        if b:
            cur += 1
            if cur > max_run:
                max_run = cur
        else:
            cur = 0
    return max_run


# ---------------------------------------------------------------------------
# 2. Per site-year growing-season QC check
# ---------------------------------------------------------------------------
# Each site-year's growing-season window is reindexed onto the FULL Mar 1 -
# Oct 31 calendar range for that year, so a calendar day that's entirely
# missing from the record (e.g. inside a data gap, or before/after the
# site's actual record starts/ends) counts as a bad day too - not just days
# that are present but flagged low-quality. This is what makes the
# "50 consecutive bad days" check meaningful even across a gap in the
# underlying record, not just within the rows that happen to exist.
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

n_total = len(summary_df)
n_passed = int(summary_df['passed'].sum())
print(f"Site-year growing-season QC check: {n_passed}/{n_total} site-years passed "
      f"(< {int(MAX_BAD_FRACTION * 100)}% low-quality days and no run of "
      f">= {MAX_CONSECUTIVE_BAD_DAYS} consecutive bad days, in all of {list(available_bio_vars)}).")
print(f"Per-site-year QC summary written to '{SUMMARY_CSV}'.")

# ---------------------------------------------------------------------------
# 3. Filter the full daily dataset (all months, not just growing season) to
#    only the site-years that passed - a failing growing season disqualifies
#    the whole year's record for that site, since it's the reliability of
#    the site-year as a unit being assessed here, not individual days.
# ---------------------------------------------------------------------------
passing_df = pd.DataFrame(list(passing_site_years), columns=['site_id', 'year'])
filtered_df = df.merge(passing_df, on=['site_id', 'year'], how='inner').drop(columns=['date', 'year'])

filtered_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nFiltered daily data written to '{OUTPUT_CSV}' "
      f"({len(filtered_df):,} rows, {filtered_df['site_id'].nunique()} sites, "
      f"{n_passed} passing site-years).")
