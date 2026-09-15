import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, linregress

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
# Answers: (1) does the RATE of GPP increase during spring green-up affect
# autumn senescence timing and progression, and (2) does the RATE of
# temperature/photoperiod decline during autumn affect how fast senescence
# itself proceeds (senescence_kinetic_i)? These are rate-of-change (slope)
# predictors, distinct from the level/cumulative predictors in the main
# autumn-phenology-correlation script.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"  # script 5 output
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"  # script 4 output (full daily FLUXNET + sparse VI)

OUTPUT_SITEYEAR_CSV = DATA_DIR / "rate_of_change_predictors_by_site_year_index.csv"
OUTPUT_CORR_CSV = DATA_DIR / "rate_of_change_correlations.csv"

VI_INDICES = ['NDVI', 'NIRv']
AUTUMN_PARAMS = ['EOS90', 'EOS50', 'EOS10', 'senescence_kinetic_i']

MIN_PAIRS_FOR_CORR = 5
# Minimum number of daily points required to fit a rate-of-change (slope)
# over a window - a slope from 2-3 noisy daily values isn't trustworthy.
MIN_POINTS_FOR_RATE = 5
# Quality gate on the Landsat phenology curve fit itself - see the main
# autumn-phenology-correlation script for the full rationale.
MIN_FIT_CORR = 0.8

# Length of the pre-EOS90 window (days) used for gpp_pre_eos90_rate below -
# the 15 days immediately before, and including, the EOS90 crossing itself.
PRE_EOS90_WINDOW_DAYS = 15

# GPP_NT_VUT_REF: nighttime-partitioning GPP (see the main correlation
# script for why NT over DT). TA_F: air temperature. Plausible ranges guard
# against leftover -9999-style sentinel values.
FLUX_VARS = {
    'gpp': {'column': 'GPP_NT_VUT_REF', 'plausible_range': (-5, 50)},
    'temperature': {'column': 'TA_F', 'plausible_range': (-60, 50)},
}

if not os.path.exists(PHENOLOGY_CSV):
    raise FileNotFoundError(f"Missing '{PHENOLOGY_CSV}'. Run the double-logistic phenology script first.")
if not os.path.exists(FLUX_CSV):
    raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run the FLUXNET + Landsat merge script first.")


# ---------------------------------------------------------------------------
# 1. Load phenology results - only successful fits, only NDVI/NIRv
# ---------------------------------------------------------------------------
pheno = pd.read_csv(PHENOLOGY_CSV)
pheno = pheno[pheno['vi_index'].isin(VI_INDICES) & (pheno['method'] == 'double_logistic')
              & (pheno['corr'] >= MIN_FIT_CORR)].copy()
print(f"Phenology rows available (NDVI/NIRv, successful fits, corr >= {MIN_FIT_CORR}): {len(pheno)}")

# ---------------------------------------------------------------------------
# 2. Load daily FLUXNET data, clean GPP/temperature, index by (site_id,
#    year) for fast repeated window lookups.
# ---------------------------------------------------------------------------
flux = pd.read_csv(FLUX_CSV)
flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
flux['date'] = flux['date'].fillna(pd.to_datetime(flux['TIMESTAMP'].astype(str), errors='coerce'))
flux = flux.dropna(subset=['date']).copy()
flux['year'] = flux['date'].dt.year
flux['doy'] = flux['date'].dt.dayofyear

flux_lookups = {}
for var_name, spec in FLUX_VARS.items():
    col, (lo, hi) = spec['column'], spec['plausible_range']
    if col not in flux.columns:
        raise ValueError(f"'{col}' not found in FLUX_CSV.")
    implausible = (flux[col] < lo) | (flux[col] > hi)
    if implausible.any():
        print(f"Excluding {int(implausible.sum())} implausible '{col}' value(s) outside [{lo}, {hi}].")
        flux.loc[implausible, col] = np.nan
    flux_lookups[var_name] = {
        key: group[['doy', col]].dropna(subset=['doy'])
        for key, group in flux.groupby(['site_id', 'year'])
    }

site_lat = flux.groupby('site_id')['lat'].first().to_dict()


def rate_flux_var(var_name, site_id, year, doy_start, doy_end):
    """RATE of change (slope, units/day) of a daily flux variable over
    [doy_start, doy_end] - fit by linear regression of the variable against
    DOY, not its level/mean. NaN if the window is degenerate, undefined, or
    has too few points (< MIN_POINTS_FOR_RATE) to fit a meaningful slope."""
    if pd.isna(doy_start) or pd.isna(doy_end) or doy_start >= doy_end:
        return np.nan
    g = flux_lookups[var_name].get((site_id, year))
    if g is None or g.empty:
        return np.nan
    col = FLUX_VARS[var_name]['column']
    mask = (g['doy'] >= doy_start) & (g['doy'] <= doy_end)
    sub = g.loc[mask, ['doy', col]].dropna()
    if len(sub) < MIN_POINTS_FOR_RATE or sub['doy'].nunique() < 2:
        return np.nan
    slope, intercept, r_value, p_value, std_err = linregress(sub['doy'], sub[col])
    return float(slope)


def photoperiod_hours(doy, lat_deg):
    """Day length (hours) at a given latitude and day-of-year, via the
    standard solar-declination approximation."""
    if pd.isna(doy) or pd.isna(lat_deg):
        return np.nan
    lat_rad = np.radians(lat_deg)
    declination_deg = -23.44 * np.cos(np.radians(360.0 / 365.0 * (doy + 10)))
    declination_rad = np.radians(declination_deg)
    cos_hour_angle = -np.tan(lat_rad) * np.tan(declination_rad)
    cos_hour_angle = np.clip(cos_hour_angle, -1.0, 1.0)  # polar day/night
    return float((24.0 / np.pi) * np.arccos(cos_hour_angle))


def photoperiod_rate(doy_start, doy_end, lat_deg):
    """RATE of change of day length (hours/day) over [doy_start, doy_end].
    Photoperiod is an exact, noise-free function of DOY and latitude, so
    the two-point slope across the window IS the rate, not a regression fit
    through noisy observations (unlike rate_flux_var)."""
    if pd.isna(doy_start) or pd.isna(doy_end) or pd.isna(lat_deg) or doy_start >= doy_end:
        return np.nan
    p_start = photoperiod_hours(doy_start, lat_deg)
    p_end = photoperiod_hours(doy_end, lat_deg)
    if pd.isna(p_start) or pd.isna(p_end):
        return np.nan
    return (p_end - p_start) / (doy_end - doy_start)


# ---------------------------------------------------------------------------
# 3. Build the per-site-year-index rate-predictor table
# ---------------------------------------------------------------------------
# gpp_greenup_rate: fit over the green-up window (SOS10->SOS90) - the
# spring GPP ramp-up rate.
# temperature/photoperiod_senescence_rate: fit over the senescence window
# (EOS90->EOS10, onset-of-decline to near-dormancy) - the autumn
# cooling/daylight-loss rate, specifically during the decline itself.
# gpp_pre_eos90_rate: fit over the PRE_EOS90_WINDOW_DAYS days immediately
# before EOS90 (inclusive). NOTE: like photoperiod_at_EOS90 elsewhere in
# this pipeline, this window's END boundary is anchored to EOS90's own DOY,
# so correlating it against EOS90 itself isn't testing an independent
# driver in the strictest sense - different site-years get different
# windows depending on their own EOS90 timing. It's still informative
# (a real, computed carbon-flux trend, not a deterministic function of
# EOS90 the way photoperiod is), but keep that anchoring in mind when
# interpreting the EOS90 correlation specifically; the correlations against
# EOS50, EOS10, and senescence_kinetic_i are more clearly independent.
records = []
for _, row in pheno.iterrows():
    site_id, year = row['site_id'], int(row['year'])
    sos10, sos90 = row.get('leaf_out_10'), row.get('leaf_out_90')
    eos90, eos10 = row.get('EOS90'), row.get('EOS10')
    lat = site_lat.get(site_id)

    rec = {
        'site_id': site_id, 'year': year, 'vi_index': row['vi_index'],
        'EOS90': eos90, 'EOS50': row.get('EOS50'), 'EOS10': eos10,
        'senescence_kinetic_i': row.get('senescence_kinetic_i'),
        'gpp_greenup_rate': rate_flux_var('gpp', site_id, year, sos10, sos90),
        'temperature_senescence_rate': rate_flux_var('temperature', site_id, year, eos90, eos10),
        'photoperiod_senescence_rate': photoperiod_rate(eos90, eos10, lat),
        # GPP trend in the 15 days leading up to, and including, EOS90 -
        # tests whether carbon-uptake trajectory right at the onset of
        # decline (rather than the whole prior season) relates to how/when
        # senescence subsequently plays out.
        'gpp_pre_eos90_rate': (rate_flux_var('gpp', site_id, year, eos90 - PRE_EOS90_WINDOW_DAYS, eos90)
                                if pd.notna(eos90) else np.nan),
    }
    records.append(rec)

analysis_df = pd.DataFrame(records)
analysis_df.to_csv(OUTPUT_SITEYEAR_CSV, index=False)
print(f"Rate-predictor table written to '{OUTPUT_SITEYEAR_CSV}' ({len(analysis_df)} rows).")

# ---------------------------------------------------------------------------
# 4. Correlate each autumn parameter against each rate predictor, per VI
#    index. gpp_greenup_rate is tested against ALL 4 autumn parameters
#    (timing AND progression); temperature/photoperiod_senescence_rate are
#    the ones specifically motivated by senescence_kinetic_i, but are still
#    tested against all 4 for completeness.
# ---------------------------------------------------------------------------
RATE_PREDICTORS = ['gpp_greenup_rate', 'temperature_senescence_rate', 'photoperiod_senescence_rate',
                    'gpp_pre_eos90_rate']

corr_rows = []
for vi in VI_INDICES:
    sub = analysis_df[analysis_df['vi_index'] == vi]
    for autumn_p in AUTUMN_PARAMS:
        for pred in RATE_PREDICTORS:
            pair = sub[[pred, autumn_p]].dropna()
            if len(pair) < MIN_PAIRS_FOR_CORR or pair[pred].std() == 0 or pair[autumn_p].std() == 0:
                continue
            r, p_value = pearsonr(pair[pred], pair[autumn_p])
            corr_rows.append({
                'vi_index': vi, 'autumn_parameter': autumn_p, 'predictor': pred,
                'n': len(pair), 'pearson_r': r, 'p_value': p_value,
            })

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUTPUT_CORR_CSV, index=False)
print(f"Rate-of-change correlation table written to '{OUTPUT_CORR_CSV}' ({len(corr_df)} rows).")

if not corr_df.empty:
    n_sig = int((corr_df['p_value'] < 0.05).sum())
    print(f"\n{n_sig}/{len(corr_df)} correlations significant at p<0.05.")
    print("\nAll rate-of-change correlations:")
    print(corr_df.sort_values(['predictor', 'autumn_parameter'])
          [['vi_index', 'autumn_parameter', 'predictor', 'n', 'pearson_r', 'p_value']].to_string(index=False))