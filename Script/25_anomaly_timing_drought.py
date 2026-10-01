"""
PIPELINE STEP 25 - Within a year, how much do positive / negative anomalies
of the environment and of GPP shift senescence, WHEN in the season do they
matter, and is it different in drought years?

1. Anomalies. The year is cut into WINDOW_DAYS-day windows every STEP_DAYS
   days (DOY FIRST_DOY - LAST_DOY). For each site, window and variable
       TA   air temperature (mean)        SW   incoming shortwave (mean)
       VPD  vapour pressure deficit (mean)   P    precipitation (sum)
       CWB  climatic water balance P - PET (sum), PET after Makkink (from SW, TA)
       GPP  gross primary production (mean)
   the window value of each year is turned into a z-score against the same
   window in the site's other years (>= MIN_YEARS years): the anomaly of
   that year at that time of year, in units of the site's own variability.

2. Effect on senescence, separately for positive and negative anomalies:
       EOS anomaly (days, year minus site mean) ~ pos + neg
   with pos = max(z, 0) and neg = max(-z, 0), so
       beta_pos = days of EOS shift per +1 SD POSITIVE anomaly
       beta_neg = days of EOS shift per 1 SD NEGATIVE anomaly
   (negative beta = earlier senescence). OLS on the within-site anomalies,
   standard errors clustered by site. One fit per window = the "when".

3. Drought years. A site-year is a drought year when its growing-season
   (DROUGHT_DOY) climatic water balance is <= DROUGHT_Z site standard
   deviations below the site mean. Every fit is repeated for drought and
   non-drought years, with a z-test of the difference; the mean EOS anomaly
   of drought vs non-drought years is reported as well, per leaf habit.

4. Events. For each site-year the strongest negative GPP window before the
   site's mean onset of senescence is located (timing and size), and the EOS
   anomaly is summarised by when that event happened.

Causality guard: a window is used for a site only if it ends before that
site's mean EOS90 (GPP: a late-season GPP drop IS the senescence, not its
cause) or before its mean target EOS (environment).

Input : phenology tables (steps 5 / 19 / 20), data/fluxnet_landsat_merged.csv
Output: data/anomaly_windows_by_site_year.csv     z-anomalies, long format
        data/drought_years.csv                    site-year drought classification
        data/anomaly_timing_effects.csv           every window x variable x sign x subset
        data/drought_eos_contrast.csv             EOS anomaly, drought vs non-drought years
        data/gpp_event_timing.csv                 strongest negative GPP event per site-year
        figure/anomaly_timing/<source>_<target>.png   heatmaps (when x what x sign)
        Output/anomaly_timing_summary.md
"""
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
import eos_common as ec

OUT_ANOM = ec.DATA_DIR / "anomaly_windows_by_site_year.csv"
OUT_DROUGHT = ec.DATA_DIR / "drought_years.csv"
OUT_EFFECTS = ec.DATA_DIR / "anomaly_timing_effects.csv"
OUT_CONTRAST = ec.DATA_DIR / "drought_eos_contrast.csv"
OUT_EVENTS = ec.DATA_DIR / "gpp_event_timing.csv"
OUT_MD = ec.ROOT / "Output" / "anomaly_timing_summary.md"
FIG_DIR = ec.FIGURE_DIR / "anomaly_timing"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_MD.parent.mkdir(parents=True, exist_ok=True)

WINDOW_DAYS, STEP_DAYS = 30, 15
FIRST_DOY, LAST_DOY = 61, 300
MIN_YEARS = 5                       # years per site x window to define an anomaly
MIN_COVERAGE = 0.8                  # valid days in a window
DROUGHT_DOY = (121, 273)            # May - September
DROUGHT_Z = -1.0
TARGETS = ['EOS10', 'EOS50']
VARIABLES = ['TA', 'SW', 'VPD', 'P', 'CWB', 'GPP']
SUM_VARS = {'P', 'CWB'}
SUBSETS = ['all', 'drought', 'non-drought']
MIN_OBS, MIN_SITES = 30, 5
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous', 'ENF': 'evergreen', 'EBF': 'evergreen'}
EVENT_BINS = [(61, 120, 'early spring (Mar-Apr)'), (121, 171, 'late spring (May - solstice)'),
              (172, 227, 'early summer (solstice - mid Aug)'), (228, 300, 'late summer / autumn')]


def makkink_pet(sw, ta):
    """Potential evapotranspiration (mm/day) from shortwave (W m-2) and air temperature (C)."""
    es = 0.6108 * np.exp(17.27 * ta / (ta + 237.3))
    delta = 4098.0 * es / (ta + 237.3) ** 2                       # kPa/C
    gamma, lam = 0.066, 2.45                                     # kPa/C, MJ/kg
    return np.clip(0.61 * delta / (delta + gamma) * (sw * 0.0864) / lam - 0.12, 0, None)


def bh_fdr(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    out = np.empty(len(p))
    out[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(out, 0, 1)


def cluster_ols(d, y, xs):
    X = sm.add_constant(d[xs].to_numpy(float))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return sm.OLS(d[y].to_numpy(float), X).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(d['site_id'])[0]})


# ------------------------------------------------------------------ 1. window anomalies
flux = ec.load_flux_daily()
flux['CWB'] = flux['P'] - makkink_pet(flux['SW'], flux['TA'])
habit = flux.drop_duplicates('site_id').set_index('site_id')['igbp'].map(LEAF_HABIT).fillna('grass/shrub')
lk = ec.FluxLookup(flux, VARIABLES)
site_years = flux[['site_id', 'year']].drop_duplicates()
starts = list(range(FIRST_DOY, LAST_DOY - WINDOW_DAYS + 2, STEP_DAYS))

recs = []
for s, y in zip(site_years['site_id'], site_years['year']):
    for st in starts:
        for v in VARIABLES:
            m = lk.window_mean(s, y, v, st, st + WINDOW_DAYS - 1, min_cov=MIN_COVERAGE)
            recs.append((s, int(y), st, v, m * WINDOW_DAYS if v in SUM_VARS else m))
anom = pd.DataFrame(recs, columns=['site_id', 'year', 'window_start', 'variable', 'value']).dropna(subset=['value'])
g = anom.groupby(['site_id', 'window_start', 'variable'])['value']
anom['z'] = (anom['value'] - g.transform('mean')) / g.transform('std')
anom = anom[(g.transform('count') >= MIN_YEARS) & np.isfinite(anom['z'])].copy()
anom['window_mid'] = anom['window_start'] + WINDOW_DAYS / 2
anom.to_csv(OUT_ANOM, index=False)
print(f"Window anomalies -> '{OUT_ANOM}' ({len(anom):,} rows, {anom['site_id'].nunique()} sites, "
      f"{len(starts)} windows of {WINDOW_DAYS} d)")

# ------------------------------------------------------------------ 3a. drought years
n_gs = DROUGHT_DOY[1] - DROUGHT_DOY[0] + 1
dr = site_years.copy()
dr['cwb_growing_season'] = [lk.window_mean(s, y, 'CWB', *DROUGHT_DOY, min_cov=MIN_COVERAGE) * n_gs
                            for s, y in zip(dr['site_id'], dr['year'])]
dr = dr.dropna(subset=['cwb_growing_season'])
gd = dr.groupby('site_id')['cwb_growing_season']
dr['cwb_z'] = (dr['cwb_growing_season'] - gd.transform('mean')) / gd.transform('std')
dr = dr[(gd.transform('count') >= MIN_YEARS) & np.isfinite(dr['cwb_z'])].copy()
dr['drought'] = dr['cwb_z'] <= DROUGHT_Z
dr['leaf_habit'] = dr['site_id'].map(habit)
dr.to_csv(OUT_DROUGHT, index=False)
print(f"Drought years -> '{OUT_DROUGHT}': {int(dr['drought'].sum())} of {len(dr)} site-years "
      f"({100 * dr['drought'].mean():.0f}%) at {dr['site_id'].nunique()} sites")

# ------------------------------------------------------------------ phenology: EOS anomalies
pheno = ec.load_phenology()[['site_id', 'year', 'vi_index'] + TARGETS + ['EOS90']]
pheno = pheno.merge(dr[['site_id', 'year', 'drought', 'cwb_z', 'leaf_habit']], on=['site_id', 'year'], how='inner')
for col in set(TARGETS + ['EOS90']):
    gp = pheno.groupby(['site_id', 'vi_index'])[col]
    pheno[f'{col}_sitemean'] = gp.transform('mean').where(gp.transform('count') >= 3)
    pheno[f'{col}_anom'] = pheno[col] - pheno[f'{col}_sitemean']
sources = sorted(pheno['vi_index'].unique())
wide = anom.pivot_table(index=['site_id', 'year', 'window_start'], columns='variable', values='z').reset_index()

# ------------------------------------------------------------------ 2 + 3b. effect by window, sign, subset
rows = []
for src in sources:
    ps = pheno[pheno['vi_index'] == src]
    for target in TARGETS:
        base = ps.dropna(subset=[f'{target}_anom']).merge(wide, on=['site_id', 'year'], how='inner')
        base['window_end'] = base['window_start'] + WINDOW_DAYS - 1
        for v in VARIABLES:
            limit = base['EOS90_sitemean'] if v == 'GPP' else base[f'{target}_sitemean']
            dv = base[base['window_end'] < limit].dropna(subset=[v])
            for st, dw in dv.groupby('window_start'):
                dw = dw.assign(pos=dw[v].clip(lower=0), neg=(-dw[v]).clip(lower=0))
                fits = {}
                for subset in SUBSETS:
                    d = dw if subset == 'all' else dw[dw['drought'] == (subset == 'drought')]
                    if len(d) < MIN_OBS or d['site_id'].nunique() < MIN_SITES or d['pos'].std() == 0 or d['neg'].std() == 0:
                        continue
                    f = cluster_ols(d, f'{target}_anom', ['pos', 'neg'])
                    fits[subset] = f
                    for k, sign in ((1, 'positive'), (2, 'negative')):
                        rows.append({'eos_source': src, 'target': target, 'variable': v, 'window_start': st,
                                     'window_mid': st + WINDOW_DAYS / 2, 'sign': sign, 'subset': subset,
                                     'beta_days_per_sd': f.params[k], 'std_err': f.bse[k], 'p_value': f.pvalues[k],
                                     'n_obs': len(d), 'n_sites': d['site_id'].nunique()})
                if 'drought' in fits and 'non-drought' in fits:   # drought vs non-drought difference
                    for k, sign in ((1, 'positive'), (2, 'negative')):
                        diff = fits['drought'].params[k] - fits['non-drought'].params[k]
                        se = np.hypot(fits['drought'].bse[k], fits['non-drought'].bse[k])
                        for r in rows[-6:]:
                            if r['sign'] == sign and r['subset'] == 'drought':
                                r['diff_vs_non_drought'] = diff
                                r['p_diff'] = float(2 * norm.sf(abs(diff / se))) if se > 0 else np.nan
effects = pd.DataFrame(rows)
if effects.empty:
    raise SystemExit("No window had enough data.")
effects['p_fdr'] = effects.groupby('subset')['p_value'].transform(bh_fdr)
effects.to_csv(OUT_EFFECTS, index=False)
n_raw, n_fdr = int((effects['p_value'] < 0.05).sum()), int((effects['p_fdr'] < 0.05).sum())
print(f"\nAnomaly timing effects -> '{OUT_EFFECTS}' ({len(effects)} cells; p < 0.05: {n_raw} uncorrected "
      f"= {100 * n_raw / len(effects):.1f}%, {n_fdr} after FDR)")

# ------------------------------------------------------------------ 3c. drought vs non-drought EOS anomaly
crows = []
for src in sources:
    for target in TARGETS:
        for hab in ['all', 'deciduous', 'evergreen', 'grass/shrub']:
            d = pheno[(pheno['vi_index'] == src)].dropna(subset=[f'{target}_anom'])
            d = d if hab == 'all' else d[d['leaf_habit'] == hab]
            if d['drought'].sum() < 8 or (~d['drought']).sum() < 8 or d['site_id'].nunique() < MIN_SITES:
                continue
            f = cluster_ols(d.assign(dr=d['drought'].astype(float)), f'{target}_anom', ['dr'])
            crows.append({'eos_source': src, 'target': target, 'leaf_habit': hab,
                          'n_drought': int(d['drought'].sum()), 'n_non_drought': int((~d['drought']).sum()),
                          'n_sites': d['site_id'].nunique(),
                          'eos_shift_in_drought_years_days': f.params[1], 'std_err': f.bse[1], 'p_value': f.pvalues[1]})
contrast = pd.DataFrame(crows)
contrast.to_csv(OUT_CONTRAST, index=False)
print(f"Drought vs non-drought EOS -> '{OUT_CONTRAST}'")
pd.set_option('display.width', 250)
print(contrast.round(2).to_string(index=False))

# ------------------------------------------------------------------ 4. strongest negative GPP event
gz = anom[anom['variable'] == 'GPP']
erows = []
for src in sources:
    ps = pheno[pheno['vi_index'] == src].dropna(subset=['EOS90_sitemean'])
    m = ps.merge(gz[['site_id', 'year', 'window_start', 'window_mid', 'z']], on=['site_id', 'year'])
    m = m[m['window_start'] + WINDOW_DAYS - 1 < m['EOS90_sitemean']]
    if m.empty:
        continue
    ev = m.loc[m.groupby(['site_id', 'year'])['z'].idxmin()].copy()
    ev['event_period'] = None
    for lo, hi, lab in EVENT_BINS:
        ev.loc[(ev['window_mid'] >= lo) & (ev['window_mid'] <= hi), 'event_period'] = lab
    ev['strong_event'] = ev['z'] <= -1.5
    erows.append(ev[['site_id', 'year', 'vi_index', 'leaf_habit', 'drought', 'window_mid', 'z', 'event_period',
                     'strong_event'] + [f'{t}_anom' for t in TARGETS]].rename(
        columns={'vi_index': 'eos_source', 'window_mid': 'event_doy', 'z': 'event_gpp_z'}))
events = pd.concat(erows, ignore_index=True)
events.to_csv(OUT_EVENTS, index=False)
strong = events[events['strong_event']]
ev_tab = strong.groupby(['eos_source', 'event_period'])[[f'{t}_anom' for t in TARGETS]].agg(['mean', 'count']).round(1)
print(f"\nStrongest negative GPP event per site-year -> '{OUT_EVENTS}'")
print("Mean EOS anomaly (days) in years with a strong negative GPP event (z <= -1.5), by when it happened:")
print(ev_tab.to_string())

# ------------------------------------------------------------------ figures
row_labels = [f'{v} {s}' for v in VARIABLES for s in ('+', '-')]
for src in sources:
    for target in TARGETS:
        sub = effects[(effects['eos_source'] == src) & (effects['target'] == target)]
        if sub.empty:
            continue
        fig, axes = plt.subplots(1, len(SUBSETS), figsize=(5.2 * len(SUBSETS), 5.2), sharey=True)
        vmax = max(5.0, float(np.nanpercentile(np.abs(sub['beta_days_per_sd']), 95)))
        for ax, subset in zip(axes, SUBSETS):
            ss = sub[sub['subset'] == subset].assign(
                row=lambda d: d['variable'] + ' ' + d['sign'].map({'positive': '+', 'negative': '-'}))
            mat = ss.pivot(index='row', columns='window_mid', values='beta_days_per_sd').reindex(row_labels) \
                .reindex(columns=[s + WINDOW_DAYS / 2 for s in starts])
            sig = ss.pivot(index='row', columns='window_mid', values='p_fdr').reindex(row_labels) \
                .reindex(columns=mat.columns)
            raw = ss.pivot(index='row', columns='window_mid', values='p_value').reindex(row_labels) \
                .reindex(columns=mat.columns)
            im = ax.imshow(mat.to_numpy(float), cmap='RdBu', vmin=-vmax, vmax=vmax, aspect='auto')
            for i in range(mat.shape[0]):
                for j in range(mat.shape[1]):
                    if np.isfinite(sig.iloc[i, j]) and sig.iloc[i, j] < 0.05:
                        ax.text(j, i, '**', ha='center', va='center', fontsize=9)
                    elif np.isfinite(raw.iloc[i, j]) and raw.iloc[i, j] < 0.05:
                        ax.text(j, i, '*', ha='center', va='center', fontsize=9, color='0.3')
            ax.set_xticks(range(len(mat.columns)), [int(c) for c in mat.columns], rotation=90, fontsize=7)
            ax.set_yticks(range(len(row_labels)), row_labels, fontsize=8)
            n = ss['n_obs'].max() if len(ss) else 0
            ax.set_title(f"{subset} years (up to {int(n) if np.isfinite(n) else 0} site-years)", fontsize=9)
            ax.set_xlabel('window centre (day of year)')
        fig.colorbar(im, ax=axes, label='EOS shift (days per 1 SD anomaly; red = earlier)', shrink=0.8)
        fig.suptitle(f"{src} {target}: effect of positive (+) / negative (-) anomalies, by timing   "
                     f"(* p<0.05, ** p<0.05 after FDR)", fontsize=10)
        fig.savefig(FIG_DIR / f"{src}_{target}.png", dpi=150, bbox_inches='tight')
        plt.close(fig)

# ------------------------------------------------------------------ text summary
lines = ["# Anomaly timing and drought - effect on senescence", "",
         "Generated by `Script/25_anomaly_timing_drought.py`. beta = days of EOS shift per 1 SD anomaly "
         "(negative = earlier senescence), separately for positive and negative anomalies.", "",
         f"- Drought year: May-September climatic water balance (P - PET) <= {DROUGHT_Z} site SD: "
         f"{int(dr['drought'].sum())} of {len(dr)} site-years ({100 * dr['drought'].mean():.0f}%).",
         f"- Cells tested: {len(effects)}; uncorrected p < 0.05 in {100 * n_raw / len(effects):.1f}% "
         f"(about 5% expected by chance); {n_fdr} after FDR.", "",
         "## Mean EOS shift in drought years vs non-drought years (days)", "",
         "| EOS source | target | leaf habit | drought / non-drought site-years | shift | p |", "|---|---|---|---|---|---|"]
for _, r in contrast.iterrows():
    lines.append(f"| {r['eos_source']} | {r['target']} | {r['leaf_habit']} | {r['n_drought']} / {r['n_non_drought']} | "
                 f"{r['eos_shift_in_drought_years_days']:+.1f} | {r['p_value']:.3f} |")
lines += ["", "## Effects significant after FDR", ""]
sig = effects[effects['p_fdr'] < 0.05].sort_values(['subset', 'eos_source', 'target', 'variable', 'window_start'])
if sig.empty:
    lines.append("None.")
else:
    lines += ["| years | EOS source | target | variable | anomaly | window (DOY) | beta | p_fdr | n |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in sig.iterrows():
        lines.append(f"| {r['subset']} | {r['eos_source']} | {r['target']} | {r['variable']} | {r['sign']} | "
                     f"{int(r['window_start'])}-{int(r['window_start']) + WINDOW_DAYS - 1} | "
                     f"{r['beta_days_per_sd']:+.1f} | {r['p_fdr']:.3f} | {r['n_obs']} |")
lines += ["", "## Consistent across EOS sources",
          "Cells (all years) where every EOS source with data has the same sign and at least two have "
          "uncorrected p < 0.05.", ""]
allc = effects[effects['subset'] == 'all']
cons = allc.groupby(['target', 'variable', 'sign', 'window_start']).agg(
    n_sources=('beta_days_per_sd', 'size'), n_neg=('beta_days_per_sd', lambda b: int((b < 0).sum())),
    n_p05=('p_value', lambda p: int((p < 0.05).sum())), median_beta=('beta_days_per_sd', 'median')).reset_index()
cons = cons[(cons['n_sources'] >= 3) & cons['n_neg'].isin([0]) | (cons['n_sources'] >= 3) & (cons['n_neg'] == cons['n_sources'])]
cons = cons[cons['n_p05'] >= 2]
if cons.empty:
    lines.append("None.")
else:
    lines += ["| target | variable | anomaly | window (DOY) | sources | median beta |", "|---|---|---|---|---|---|"]
    for _, r in cons.iterrows():
        lines.append(f"| {r['target']} | {r['variable']} | {r['sign']} | {int(r['window_start'])}-"
                     f"{int(r['window_start']) + WINDOW_DAYS - 1} | {r['n_sources']} ({r['n_p05']} with p<0.05) | "
                     f"{r['median_beta']:+.1f} |")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print(f"\nConsistent across sources (all years): {len(cons)} cells")
if len(cons):
    print(cons.round(1).to_string(index=False))
print(f"\nSummary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")
