"""
PIPELINE STEP 42 - Build one Markdown report with the results of every
analysis step and the corresponding figures: Output/results_report.md.

The report only reads the tables and figures written by steps 21-41, so it
is always in line with the last pipeline run. Tables are cut to the rows
that matter (the full tables are the CSV files named in each section);
figures are linked, not copied (paths relative to Output/).

Input : data/*.csv and figure/**/*.png of steps 21-41
Output: Output/results_report.md
"""
import os
from datetime import date

import numpy as np
import pandas as pd
import eos_common as ec

OUT_MD = ec.ROOT / "Output" / "results_report.md"
OUT_MD.parent.mkdir(parents=True, exist_ok=True)
D, F = ec.DATA_DIR, ec.FIGURE_DIR
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous', 'ENF': 'evergreen', 'EBF': 'evergreen'}
SRC_ORDER = ['GCC', 'NDVI_tower', 'NDVI', 'NIRv', 'GPP']
SRC_NAME = {'GCC': 'PhenoCam (GCC)', 'NDVI_tower': 'Tower NDVI', 'NDVI': 'Satellite NDVI', 'NIRv': 'Satellite NIRv',
            'GPP': 'GPP-derived'}
lines = []


def read(name):
    p = D / name
    return pd.read_csv(p) if os.path.exists(p) else None


def add(*text):
    lines.extend(text)


def table(df, fmt=None, max_rows=60):
    """Markdown table; fmt maps column -> format string."""
    if df is None or len(df) == 0:
        add("_No data._", "")
        return
    fmt = fmt or {}
    df = df.head(max_rows)
    add("| " + " | ".join(str(c) for c in df.columns) + " |", "|" + "---|" * len(df.columns))
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append("-" if pd.isna(v) else format(v, fmt.get(c, '.2f')))
            else:
                cells.append(str(v))
        add("| " + " | ".join(cells) + " |")
    add("")


def figure(rel, caption):
    if os.path.exists(F / rel):
        add(f"![{caption}](../figure/{rel})", "", f"*{caption}*", "")
    else:
        add(f"_Figure not found: figure/{rel}_", "")


def order_src(df, col):
    return df.assign(_o=df[col].map({s: i for i, s in enumerate(SRC_ORDER)})).sort_values('_o').drop(columns='_o')


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''


add("# Phenology - carbon pipeline: results", "",
    f"Generated on {date.today().isoformat()} by `Script/42_results_report.py` from the outputs of the last pipeline run.",
    "", "Effects are days of shift in the end of season (EOS) per +1 SD of the predictor unless stated otherwise; "
    "negative = earlier senescence. EOS90 / EOS50 / EOS10 = day of year when greenness has fallen to 90 / 50 / 10 % "
    "of its seasonal amplitude (onset, middle, end of senescence).", "")

# =============================================================== 1. data
add("## 1. Data", "")
merged = pd.read_csv(ec.FLUX_CSV, usecols=['site_id', 'igbp', 'TIMESTAMP'])
merged['year'] = merged['TIMESTAMP'] // 10000
sy = merged.drop_duplicates(['site_id', 'year'])
qc = read("site_year_growing_season_qc_summary.csv")
add(f"- Flux sites passing quality control: **{sy['site_id'].nunique()}**, with **{len(sy)}** site-years "
    f"({int(sy['year'].min())}-{int(sy['year'].max())})"
    + (f"; {int(qc['passed'].sum())} of {len(qc)} site-years passed the growing-season QC." if qc is not None else "."))
veg = sy.drop_duplicates('site_id')['igbp'].value_counts()
add("- Sites by vegetation type (IGBP): " + ", ".join(f"{k} {v}" for k, v in veg.items()) + ".")
hls = read("fluxnet_all_highlat_landsat_indices.csv")
if hls is not None:
    hls['year'] = pd.to_datetime(hls['date']).dt.year
    n_clear = hls.dropna(subset=['NDVI']).groupby(['site_id', 'year']).size().median()
    add(f"- Satellite (HLS): {len(hls):,} image records at {hls['site_id'].nunique()} sites "
        f"({hls['sensor'].value_counts().to_dict()}); median {n_clear:.0f} clear images per site-year.")
season = read("site_season_type.csv")
if season is not None:
    add(f"- Season type: {season['season_type'].value_counts().to_dict()} "
        "(dry-summer = seasonal GPP peaks before 1 June).")
add("")

# =============================================================== 2. EOS sources
add("## 2. End-of-season dates from four independent sources", "",
    "Steps 23-25 fit the same curve to satellite NDVI / NIRv, tower broadband NDVI and PhenoCam greenness. "
    "Tables: `data/eos_source_summary.csv`, `data/eos_source_agreement.csv`.", "")
summ = read("eos_source_summary.csv")
if summ is not None:
    t = summ.pivot(index='source', columns='level', values='median_within_site_sd_days').reindex(columns=['EOS90', 'EOS50', 'EOS10'])
    n = summ[summ['level'] == 'EOS10'].set_index('source')[['n_site_years', 'n_sites', 'first_year', 'last_year']]
    med = summ.pivot(index='source', columns='level', values='median_doy').reindex(columns=['EOS90', 'EOS50', 'EOS10'])
    t = n.join(med.add_prefix('median DOY ')).join(t.add_prefix('within-site SD ')).reset_index()
    t = order_src(t, 'source')
    t['source'] = t['source'].map(SRC_NAME)
    add("**Coverage, typical dates and year-to-year variation within a site (days)**", "")
    table(t, {c: '.0f' for c in t.columns if 'median' in c or c in ('n_site_years', 'n_sites', 'first_year', 'last_year')}
          | {c: '.1f' for c in t.columns if 'SD' in c})
agree = read("eos_source_agreement.csv")
if agree is not None:
    a = agree[agree['leaf_habit'] == 'all']
    t = a.pivot_table(index=['source_a', 'source_b'], columns='level', values='r_within').reindex(columns=['EOS90', 'EOS50', 'EOS10'])
    t = t.join(a[a['level'] == 'EOS10'].set_index(['source_a', 'source_b'])[['n_site_years', 'n_sites']]).reset_index()
    t['pair'] = t['source_a'].map(SRC_NAME) + " vs " + t['source_b'].map(SRC_NAME)
    add("**Do the sources agree on early and late years? Within-site correlation of EOS anomalies**", "")
    table(t[['pair', 'EOS90', 'EOS50', 'EOS10', 'n_site_years', 'n_sites']], {'n_site_years': '.0f', 'n_sites': '.0f'})
figure("eos_source_agreement/r_within_EOS50.png", "Year-to-year agreement between EOS sources, EOS50 (step 31)")
figure("eos_source_agreement/anomaly_scatter_vs_GCC.png", "EOS anomalies of each source against PhenoCam (step 31)")
figure("vi_eos_vs_gpp_eos/EOS10.png", "EOS10 of each source against the GPP-derived EOS10 (step 30)")

# =============================================================== 3. main hypothesis
add("## 3. Are Zani and Lu both right? The split-GPP hypothesis (step 41)", "",
    "- **H1** cumulative GPP from leaf-out (SOS10) to the summer solstice shifts EOS **earlier** (Zani et al. 2020)",
    "- **H2** cumulative GPP from the solstice to EOS shifts EOS **later**",
    "- **H3** the two cancel, so GPP over the whole season has **no effect** (Lu et al. 2022)", "",
    "Headline model: strictly within sites (each year minus its site's mean; standard errors clustered by site). "
    "'Controlled' adds the year's leaf-out date and the air temperature before and after the solstice. "
    "Table: `data/split_gpp_cancellation_test.csv`; written summary: `Output/split_gpp_cancellation_summary.md`.", "")
h = read("split_gpp_cancellation_test.csv")
if h is not None:
    j = h[h['model'] == 'joint']
    cnt = j.groupby(['target', 'anchor']).agg(tests=('H1_pre_negative', 'size'), H1=('H1_pre_negative', 'sum'),
                                              H2=('H2_post_positive', 'sum'), H3=('H3_total_no_effect', 'sum'),
                                              all_three=('hypothesis_supported', 'sum')).reset_index()
    cnt['window version'] = cnt['anchor'].map({'fixed': 'ends at site mean EOS (unbiased)',
                                               'fixed_before_senescence': 'ends at site mean EOS90, before senescence',
                                               'year': 'ends at same-year EOS (window-length effect built in)'})
    add("**In how many tests does each part hold?** (EOS sources x site groups)", "")
    table(cnt[['target', 'window version', 'tests', 'H1', 'H2', 'H3', 'all_three']].astype({c: int for c in ['tests', 'H1', 'H2', 'H3', 'all_three']}))

    def h_table(target, anchor, groups):
        d = j[(j['target'] == target) & (j['anchor'] == anchor) & j['group'].isin(groups)].copy()
        d['group'] = pd.Categorical(d['group'], groups, ordered=True)
        d = order_src(d, 'eos_source').sort_values('group', kind='stable')
        out = pd.DataFrame({
            'sites': d['group'].astype(str), 'EOS source': d['eos_source'].map(SRC_NAME),
            'n (sites)': [f"{int(a)} ({int(b)})" for a, b in zip(d['n_within'], d['n_sites_within'])],
            'pre-solstice GPP': [f"{b:+.1f}{stars(p)}" for b, p in zip(d['b_pre_days_per_sd'], d['p_pre_negative'])],
            'post-solstice GPP': [f"{b:+.1f}{stars(p)}" for b, p in zip(d['b_post_days_per_sd'], d['p_post_positive'])],
            'whole season [95% CI]': [f"{b:+.1f} [{lo:+.1f}, {hi:+.1f}]" for b, lo, hi in
                                      zip(d['b_total_days_per_sd'], d['total_ci_lo'], d['total_ci_hi'])],
            'all three hold': np.where(d['hypothesis_supported'].astype(bool), 'yes', 'no'),
            'controlled: pre / post / whole': [
                "-" if pd.isna(a) else f"{a:+.1f}{stars(pa)} / {b:+.1f}{stars(pb)} / {c:+.1f}"
                for a, pa, b, pb, c in zip(d.get('ctrl_b_pre_days_per_sd'), d.get('ctrl_p_pre_negative'),
                                           d.get('ctrl_b_post_days_per_sd'), d.get('ctrl_p_post_positive'),
                                           d.get('ctrl_b_total_days_per_sd'))]})
        table(out)

    groups = ['all', 'summer-green', 'dry-summer', 'deciduous', 'evergreen', 'grass/shrub']
    add("### 3.1 Onset of senescence (EOS90), windows ending at the site mean EOS90", "",
        "Days per +1 within-site SD of cumulative GPP. Stars: one-sided test in the direction of the hypothesis "
        "(* p<0.05, ** p<0.01, *** p<0.001).", "")
    h_table('EOS90', 'fixed', groups)
    figure("split_gpp_cancellation/coefficients_EOS90_all.png", "EOS90, all sites: pre-solstice, post-solstice and whole-season GPP, for each window version")
    figure("split_gpp_cancellation/coefficients_EOS90_deciduous.png", "EOS90, deciduous sites")
    for src in ['NDVI', 'NIRv', 'GCC', 'NDVI_tower']:
        figure(f"split_gpp_cancellation/scatter_EOS90_{src}.png",
               f"{SRC_NAME[src]}: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)")
    add("### 3.2 End of senescence (EOS10), windows ending at the site mean EOS10", "")
    h_table('EOS10', 'fixed', groups)
    figure("split_gpp_cancellation/coefficients_EOS10_all.png", "EOS10, all sites")
    add("### 3.3 Same-year windows and the window-length effect", "",
        "When the post-solstice window ends at the same year's EOS, a later EOS makes the window longer and its "
        "cumulative GPP larger by construction. The last column is the post-solstice slope produced by window length "
        "alone (the site's average GPP curve, no year-specific GPP).", "")
    d = order_src(j[(j['target'] == 'EOS90') & (j['anchor'] == 'year') & (j['group'] == 'all')], 'eos_source')
    cols = {'eos_source': 'EOS source', 'b_pre_days_per_sd': 'pre', 'b_post_days_per_sd': 'post',
            'b_total_days_per_sd': 'whole season', 'b_post_length_only_days_per_sd': 'post from window length alone'}
    t = d[[c for c in cols if c in d.columns]].rename(columns=cols)
    t['EOS source'] = t['EOS source'].map(SRC_NAME)
    table(t, {c: '+.1f' for c in t.columns if c != 'EOS source'})

# =============================================================== 4. 1:1 plots
add("## 4. One-to-one plots: pooled vs within-site (step 40)", "",
    "The same pairs drawn twice. Pooled plots mix differences between sites with year-to-year changes; the "
    "within-site plots keep only the latter and correspond to the tests of section 3. "
    "All pairs: `figure/predictor_correlations/` and `figure/predictor_correlations_within_site/`.", "")
corr = read("autumn_phenology_correlations.csv")
if corr is not None:
    keep = ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90', 'gpp_solstice_to_eos90_fixed', 'total_gpp_growing_season']
    t = corr[(corr['autumn_parameter'] == 'EOS90') & corr['predictor'].isin(keep)] \
        .pivot(index='predictor', columns='vi_index', values='pearson_r').reindex(keep).reset_index()
    add("**Pooled correlation with EOS90**", "")
    table(t)
for pred, lab in (('gpp_sos10_to_solstice', 'GPP from SOS10 to the solstice'),
                  ('gpp_solstice_to_eos90', 'GPP from the solstice to the same-year EOS90'),
                  ('gpp_solstice_to_eos90_fixed', 'GPP from the solstice to the site mean EOS90'),
                  ('total_gpp_growing_season', 'GPP over the growing season')):
    figure(f"predictor_correlations/EOS90_vs_{pred}.png", f"Pooled: EOS90 vs {lab}")
    figure(f"predictor_correlations_within_site/EOS90_vs_{pred}.png", f"Within sites: EOS90 vs {lab}")

# =============================================================== 5. sink: CUE and NPP
add("## 5. Carbon use efficiency and NPP (steps 26, 32)", "",
    "CUE after Luo et al. (annual), extended to 30-day sliding windows (seasonal / daily). "
    "NPP = CUE x GPP. Tables: `data/cue_luo2025_site_year.csv`, `data/cue_seasonal_means_by_site_year.csv`.", "")
cue = read("cue_luo2025_site_year.csv")
means = read("cue_seasonal_means_by_site_year.csv")
if cue is not None:
    add(f"- Annual CUE: {int(cue['CUE'].notna().sum())} site-years at {cue.loc[cue['CUE'].notna(), 'site_id'].nunique()} sites, "
        f"mean {cue['CUE'].mean():.2f} (SD {cue['CUE'].std():.2f}).")
if means is not None:
    add(f"- Mean of the daily CUE before the solstice {means['CUE_pre'].mean():.2f}, after {means['CUE_post'].mean():.2f}.", "")
ca = read("cue_scale_agreement.csv")
if ca is not None:
    t = ca[ca['leaf_habit'] == 'all'][['a', 'b', 'n_site_years', 'r_within', 'r_pooled']]
    add("**How much do the time scales share? (within-site correlation)**", "")
    table(t, {'n_site_years': '.0f'})
figure("cue_seasonal/seasonal_curve.png", "Seasonal course of CUE by plant type (step 32)")
figure("cue_seasonal/annual_vs_seasonal.png", "Annual CUE against pre- and post-solstice CUE, within sites (step 32)")

# =============================================================== 6. window scan
add("## 6. Which carbon window relates to which EOS? (step 33)", "",
    "Single-predictor mixed models, windows ending at the site mean EOS (fixed anchors). "
    "Table: `data/eos_window_scan.csv`.", "")
scan = read("eos_window_scan.csv")
if scan is not None:
    s = scan[(scan['anchor'] == 'fixed') & scan['window'].isin(['SOS_to_SOL', 'SOL_to_EOS10'])].copy()
    s['cell'] = [f"{b:+.1f}{stars(p)}" for b, p in zip(s['beta_days_per_sd'], s['p_value'])]
    s['predictor'] = s['carbon'] + ' ' + s['metric']
    want = ['GPP cum', 'NPP cum', 'NPPd cum', 'GPP mean', 'CUEd mean']
    for target in ['EOS10', 'EOS50']:
        for win, lab in (('SOS_to_SOL', 'leaf-out to solstice'), ('SOL_to_EOS10', 'solstice to mean EOS10')):
            t = s[(s['target'] == target) & (s['window'] == win)].pivot(index='predictor', columns='vi_index', values='cell') \
                .reindex(want).reindex(columns=[c for c in SRC_ORDER if c in s['vi_index'].unique()]).reset_index()
            t.columns = [SRC_NAME.get(c, c) for c in t.columns]
            add(f"**{target}, window {lab}** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)", "")
            table(t.fillna('-'))
for src in ['GCC', 'NDVI']:
    figure(f"eos_window_scan/fixed_{src}_GPP.png", f"{SRC_NAME[src]}: GPP windows against EOS, fixed anchors (step 33)")
    figure(f"eos_window_scan/fixed_{src}_NPPd.png", f"{SRC_NAME[src]}: NPP (daily CUE) windows against EOS, fixed anchors (step 33)")

# =============================================================== 7. sliding scan
add("## 7. When around the solstice is the relation strongest? (step 34)", "",
    "Mean flux in 15- and 30-day windows starting 75 days before to 75 days after the solstice. "
    "Table: `data/eos_solstice_sliding_scan.csv`; the most negative window per case:", "")
mins = read("eos_solstice_sliding_scan_minima.csv")
if mins is not None:
    t = mins[(mins['length_days'] == 30) & (mins['target'] == 'EOS10')][
        ['carbon', 'vi_index', 'offset_from_solstice_days', 'beta_days_per_sd', 'p_value', 'n_obs']]
    t = t.rename(columns={'vi_index': 'EOS source', 'offset_from_solstice_days': 'window start (days from solstice)'})
    t['EOS source'] = t['EOS source'].map(SRC_NAME)
    table(t, {'window start (days from solstice)': '.0f', 'p_value': '.3f', 'n_obs': '.0f', 'beta_days_per_sd': '+.1f'})
for cv in ['GPP', 'NPPd', 'CUEd']:
    figure(f"eos_solstice_sliding_scan/{cv}_L30.png", f"{cv}: effect on EOS of 30-day windows by start date relative to the solstice (step 34)")

# =============================================================== 8. env vs carbon
add("## 8. Does carbon explain EOS beyond climate? (step 35)", "",
    "Nested mixed models on identical rows. M0 climate only; + leaf-out date (SOS); + source (GPP); + sink (NPP). "
    "Tables: `data/eos_env_vs_carbon_comparison.csv`, `..._cv.csv`, `..._lrt.csv`.", "")
cmp_, cv = read("eos_env_vs_carbon_comparison.csv"), read("eos_env_vs_carbon_cv.csv")
if cmp_ is not None and cv is not None:
    m = cmp_.merge(cv, on=['vi_index', 'target', 'model'])
    for target in ['EOS10', 'EOS50']:
        t = m[m['target'] == target][['vi_index', 'model', 'n_obs', 'aic', 'r2_marginal', 'delta_r2_marginal_vs_env', 'loso_r2']]
        t = order_src(t, 'vi_index').rename(columns={'vi_index': 'EOS source', 'r2_marginal': 'R2 (fixed effects)',
                                                      'delta_r2_marginal_vs_env': 'gain over climate',
                                                      'loso_r2': 'R2 leave-one-site-out'})
        t['EOS source'] = t['EOS source'].map(SRC_NAME)
        add(f"**{target}**", "")
        table(t, {'n_obs': '.0f', 'aic': '.0f'}, max_rows=80)

# =============================================================== 9. path analysis, rate vs cumulative
add("## 9. Pathways and rate vs cumulative uptake (steps 36, 37)", "")
for src in ['GCC', 'NDVI']:
    figure(f"eos_path_analysis/{src}_EOS10_GPP.png", f"{SRC_NAME[src]}: climate -> GPP -> EOS10 path coefficients (step 36)")
rc = read("eos_rate_vs_cumulative.csv")
if rc is not None:
    t = rc[(rc['anchor'] == 'fixed') & (rc['group'].str.upper() == 'ALL') & (rc['window'] == 'SOS_to_SOL') & (rc['carbon'] == 'GPP')]
    cols = [c for c in ['vi_index', 'target', 'n_obs', 'beta_rate', 'p_rate', 'beta_cum', 'p_cum'] if c in t.columns]
    t = order_src(t[cols], 'vi_index').rename(columns={'vi_index': 'EOS source'})
    t['EOS source'] = t['EOS source'].map(SRC_NAME)
    add("**GPP rate (mean) vs cumulative GPP, leaf-out to solstice, all sites** (`data/eos_rate_vs_cumulative.csv`)", "")
    table(t, {'n_obs': '.0f', 'beta_rate': '+.1f', 'beta_cum': '+.1f', 'p_rate': '.3f', 'p_cum': '.3f'})

# =============================================================== 10. plant type
add("## 10. Results by plant type, with multiple-testing correction (step 38)", "",
    "Tables: `data/eos_results_by_leaf_habit.csv`, `data/eos_results_consistency.csv`.", "")
cells, cons = read("eos_results_by_leaf_habit.csv"), read("eos_results_consistency.csv")
if cells is not None:
    n_raw, n_fdr = int((cells['p_value'] < 0.05).sum()), int((cells['p_fdr'] < 0.05).sum())
    add(f"- {len(cells)} tests; {n_raw} ({100 * n_raw / len(cells):.1f}%) with uncorrected p < 0.05 (about 5% expected by chance); "
        f"{n_fdr} significant after FDR correction.", "")
if cons is not None:
    rb = cons[cons['robust'].astype(bool)][['leaf_habit', 'target', 'model', 'window', 'predictor', 'direction', 'median_beta', 'sources']]
    add("**Effects with the same sign in every EOS source, both ground sources agreeing, at least one significant after FDR**", "")
    table(rb, {'median_beta': '+.1f'})
figure("eos_results_by_leaf_habit/EOS10_raw.png", "EOS10: every predictor by plant type and EOS source (step 38)")
figure("eos_results_by_leaf_habit/EOS50_adj_T.png", "EOS50, air temperature controlled (step 38)")

# =============================================================== 11. anomalies and drought
add("## 11. Timing of anomalies and drought years (step 39)", "",
    "Tables: `data/anomaly_timing_effects.csv`, `data/drought_years.csv`, `data/drought_eos_contrast.csv`.", "")
dr, dc, eff = read("drought_years.csv"), read("drought_eos_contrast.csv"), read("anomaly_timing_effects.csv")
if dr is not None:
    add(f"- Drought years (May-September water balance at least 1 SD below the site mean): {int(dr['drought'].sum())} "
        f"of {len(dr)} site-years ({100 * dr['drought'].mean():.0f}%).", "")
if dc is not None:
    t = dc[dc['leaf_habit'].isin(['all', 'deciduous'])].copy()
    t['drought / other site-years'] = t['n_drought'].astype(str) + " / " + t['n_non_drought'].astype(str)
    t = order_src(t, 'eos_source').rename(columns={'eos_source': 'EOS source', 'leaf_habit': 'sites',
                                                    'eos_shift_in_drought_years_days': 'EOS shift in drought years (days)'})
    t['EOS source'] = t['EOS source'].map(SRC_NAME)
    add("**Mean EOS shift in drought years against the other years**", "")
    table(t[['EOS source', 'target', 'sites', 'drought / other site-years', 'EOS shift in drought years (days)', 'p_value']],
          {'EOS shift in drought years (days)': '+.1f', 'p_value': '.3f'})
if eff is not None:
    n_raw = int((eff['p_value'] < 0.05).sum())
    add(f"- Anomaly windows: {len(eff)} tests, {100 * n_raw / len(eff):.1f}% with uncorrected p < 0.05; after FDR: "
        + ", ".join(f"{k} years {v}" for k, v in eff[eff['p_fdr'] < 0.05].groupby('subset').size().items()) + ".", "")
for src in ['NDVI', 'GCC']:
    figure(f"anomaly_timing/{src}_EOS10.png", f"{SRC_NAME[src]}: effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)")

add("---", "", "Per-step logs are in `logs/`. Method notes are in `README.md` and in the docstring of each script.", "")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
n_fig = sum(1 for line in lines if line.startswith("!["))
print(f"Report -> '{OUT_MD}' ({len(lines)} lines, {n_fig} figures, "
      f"{sum(1 for line in lines if line.startswith('_Figure not found'))} figures missing)")
