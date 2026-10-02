"""
PIPELINE STEP 41 - Are Zani et al. (2020) and Lu et al. (2022) both right?
Formal test of the split-GPP cancellation hypothesis, for EOS90 and EOS10:

    H1  cumulative GPP from SOS10 to the summer solstice (pre) shifts EOS EARLIER
        (Zani: more early-season productivity -> earlier senescence)
    H2  cumulative GPP from the solstice to EOS (post) shifts EOS LATER
    H3  the two cancel, so cumulative GPP over the whole season (total = pre +
        post, SOS10 -> EOS) has NO effect on EOS
        (Lu: growing-season productivity does not advance senescence)

HEADLINE MODEL - strictly within sites. Every variable is the year's value
minus the site's own mean (sites with >= MIN_YEARS_WITHIN years), predictors
are scaled to their within-site SD, and standard errors are clustered by site:
    joint     dEOS ~ dpre + dpost         -> b_pre, b_post          (H1, H2)
    total     dEOS ~ dtotal               -> b_total                (H3)
    controls  the same two models with the year's leaf-out date (SOS10) and
              the mean air temperature before and after the solstice added
              (ctrl_ columns): a warm, early spring raises pre-solstice GPP
              and could advance senescence by itself.
Slopes are days of EOS shift per +1 within-site SD of GPP. The hypothesis is
about year-to-year changes at a site, and only this model isolates them: a
pooled scatter plot, and to a lesser degree a random-intercept mixed model
(lme_ columns, kept for comparison; GPP in units of 100 gC m-2), also pick up
differences BETWEEN sites (productive sites vs unproductive ones, dry-summer
vs summer-green), which inflates the slopes.
H3 needs more than "b_total is not significant" (absence of evidence): it is
tested as EQUIVALENCE (two one-sided tests) - b_total must lie significantly
inside +/- SESOI_DAYS_PER_SD. The hypothesis is called supported only if H1,
H2 and H3 all hold. Plain correlations are also reported pooled, between
sites (site means) and within sites.

Window versions ("anchor"):
    year    post and total end at the SAME YEAR's EOS, as the hypothesis is
            worded. But then a later EOS makes the window longer, and a longer
            window accumulates more GPP whatever the vegetation does: a
            positive post slope is built in.
    fixed   post and total end at the site's MEAN EOS (same end date every
            year). No built-in window-length relation.
    fixed_before_senescence   (EOS10 target only)
            post and total end at the site's MEAN EOS90, the usual onset of
            the decline. The fixed EOS10 window still covers the weeks in
            which the canopy is senescing, and in a year when senescence is
            late the canopy keeps photosynthesising in those weeks - late EOS
            then CAUSES high post-solstice GPP, not the reverse. Ending the
            window before the decline starts removes that reverse path.
    calendar   pre = the CAL_PRE days before the solstice, post = the CAL_POST
            days from the solstice, total = both. The other versions start
            at the year's leaf-out, so an early spring lengthens the pre
            window: pre-solstice GPP then partly IS leaf-out date. Calendar
            windows depend on neither leaf-out nor EOS.
    rate    mean daily GPP over the 'fixed' windows instead of the sum:
            uptake per day, free of window length (but still averaged from
            the year's leaf-out).
How large the built-in part is, is measured with a LENGTH-ONLY NULL: the
year-anchored windows are filled with the site's average seasonal GPP curve
instead of that year's GPP, so the predictor varies ONLY through window
length. The slope it produces (b_post_length_only) is pure artefact; H2 under
the year anchor is real only to the extent that b_post exceeds it.

Site groups:
    season type   summer-green (seasonal GPP peak on or after SPRING_PEAK_DOY)
                  vs dry-summer (GPP peaks in spring, canopy senesces in early
                  summer when the soil dries out: Californian grassland and oak
                  savanna, sagebrush steppe). In dry-summer sites early
                  senescence is driven by water, and EOS can fall before the
                  solstice, in which case the post window does not exist.
    plant type    deciduous / evergreen / grass-shrub

Input : data/eos_window_predictors_{fixed,year}_anchor.csv (step 27), daily flux
Output: data/split_gpp_cancellation_test.csv
        data/site_season_type.csv
        figure/split_gpp_cancellation/coefficients_<target>_<group>.png
        figure/split_gpp_cancellation/scatter_<target>_<source>.png   (within-site anomalies)
        Output/split_gpp_cancellation_summary.md
"""
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
import eos_common as ec

OUT_CSV = ec.DATA_DIR / "split_gpp_cancellation_test.csv"
OUT_SEASON = ec.DATA_DIR / "site_season_type.csv"
OUT_MD = ec.ROOT / "Output" / "split_gpp_cancellation_summary.md"
FIG_DIR = ec.FIGURE_DIR / "split_gpp_cancellation"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_MD.parent.mkdir(parents=True, exist_ok=True)
ec.require(ec.WINDOW_FIXED_CSV, ec.WINDOW_YEAR_CSV, hint="Run 27_window_predictors.py first.")

TARGETS = ['EOS90', 'EOS10']
PRE = 'GPP_cum__SOS_to_SOL'
UNIT = 100.0                    # gC m-2 per unit of the predictors
SESOI_DAYS_PER_SD = 2.0         # smallest effect of interest for "no effect" (H3)
ALPHA = 0.05
MIN_OBS, MIN_SITES = 30, 5
MIN_YEARS_WITHIN = 3
COVARS = ['SOS', 'TA_pre', 'TA_post']     # leaf-out date, air temperature before / after the solstice
CAL_PRE, CAL_POST = 60, 45      # calendar windows: days before / from the solstice
SPRING_PEAK_DOY = 152           # seasonal GPP peak before 1 June -> dry-summer site
GROUPS = ['all', 'summer-green', 'dry-summer', 'deciduous', 'evergreen', 'grass/shrub']
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous', 'ENF': 'evergreen', 'EBF': 'evergreen'}


def variants(target):
    """(anchor name, description, window file, pre column, post column, total column)."""
    v = [('fixed', f'windows end at the site mean {target}', ec.WINDOW_FIXED_CSV, PRE,
          f'GPP_cum__SOL_to_{target}', f'GPP_cum__SOS_to_{target}')]
    if target == 'EOS10':
        v.append(('fixed_before_senescence', 'windows end at the site mean EOS90 (before the decline starts)',
                  ec.WINDOW_FIXED_CSV, PRE, 'GPP_cum__SOL_to_EOS90', 'GPP_cum__SOS_to_EOS90'))
    v.append(('calendar', f'calendar windows: {CAL_PRE} d before / {CAL_POST} d from the solstice (independent of leaf-out)',
              ec.WINDOW_FIXED_CSV, 'GPP_cal__pre', 'GPP_cal__post', 'GPP_cal__total'))
    v.append(('rate', f'mean daily GPP (rate) over the windows ending at the site mean {target}',
              ec.WINDOW_FIXED_CSV, 'GPP_mean__SOS_to_SOL', f'GPP_mean__SOL_to_{target}', f'GPP_mean__SOS_to_{target}'))
    v.append(('year', f'windows end at the same-year {target} (as worded; positive post slope built in)',
              ec.WINDOW_YEAR_CSV, PRE, f'GPP_cum__SOL_to_{target}', f'GPP_cum__SOS_to_{target}'))
    return v


def lme(d, formula):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            return smf.mixedlm(formula, d, groups=d['site_id']).fit(reml=True)
        except Exception:
            return None


def fe(fit, name):
    return float(fit.fe_params[name]), float(fit.bse_fe[name])


def within_fit(d, y, xs):
    """OLS on within-site anomalies (year minus site mean), predictors scaled to their
    within-site SD, SE clustered by site. Returns {x: (beta, se, p two-sided)}, n, n_sites."""
    return ec.within_fit(d, y, xs, MIN_YEARS_WITHIN, MIN_OBS, MIN_SITES)


def one_sided(b, se, negative):
    return float(norm.cdf(b / se)) if negative else float(norm.sf(b / se))


def corr3(d, x, y):
    """Pooled, between-site and within-site correlation."""
    out = {'pooled': np.nan, 'between': np.nan, 'within': np.nan}
    if len(d) > 2 and d[x].std() > 0:
        out['pooled'] = float(np.corrcoef(d[x], d[y])[0, 1])
    m = d.groupby('site_id')[[x, y]].mean()
    if len(m) > 2 and m[x].std() > 0 and m[y].std() > 0:
        out['between'] = float(np.corrcoef(m[x], m[y])[0, 1])
    dm = d[[x, y]] - d.groupby('site_id')[[x, y]].transform('mean')
    if len(dm) > 2 and dm[x].std() > 0 and dm[y].std() > 0:
        out['within'] = float(np.corrcoef(dm[x], dm[y])[0, 1])
    return out


# ---------------------------------------------------------------- site season type, length-only null
flux = ec.load_flux_daily()
clim = flux.groupby(['site_id', 'doy'])['GPP'].mean()            # site's average seasonal GPP curve
clim_cum, season_rows = {}, []
for s, g in clim.groupby(level=0):
    curve = g.droplevel(0).reindex(range(1, 367))
    clim_cum[s] = curve.fillna(0).cumsum()
    smooth = curve.rolling(31, center=True, min_periods=10).mean()
    if smooth.notna().any():
        peak = int(smooth.idxmax())
        season_rows.append({'site_id': s, 'gpp_peak_doy': peak,
                            'season_type': 'dry-summer' if peak < SPRING_PEAK_DOY else 'summer-green'})
season = pd.DataFrame(season_rows)
igbp = flux.drop_duplicates('site_id')[['site_id', 'igbp']]
season = season.merge(igbp, on='site_id', how='left')
season.to_csv(OUT_SEASON, index=False)
print(f"Site season type -> '{OUT_SEASON}': {season['season_type'].value_counts().to_dict()}")
print("  dry-summer sites:", ' '.join(f"{r.site_id}({r.igbp},{r.gpp_peak_doy})"
                                      for r in season[season['season_type'] == 'dry-summer'].itertuples()))
season_type = season.set_index('site_id')['season_type']


def clim_sum(site, start, end):
    """Cumulative climatological GPP over [start, end] (DOY) - depends only on the window."""
    if site not in clim_cum or not (np.isfinite(start) and np.isfinite(end)) or end < start:
        return np.nan
    c = clim_cum[site]
    s, e = int(np.ceil(start)), int(np.floor(end))
    return float(c.loc[min(e, 366)] - (c.loc[s - 1] if s > 1 else 0.0))


def in_group(w, group):
    if group == 'all':
        return w
    col = 'season_type' if group in ('summer-green', 'dry-summer') else 'leaf_habit'
    return w[w[col] == group]


lk = ec.FluxLookup(flux, ['GPP'])
csv_cache, tables, anchors = {}, {}, {}
for target in TARGETS:
    for anchor, desc, path, pre_col, post_col, total_col in variants(target):
        if path not in csv_cache:
            c = pd.read_csv(path)
            c['GPP_cal__pre'] = [lk.window_mean(s, y, 'GPP', so - CAL_PRE, so - 1) * CAL_PRE
                                 for s, y, so in zip(c['site_id'], c['year'], c['SOL'])]
            c['GPP_cal__post'] = [lk.window_mean(s, y, 'GPP', so, so + CAL_POST - 1) * CAL_POST
                                  for s, y, so in zip(c['site_id'], c['year'], c['SOL'])]
            c['GPP_cal__total'] = c['GPP_cal__pre'] + c['GPP_cal__post']
            csv_cache[path] = c
        w = csv_cache[path].copy()
        w['leaf_habit'] = w['igbp'].map(LEAF_HABIT).fillna('grass/shrub')
        w['season_type'] = w['site_id'].map(season_type)
        w['pre'], w['post'], w['total'] = w[pre_col] / UNIT, w[post_col] / UNIT, w[total_col] / UNIT
        w['TA_pre'], w['TA_post'] = w.get('TA_mean__SOS_to_SOL', np.nan), w.get('TA_mean__SOL_to_EOS10', np.nan)
        if anchor == 'year':
            w['post_len'] = [clim_sum(s, a, b) / UNIT for s, a, b in zip(w['site_id'], w['SOL'], w[f'{target}_anchor'])]
            w['pre_len'] = [clim_sum(s, a, b) / UNIT for s, a, b in zip(w['site_id'], w['SOS'], w['SOL'])]
        tables[(target, anchor)] = w
        anchors[(target, anchor)] = desc
sources = sorted(csv_cache[ec.WINDOW_FIXED_CSV]['vi_index'].unique())

# ---------------------------------------------------------------- tests
rows = []
for (target, anchor), w in tables.items():
    for group in GROUPS:
        wg = in_group(w, group)
        for src in sources:
            ws = wg[wg['vi_index'] == src]
            base = {'target': target, 'anchor': anchor, 'group': group, 'eos_source': src}
            # correlations use every row available for that predictor (as in the 1:1 plots)
            for term in ('pre', 'post', 'total'):
                c = corr3(ws[[target, 'site_id', term]].dropna(), term, target)
                base.update({f'r_{term}_{k}': v for k, v in c.items()})
            cols = [target, 'site_id', 'pre', 'post', 'total'] + (['pre_len', 'post_len'] if anchor == 'year' else [])
            d = ws[cols].dropna().reset_index(drop=True)
            base.update({'n_obs': len(d), 'n_sites': d['site_id'].nunique(),
                         'n_obs_pre_only': int(ws[[target, 'pre']].dropna().shape[0])})
            if len(d) < MIN_OBS or d['site_id'].nunique() < MIN_SITES:
                wf = within_fit(ws, target, ['pre'])         # e.g. dry-summer: EOS before the solstice, no post window
                if wf is not None:
                    (b, se, _), n_w, ns_w = wf[0]['pre'], wf[1], wf[2]
                    base.update({'model': 'pre only (no post window)', 'n_within': n_w, 'n_sites_within': ns_w,
                                 'b_pre_days_per_sd': b, 'se_pre_days_per_sd': se,
                                 'p_pre_negative': one_sided(b, se, True), 'H1_pre_negative': one_sided(b, se, True) < ALPHA})
                    rows.append(base)
                continue
            joint, tot = lme(d, f'{target} ~ pre + post'), lme(d, f'{target} ~ total')
            if joint is None or tot is None:
                continue
            b_pre, se_pre = fe(joint, 'pre')
            b_post, se_post = fe(joint, 'post')
            b_tot, se_tot = fe(tot, 'total')
            cov = joint.cov_params().loc[['pre', 'post'], ['pre', 'post']].to_numpy()
            se_sum = float(np.sqrt(cov.sum()))                     # se of b_pre + b_post
            sd_tot = float(d['total'].std())
            b_tot_sd, se_tot_sd = b_tot * sd_tot, se_tot * sd_tot  # days per +1 SD of total GPP
            p_tost = max(norm.sf((b_tot_sd + SESOI_DAYS_PER_SD) / se_tot_sd),
                         norm.cdf((b_tot_sd - SESOI_DAYS_PER_SD) / se_tot_sd))
            r = {**base, 'model': 'joint',
                 'lme_b_pre': b_pre, 'lme_se_pre': se_pre, 'lme_p_pre_negative': float(norm.cdf(b_pre / se_pre)),
                 'lme_b_post': b_post, 'lme_se_post': se_post, 'lme_p_post_positive': float(norm.sf(b_post / se_post)),
                 'lme_p_mirror_sum_zero': float(2 * norm.sf(abs((b_pre + b_post) / se_sum))),
                 'lme_b_total_days_per_sd': b_tot_sd, 'lme_se_total_days_per_sd': se_tot_sd,
                 'lme_p_total_two_sided': float(2 * norm.sf(abs(b_tot / se_tot))),
                 'lme_p_total_equivalent_to_zero': float(p_tost),
                 'lme_b_pre_days_per_sd': b_pre * d['pre'].std(), 'lme_b_post_days_per_sd': b_post * d['post'].std()}

            # headline: strictly within sites
            wj, wt = within_fit(ws, target, ['pre', 'post']), within_fit(ws, target, ['total'])
            if wj is None or wt is None:
                continue
            (bp, sp, _), (bq, sq, _), (bt, st, pt) = wj[0]['pre'], wj[0]['post'], wt[0]['total']
            wd = ws[['site_id', 'pre', 'post']].dropna()
            wdm = wd[['pre', 'post']] - wd.groupby('site_id')[['pre', 'post']].transform('mean')
            r.update({'n_within': wj[1], 'n_sites_within': wj[2],
                      'b_pre_days_per_sd': bp, 'se_pre_days_per_sd': sp, 'p_pre_negative': one_sided(bp, sp, True),
                      'b_post_days_per_sd': bq, 'se_post_days_per_sd': sq, 'p_post_positive': one_sided(bq, sq, False),
                      'b_total_days_per_sd': bt, 'se_total_days_per_sd': st, 'p_total_two_sided': pt,
                      'total_ci_lo': bt - 1.96 * st, 'total_ci_hi': bt + 1.96 * st,
                      'p_total_equivalent_to_zero': float(max(norm.sf((bt + SESOI_DAYS_PER_SD) / st),
                                                               norm.cdf((bt - SESOI_DAYS_PER_SD) / st))),
                      'r_pre_post_within': float(np.corrcoef(wdm['pre'], wdm['post'])[0, 1])})
            # with leaf-out date and air temperature controlled
            covars = [c for c in COVARS if ws[c].notna().any()]
            cj, ct = within_fit(ws, target, ['pre', 'post'] + covars), within_fit(ws, target, ['total'] + covars)
            if cj is not None and ct is not None:
                (cbp, csp, _), (cbq, csq, _), (cbt, cst, cpt) = cj[0]['pre'], cj[0]['post'], ct[0]['total']
                r.update({'ctrl_n': cj[1], 'ctrl_b_pre_days_per_sd': cbp, 'ctrl_se_pre_days_per_sd': csp,
                          'ctrl_p_pre_negative': one_sided(cbp, csp, True),
                          'ctrl_b_post_days_per_sd': cbq, 'ctrl_se_post_days_per_sd': csq,
                          'ctrl_p_post_positive': one_sided(cbq, csq, False),
                          'ctrl_b_total_days_per_sd': cbt, 'ctrl_se_total_days_per_sd': cst, 'ctrl_p_total_two_sided': cpt})
                for c in covars:
                    r[f'ctrl_b_{c}_days_per_sd'], _, r[f'ctrl_p_{c}'] = cj[0][c]
            if anchor == 'year':                    # slope produced by window length alone
                nl = within_fit(ws, target, ['pre_len', 'post_len'])
                if nl is not None:
                    r['b_post_length_only_days_per_sd'], r['se_post_length_only_days_per_sd'], _ = nl[0]['post_len']
            r['H1_pre_negative'] = r['p_pre_negative'] < ALPHA
            r['H2_post_positive'] = r['p_post_positive'] < ALPHA
            r['H3_total_no_effect'] = r['p_total_equivalent_to_zero'] < ALPHA
            r['hypothesis_supported'] = r['H1_pre_negative'] and r['H2_post_positive'] and r['H3_total_no_effect']
            rows.append(r)
res = pd.DataFrame(rows)
if res.empty:
    raise SystemExit("No EOS source x group had enough data.")
res.to_csv(OUT_CSV, index=False)
joint_res = res[res['model'] == 'joint']
print(f"\n{len(res)} tests ({len(joint_res)} with all three windows) -> '{OUT_CSV}'")

pd.set_option('display.width', 260)
show = ['group', 'eos_source', 'n_within', 'n_sites_within', 'b_pre_days_per_sd', 'p_pre_negative',
        'b_post_days_per_sd', 'p_post_positive', 'b_total_days_per_sd', 'total_ci_lo', 'total_ci_hi',
        'hypothesis_supported', 'ctrl_b_pre_days_per_sd', 'ctrl_p_pre_negative', 'ctrl_b_post_days_per_sd',
        'ctrl_p_post_positive', 'ctrl_b_total_days_per_sd']
show = [c for c in show if c in res.columns]
for (target, anchor), desc in anchors.items():
    a = res[(res['target'] == target) & (res['anchor'] == anchor)]
    if a.empty:
        continue
    print(f"\n=== {target} | {anchor}: {desc}")
    print(a[show + (['b_post_length_only_days_per_sd'] if anchor == 'year' and 'b_post_length_only_days_per_sd' in a else [])]
          .round(2).to_string(index=False))

# ---------------------------------------------------------------- figure 1: coefficients
colors = {'NDVI': '#2b6cb0', 'NIRv': '#805ad5', 'NDVI_tower': '#dd6b20', 'GCC': '#2f855a'}
for target in TARGETS:
    t_anchors = [a for (t, a) in anchors if t == target]
    terms = [('b_pre_days_per_sd', 'se_pre_days_per_sd', None, 'pre\n(SOS10 -> SOL)'),
             ('b_post_days_per_sd', 'se_post_days_per_sd', None, f'post\n(SOL -> {target})'),
             ('b_total_days_per_sd', 'se_total_days_per_sd', None, f'total\n(SOS10 -> {target})')]
    for group in GROUPS:
        sub = joint_res[(joint_res['target'] == target) & (joint_res['group'] == group)]
        if sub.empty:
            continue
        fig, axes = plt.subplots(1, len(t_anchors), figsize=(5 * len(t_anchors), 4.2), sharey=True, squeeze=False)
        for ax, anchor in zip(axes[0], t_anchors):
            sa = sub[sub['anchor'] == anchor]
            for si, src in enumerate(sources):
                r = sa[sa['eos_source'] == src]
                if r.empty:
                    continue
                r = r.iloc[0]
                x = np.arange(len(terms)) + (si - (len(sources) - 1) / 2) * 0.17
                y, err = [], []
                for bsd, se, _, _ in terms:
                    y.append(r[bsd])
                    err.append(1.96 * r[se])
                ax.errorbar(x, y, yerr=err, fmt='o', ms=5, lw=1.2, capsize=2, color=colors.get(src, 'gray'),
                            label=f"{src} (n={int(r['n_within'])})")
            ax.axhline(0, color='k', lw=0.7)
            ax.axhspan(-SESOI_DAYS_PER_SD, SESOI_DAYS_PER_SD, color='0.9', zorder=0)
            ax.set_xticks(range(len(terms)), [t[3] for t in terms])
            ax.set_title(anchors[(target, anchor)], fontsize=8)
            ax.legend(fontsize=7)
        axes[0][0].set_ylabel(f'{target} shift (days per +1 within-site SD of GPP)')
        fig.suptitle(f"Split-GPP cancellation, within sites, {target} - {group} sites "
                     f"(bars: 95% CI; grey band: +/- {SESOI_DAYS_PER_SD:g} d = 'no effect')", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"coefficients_{target}_{group.replace('/', '_')}.png", dpi=150)
        plt.close(fig)

# ---------------------------------------------------------------- figure 2: scatter of within-site anomalies
for target in TARGETS:
    t_anchors = [a for (t, a) in anchors if t == target]
    for src in sources:
        fig, axes = plt.subplots(len(t_anchors), 3, figsize=(12, 3.5 * len(t_anchors)), sharey=True, squeeze=False)
        for i, anchor in enumerate(t_anchors):
            w = tables[(target, anchor)]
            ws = w[w['vi_index'] == src]
            end = 'EOS90' if anchor == 'fixed_before_senescence' else target
            for j, (col, lab) in enumerate((('pre', 'pre: SOS10 -> SOL'), ('post', f'post: SOL -> {end}'),
                                            ('total', f'total: SOS10 -> {end}'))):
                ax = axes[i, j]
                ax.axhline(0, color='0.85', lw=0.8)
                ax.axvline(0, color='0.85', lw=0.8)
                d = ws[[target, 'site_id', col, 'season_type']].dropna(subset=[target, col])
                dm = d[[target, col]] - d.groupby('site_id')[[target, col]].transform('mean')
                if len(dm) > 2 and dm[col].std() > 0:
                    x, y = dm[col].to_numpy() * UNIT, dm[target].to_numpy()
                    dry = (d['season_type'] == 'dry-summer').to_numpy()
                    ax.scatter(x[~dry], y[~dry], s=10, alpha=0.5, color=colors.get(src, 'gray'), label='summer-green')
                    ax.scatter(x[dry], y[dry], s=14, alpha=0.7, color='#c05621', marker='^', label='dry-summer')
                    slope, icpt = np.polyfit(x, y, 1)
                    xs = np.linspace(x.min(), x.max(), 50)
                    ax.plot(xs, slope * xs + icpt, color='#c53030', lw=1.5)
                    per = (f"{slope:+.1f} d per gC m-2 d-1" if anchor == 'rate' else f"{slope * 100:+.1f} d per 100 gC m-2")
                    ax.set_title(f"{anchor} - {lab}\nr = {np.corrcoef(x, y)[0, 1]:.2f}, {per}, n = {len(x)}", fontsize=9)
                    if i == 0 and j == 0:
                        ax.legend(fontsize=7)
                ax.set_xlabel('GPP anomaly (gC m-2 d-1, year minus site mean)' if anchor == 'rate'
                              else 'GPP anomaly (gC m-2, year minus site mean)')
                if j == 0:
                    ax.set_ylabel(f'{target} anomaly (days)')
        fig.suptitle(f"{src}: {target} vs cumulative GPP, within-site anomalies (one row per window version)", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"scatter_{target}_{src}.png", dpi=150)
        plt.close(fig)

# ---------------------------------------------------------------- written summary
yes = {True: 'yes', False: 'no'}


def num(v, fmt='+.1f'):
    return format(v, fmt) if pd.notna(v) else '-'


lines = ["# Split-GPP cancellation hypothesis (Zani and Lu both right?) - test results", "",
         "Generated by `Script/41_split_gpp_cancellation_test.py`.", "",
         "- H1: cumulative GPP SOS10 -> solstice shifts EOS earlier (b_pre < 0)",
         "- H2: cumulative GPP solstice -> EOS shifts EOS later (b_post > 0)",
         f"- H3: total GPP SOS10 -> EOS has no effect (b_total significantly inside +/- {SESOI_DAYS_PER_SD:g} days per SD)", "",
         "All slopes are WITHIN-SITE: days of EOS shift per +1 within-site SD of cumulative GPP (year minus site "
         "mean, SE clustered by site). 'controlled' adds the year's leaf-out date and the air temperature before "
         "and after the solstice.", ""]
for (target, anchor), desc in anchors.items():
    a = joint_res[(joint_res['target'] == target) & (joint_res['anchor'] == anchor)]
    if a.empty:
        continue
    lines += [f"## {target} - {anchor}: {desc}", "",
              "| group | EOS source | n (sites) | pre | H1 | post | H2 | total [95% CI] | H3 | all three | controlled: pre / post / total |"
              + (" post from window length alone |" if anchor == 'year' else ""),
              "|---|---|---|---|---|---|---|---|---|---|---|" + ("---|" if anchor == 'year' else "")]
    for _, r in a.iterrows():
        ctrl = "-"
        if pd.notna(r.get('ctrl_b_pre_days_per_sd', np.nan)):
            ctrl = (f"{r['ctrl_b_pre_days_per_sd']:+.1f} (p={r['ctrl_p_pre_negative']:.3f}) / "
                    f"{r['ctrl_b_post_days_per_sd']:+.1f} (p={r['ctrl_p_post_positive']:.3f}) / "
                    f"{r['ctrl_b_total_days_per_sd']:+.1f}")
        line = (f"| {r['group']} | {r['eos_source']} | {int(r['n_within'])} ({int(r['n_sites_within'])}) | "
                f"{r['b_pre_days_per_sd']:+.1f} | {yes[bool(r['H1_pre_negative'])]} (p={r['p_pre_negative']:.3f}) | "
                f"{r['b_post_days_per_sd']:+.1f} | {yes[bool(r['H2_post_positive'])]} (p={r['p_post_positive']:.3f}) | "
                f"{r['b_total_days_per_sd']:+.1f} [{r['total_ci_lo']:+.1f}, {r['total_ci_hi']:+.1f}] | "
                f"{yes[bool(r['H3_total_no_effect'])]} | **{yes[bool(r['hypothesis_supported'])]}** | {ctrl} |")
        if anchor == 'year':
            line += f" {num(r.get('b_post_length_only_days_per_sd', np.nan))} |"
        lines.append(line)
    lines += ["", f"H1 holds in {int(a['H1_pre_negative'].sum())}, H2 in {int(a['H2_post_positive'].sum())}, "
              f"H3 in {int(a['H3_total_no_effect'].sum())}, all three in {int(a['hypothesis_supported'].sum())} of {len(a)} tests.", ""]
pre_only = res[res['model'] != 'joint']
if len(pre_only):
    lines += ["## Groups without a post-solstice window (EOS before the solstice): pre-solstice GPP only", "",
              "| target | anchor | group | EOS source | n | pre | H1 |", "|---|---|---|---|---|---|---|"]
    for _, r in pre_only.iterrows():
        lines.append(f"| {r['target']} | {r['anchor']} | {r['group']} | {r['eos_source']} | {int(r['n_within'])} | "
                     f"{r['b_pre_days_per_sd']:+.1f} | {yes[bool(r['H1_pre_negative'])]} (p={r['p_pre_negative']:.3f}) |")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print()
for (target, anchor) in anchors:
    a = joint_res[(joint_res['target'] == target) & (joint_res['anchor'] == anchor)]
    if len(a):
        print(f"[{target} | {anchor}] H1 in {int(a['H1_pre_negative'].sum())}/{len(a)}, H2 in {int(a['H2_post_positive'].sum())}/{len(a)}, "
              f"H3 in {int(a['H3_total_no_effect'].sum())}/{len(a)}, all three in {int(a['hypothesis_supported'].sum())}/{len(a)}")
print(f"\nSummary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")
