"""
PIPELINE STEP 44 - WHY does high pre-solstice GPP go with an earlier EOS90?
Three candidate explanations, each tested within sites (year minus site mean,
sites with >= 3 years, SE clustered by site; slopes in days per +1
within-site SD).

A. Leaf-out (deciduous forests). An early spring both raises pre-solstice
   GPP and could advance senescence by itself (a leaf has a limited life
   span). Leaf-out date (SOS10) and early GPP are put in the same model, for
   three versions of early GPP:
       sos    cumulative GPP, leaf-out -> solstice   (contains leaf-out date:
              an earlier leaf-out makes the window longer)
       cal    cumulative GPP, CAL_PRE days before the solstice
       rate   mean daily GPP, leaf-out -> solstice
   Reported: each alone, both together, how strongly the two are correlated
   within sites, and the share of the year-to-year EOS90 variance that only
   one of them explains (R2 of the joint model minus R2 of the other alone).

B. Water (dry-summer sites). Where the canopy dries out in early summer, a
   wet spring could raise GPP and move senescence at the same time. Climatic
   water balance (P - PET, PET after Makkink) before the solstice is added
   to the model, with pre-solstice air temperature. Two measures of water:
       wb_pre   climatic water balance, P - PET summed over the WB_PRE days
                before the solstice (every site)
       rew_pre  plant-available soil water: measured soil water content of
                the shallowest sensor (SWC_F_MDS_1, step 11), expressed as
                relative extractable water REW = (SWC - dry) / (wet - dry)
                with the site's own REW_LOW / REW_HIGH quantiles as dry and
                wet end; mean over the same WB_PRE days. Only sites that
                measure soil water.
       rew_post the same over the CAL_POST days from the solstice

C. Sink instead of source. If senescence responds to how much the plant
   could GROW rather than to how much it photosynthesised, sink variables
   should explain EOS90 better than GPP. Variables, daily:
       NPPd   daily-CUE x GPP                                   (step 26)
       Ra     autotrophic respiration = (1 - CUE) x GPP
       Rg     growth respiration      = gR x CUE x GPP   (gR: the year's
              growth-respiration coefficient from step 26)
       Rm     maintenance respiration = Ra - Rg
       Rm_model  maintenance respiration from the fitted model of step 26:
              temperature response x biomass built since 1 January. Unlike
              the three above it does not contain the day's GPP.
       SI     sink-limitation index = fT x fW, independent of the fluxes:
              fT = 0 below SINK_T_MIN, rising linearly to 1 at SINK_T_OPT
                   (cambial growth stops near 5 C);
              fW = relative filling of a soil-water bucket of BUCKET_MM mm
                   driven by P and Makkink PET (growth needs turgor).
              A proxy for when tissue growth is possible, not a measurement.
   NPPd, Ra, Rg and Rm are GPP multiplied by a factor, so they are strongly
   correlated with GPP by construction (r reported); SI is not.
   Each variable replaces GPP in the calendar-window model of step 41:
       dEOS90 ~ dpre + dpost     (CAL_PRE days before / CAL_POST days from the solstice)

D. Spring temperature instead of GPP. A warm spring raises GPP and could
   advance senescence by itself (faster development). Air temperature and
   GPP of the CAL_PRE days before the solstice are put in the same model,
   then with radiation and water balance added: the GPP slope of that last
   model is the effect of the GPP that the weather does not explain.
   The direct test uses the years in which the two DIVERGE. Each year is
   classed by the sign of its within-site temperature and GPP anomaly:
       warm + high GPP, cool + low GPP    both explanations agree
       warm + low GPP                     early EOS90 if temperature drives it
       cool + high GPP                    early EOS90 if GPP drives it
   and the mean EOS90 anomaly of each class is reported.

For A, B and D the model is also fitted on all EOS sources stacked (anomalies
relative to each site x source mean, SE clustered by site), which uses every
site-year once per source: more power, same question.

Input : phenology tables (steps 23, 24, 25), daily flux, data/cue_luo2025_site_year.csv,
        data/site_season_type.csv (step 41)
Output: data/mechanism_leafout_vs_gpp.csv
        data/mechanism_water.csv
        data/mechanism_sink.csv
        data/mechanism_temperature.csv
        data/mechanism_temperature_quadrants.csv
        figure/mechanism/leafout_vs_gpp.png, water.png, sink.png, temperature.png
        Output/mechanism_summary.md
"""
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

OUT_A = ec.DATA_DIR / "mechanism_leafout_vs_gpp.csv"
OUT_B = ec.DATA_DIR / "mechanism_water.csv"
OUT_C = ec.DATA_DIR / "mechanism_sink.csv"
OUT_D = ec.DATA_DIR / "mechanism_temperature.csv"
OUT_DQ = ec.DATA_DIR / "mechanism_temperature_quadrants.csv"
OUT_MD = ec.ROOT / "Output" / "mechanism_summary.md"
SEASON_CSV = ec.DATA_DIR / "site_season_type.csv"
CUE_CSV = ec.DATA_DIR / "cue_luo2025_site_year.csv"
FIG_DIR = ec.FIGURE_DIR / "mechanism"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_MD.parent.mkdir(parents=True, exist_ok=True)

TARGET = 'EOS90'
CAL_PRE, CAL_POST, WB_PRE = 60, 45, 90
REW_LOW, REW_HIGH = 0.02, 0.98          # site quantiles of soil water content taken as dry and wet end
SINK_T_MIN, SINK_T_OPT, BUCKET_MM = 5.0, 18.0, 150.0
MIN_YEARS, MIN_OBS, MIN_SITES = 3, 25, 5
LEAF_HABIT = {'DBF': 'deciduous', 'DNF': 'deciduous', 'MF': 'deciduous', 'ENF': 'evergreen', 'EBF': 'evergreen'}
SRC_ORDER = ['GCC', 'NDVI_tower', 'NDVI', 'NIRv']
NAME = {'GCC': 'PhenoCam', 'NDVI_tower': 'tower NDVI', 'NDVI': 'satellite NDVI', 'NIRv': 'satellite NIRv',
        'stacked': 'all sources stacked'}
COLORS = {'NDVI': '#2b6cb0', 'NIRv': '#805ad5', 'NDVI_tower': '#dd6b20', 'GCC': '#2f855a', 'stacked': 'k'}

# ---------------------------------------------------------------- daily variables
flux = ec.load_flux_daily().sort_values(['site_id', 'year', 'doy']).reset_index(drop=True)
flux['PET'] = ec.makkink_pet(flux['SW'], flux['TA'])
flux['WB'] = flux['P'] - flux['PET']

sink_vars = ['GPP']
if 'CUEd' in flux.columns:
    cue = flux['CUEd'].clip(0, 1)
    flux['Ra'] = (1 - cue) * flux['GPP']
    sink_vars += ['NPPd', 'Ra']
    if os.path.exists(CUE_CSV):
        gr = pd.read_csv(CUE_CSV)[['site_id', 'year', 'gR']]
        flux = flux.merge(gr, on=['site_id', 'year'], how='left')
        flux['Rg'] = (flux['gR'] * cue * flux['GPP']).clip(upper=flux['Ra'])
        flux['Rm'] = flux['Ra'] - flux['Rg']
        sink_vars += ['Rg', 'Rm']
if 'Rm_model' in flux.columns and flux['Rm_model'].notna().any():
    sink_vars += ['Rm_model']


def bucket(g):
    """Relative soil-water filling (0-1) of a simple bucket, run continuously through a site's record."""
    p, pet = g['P'].fillna(0).to_numpy(float), g['PET'].fillna(0).to_numpy(float)
    w, out = BUCKET_MM, np.empty(len(g))
    for i in range(len(g)):
        w = min(BUCKET_MM, max(0.0, w + p[i] - pet[i] * w / BUCKET_MM))
        out[i] = w / BUCKET_MM
    return pd.Series(out, index=g.index)


flux['fW'] = pd.concat([bucket(g) for _, g in flux.groupby('site_id', sort=False)]).where(flux['P'].notna() & flux['PET'].notna())
flux['fT'] = ((flux['TA'] - SINK_T_MIN) / (SINK_T_OPT - SINK_T_MIN)).clip(0, 1)
flux['SI'] = flux['fT'] * flux['fW']
sink_vars += ['SI', 'fT', 'fW']
if 'SWC' in flux.columns and flux['SWC'].notna().any():
    # plant-available water: soil water content relative to the site's own dry and wet ends
    q = flux.groupby('site_id')['SWC'].quantile([REW_LOW, REW_HIGH]).unstack()
    lo, hi = flux['site_id'].map(q[REW_LOW]), flux['site_id'].map(q[REW_HIGH])
    flux['REW'] = ((flux['SWC'] - lo) / (hi - lo).where(hi > lo)).clip(0, 1)
lk = ec.FluxLookup(flux, sink_vars + ['WB', 'TA', 'SW', 'REW'])

pheno = ec.load_phenology()
habit = flux.drop_duplicates('site_id').set_index('site_id')['igbp'].map(LEAF_HABIT).fillna('grass/shrub')
season = pd.read_csv(SEASON_CSV).set_index('site_id')['season_type'] if os.path.exists(SEASON_CSV) else pd.Series(dtype=str)
pheno['leaf_habit'] = pheno['site_id'].map(habit)
pheno['season_type'] = pheno['site_id'].map(season)
sources = [s for s in SRC_ORDER if s in set(pheno['vi_index'])]

sol = pheno['year'].map(ec.solstice_doy).astype(float)
rows = list(zip(pheno['site_id'], pheno['year'], sol))


def win(col, a, b, total=False):
    """Window mean (or sum) per phenology row; a, b are offsets from the solstice."""
    n = b - a + 1
    return [lk.window_mean(s, y, col, so + a, so + b) * (n if total else 1) for s, y, so in rows]


pheno['SOS'] = pheno['leaf_out_10']
pheno['gpp_cal'] = win('GPP', -CAL_PRE, -1, total=True)
pheno['gpp_rate'] = [lk.window_mean(s, y, 'GPP', a, so) for (s, y, so), a in zip(rows, pheno['SOS'])]
pheno['gpp_sos'] = pheno['gpp_rate'] * [ec.FluxLookup.n_days(a, so) for (_, _, so), a in zip(rows, pheno['SOS'])]
pheno['wb_pre'] = win('WB', -WB_PRE, -1, total=True)
pheno['wb_post'] = win('WB', 0, CAL_POST - 1, total=True)
WATER_VARS = {'wb_pre': 'water balance (P - PET)'}
if 'REW' in flux.columns:
    pheno['rew_pre'] = win('REW', -WB_PRE, -1)
    pheno['rew_post'] = win('REW', 0, CAL_POST - 1)
    WATER_VARS['rew_pre'] = 'plant-available soil water'
    WATER_VARS['rew_post'] = 'plant-available soil water, after the solstice'
    print(f"Soil water: {int(pheno.drop_duplicates(['site_id', 'year'])['rew_pre'].notna().sum())} site-years at "
          f"{pheno.loc[pheno['rew_pre'].notna(), 'site_id'].nunique()} sites have measured soil water before the solstice.")
pheno['ta_pre'] = win('TA', -CAL_PRE, -1)
pheno['sw_pre'] = win('SW', -CAL_PRE, -1)
pheno['wb_cal'] = win('WB', -CAL_PRE, -1, total=True)
for v in sink_vars:
    pheno[f'{v}_pre'] = win(v, -CAL_PRE, -1)
    pheno[f'{v}_post'] = win(v, 0, CAL_POST - 1)


def wfit(d, xs, unit='site_id'):
    """Within-`unit` OLS of TARGET on xs, SE clustered by site. Returns dict of (b, se, p), n, sites, R2."""
    d = d[['site_id', unit, TARGET] + xs].loc[:, lambda f: ~f.columns.duplicated()].dropna()
    d = d[d.groupby(unit)[unit].transform('size') >= MIN_YEARS]
    if len(d) < MIN_OBS or d['site_id'].nunique() < MIN_SITES:
        return None
    dm = d[[TARGET] + xs] - d.groupby(unit)[[TARGET] + xs].transform('mean')
    sd = dm[xs].std()
    if (sd == 0).any() or sd.isna().any():
        return None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        f = sm.OLS(dm[TARGET].to_numpy(float), sm.add_constant((dm[xs] / sd).to_numpy(float))).fit(
            cov_type='cluster', cov_kwds={'groups': pd.factorize(d['site_id'])[0]})
    z = dm[xs] / sd
    return {'coef': {x: (float(f.params[i + 1]), float(f.bse[i + 1]), float(f.pvalues[i + 1])) for i, x in enumerate(xs)},
            'n_obs': len(d), 'n_sites': d['site_id'].nunique(), 'r2': float(f.rsquared),
            'corr': z.corr() if len(xs) > 1 else None}


pheno['unit'] = pheno['site_id'] + '|' + pheno['vi_index']


def samples(group_col, group):
    """(label, frame, unit column) for each EOS source and for all sources stacked."""
    g = pheno if group == 'all' else pheno[pheno[group_col] == group]
    out = [(src, g[g['vi_index'] == src], 'site_id') for src in sources]
    out.append(('stacked', g, 'unit'))
    return out


def put(r, prefix, fit, x):
    b, se, p = fit['coef'][x]
    r.update({f'{prefix}_b': b, f'{prefix}_se': se, f'{prefix}_p': p})


# ---------------------------------------------------------------- A. leaf-out vs early GPP
a_rows = []
for group in ('deciduous', 'evergreen', 'grass/shrub', 'all'):
    for src, d, unit in samples('leaf_habit', group):
        for gv in ('gpp_sos', 'gpp_cal', 'gpp_rate'):
            dd = d[['site_id', unit, TARGET, 'SOS', gv]].loc[:, lambda f: ~f.columns.duplicated()].dropna()   # identical rows for the three models
            fs, fg, fj = wfit(dd, ['SOS'], unit), wfit(dd, [gv], unit), wfit(dd, [gv, 'SOS'], unit)
            if fs is None or fg is None or fj is None:
                continue
            r = {'group': group, 'eos_source': src, 'gpp_version': gv.replace('gpp_', ''), 'n_obs': fj['n_obs'],
                 'n_sites': fj['n_sites'], 'r_gpp_leafout_within': float(fj['corr'].iloc[0, 1]),
                 'r2_leafout': fs['r2'], 'r2_gpp': fg['r2'], 'r2_both': fj['r2'],
                 'r2_unique_gpp': fj['r2'] - fs['r2'], 'r2_unique_leafout': fj['r2'] - fg['r2']}
            put(r, 'leafout_alone', fs, 'SOS')
            put(r, 'gpp_alone', fg, gv)
            put(r, 'gpp_joint', fj, gv)
            put(r, 'leafout_joint', fj, 'SOS')
            a_rows.append(r)
A = pd.DataFrame(a_rows)
A.to_csv(OUT_A, index=False)
pd.set_option('display.width', 250)
print(f"A. Leaf-out date vs early GPP -> '{OUT_A}'")
show = ['group', 'eos_source', 'gpp_version', 'n_obs', 'n_sites', 'r_gpp_leafout_within', 'gpp_alone_b', 'gpp_alone_p',
        'leafout_alone_b', 'leafout_alone_p', 'gpp_joint_b', 'gpp_joint_p', 'leafout_joint_b', 'leafout_joint_p',
        'r2_unique_gpp', 'r2_unique_leafout']
print(A[A['group'].isin(['deciduous', 'all'])][show].round(3).to_string(index=False))

# ---------------------------------------------------------------- B. water
b_rows = []
for group in ('dry-summer', 'summer-green', 'all'):
    for src, d, unit in samples('season_type', group):
        for wv in WATER_VARS:
            for gv in ('gpp_sos', 'gpp_cal'):
                dd = d[['site_id', unit, TARGET, gv, wv, 'ta_pre']].loc[:, lambda f: ~f.columns.duplicated()].dropna()
                fg, fw = wfit(dd, [gv], unit), wfit(dd, [wv], unit)
                fj, ft = wfit(dd, [gv, wv], unit), wfit(dd, [gv, wv, 'ta_pre'], unit)
                if fg is None or fw is None or fj is None or ft is None:
                    continue
                r = {'group': group, 'eos_source': src, 'water': WATER_VARS[wv], 'gpp_version': gv.replace('gpp_', ''),
                     'n_obs': fj['n_obs'], 'n_sites': fj['n_sites'], 'r_gpp_water_within': float(fj['corr'].iloc[0, 1]),
                     'r2_gpp': fg['r2'], 'r2_water': fw['r2'], 'r2_both': fj['r2'],
                     'r2_unique_gpp': fj['r2'] - fw['r2'], 'r2_unique_water': fj['r2'] - fg['r2']}
                put(r, 'gpp_alone', fg, gv)
                put(r, 'water_alone', fw, wv)
                put(r, 'gpp_with_water', fj, gv)
                put(r, 'water_with_gpp', fj, wv)
                put(r, 'gpp_with_water_T', ft, gv)
                put(r, 'water_with_gpp_T', ft, wv)
                put(r, 'T_with_gpp_water', ft, 'ta_pre')
                b_rows.append(r)
B = pd.DataFrame(b_rows)
B.to_csv(OUT_B, index=False)
print(f"\nB. Water vs pre-solstice GPP -> '{OUT_B}'")
show = ['group', 'eos_source', 'water', 'gpp_version', 'n_obs', 'n_sites', 'r_gpp_water_within', 'gpp_alone_b', 'gpp_alone_p',
        'water_alone_b', 'water_alone_p', 'gpp_with_water_b', 'gpp_with_water_p', 'water_with_gpp_b', 'water_with_gpp_p',
        'gpp_with_water_T_b', 'gpp_with_water_T_p']
print(B[show].round(3).to_string(index=False))

# ---------------------------------------------------------------- C. sink variables in place of GPP
c_rows = []
for group in ('all', 'deciduous', 'evergreen'):
    for src, d, unit in samples('leaf_habit', group):
        if src == 'stacked':
            continue
        for v in sink_vars:
            fj = wfit(d, [f'{v}_pre', f'{v}_post'], unit)
            if fj is None:
                continue
            r = {'group': group, 'eos_source': src, 'variable': v, 'n_obs': fj['n_obs'], 'n_sites': fj['n_sites'], 'r2': fj['r2']}
            put(r, 'pre', fj, f'{v}_pre')
            put(r, 'post', fj, f'{v}_post')
            if v != 'GPP':                       # same rows: GPP alone, and the variable next to GPP
                dd = d[['site_id', TARGET, 'GPP_pre', 'GPP_post', f'{v}_pre', f'{v}_post']].dropna()
                dm = dd[['GPP_pre', f'{v}_pre']] - dd.groupby('site_id')[['GPP_pre', f'{v}_pre']].transform('mean')
                r['r_with_gpp_pre_within'] = float(dm.corr().iloc[0, 1]) if len(dm) > 2 else np.nan
                fg, fb = wfit(dd, ['GPP_pre', 'GPP_post']), wfit(dd, ['GPP_pre', f'{v}_pre'])
                if fg is not None and fb is not None:
                    r['r2_gpp_same_rows'] = fg['r2']
                    put(r, 'pre_next_to_gpp', fb, f'{v}_pre')
                    put(r, 'gpp_pre_next_to_it', fb, 'GPP_pre')
            c_rows.append(r)
C = pd.DataFrame(c_rows)
C.to_csv(OUT_C, index=False)
print(f"\nC. Sink variables in place of GPP (calendar windows) -> '{OUT_C}'")
show = [c for c in ['eos_source', 'variable', 'n_obs', 'n_sites', 'pre_b', 'pre_p', 'post_b', 'post_p', 'r2', 'r2_gpp_same_rows',
                    'r_with_gpp_pre_within', 'pre_next_to_gpp_b', 'pre_next_to_gpp_p', 'gpp_pre_next_to_it_b',
                    'gpp_pre_next_to_it_p'] if c in C.columns]
print(C[C['group'] == 'all'][show].round(3).to_string(index=False))

# ---------------------------------------------------------------- D. spring temperature vs GPP
QUADRANTS = [('warm + high GPP', 1, 1, 'both say earlier'), ('cool + low GPP', -1, -1, 'both say later'),
             ('warm + low GPP', 1, -1, 'earlier if temperature drives it'),
             ('cool + high GPP', -1, 1, 'earlier if GPP drives it')]
d_rows, q_rows = [], []
for group in ('all', 'summer-green', 'deciduous', 'evergreen'):
    col = 'season_type' if group == 'summer-green' else 'leaf_habit'
    for src, d, unit in samples(col, group):
        keep = list(dict.fromkeys(['site_id', unit, TARGET, 'gpp_cal', 'ta_pre', 'sw_pre', 'wb_cal']))
        dd = d[keep].dropna()
        fg, ft = wfit(dd, ['gpp_cal'], unit), wfit(dd, ['ta_pre'], unit)
        fj, fa = wfit(dd, ['gpp_cal', 'ta_pre'], unit), wfit(dd, ['gpp_cal', 'ta_pre', 'sw_pre', 'wb_cal'], unit)
        if fg is None or ft is None or fj is None or fa is None:
            continue
        r = {'group': group, 'eos_source': src, 'n_obs': fj['n_obs'], 'n_sites': fj['n_sites'],
             'r_gpp_temperature_within': float(fj['corr'].iloc[0, 1]), 'r2_gpp': fg['r2'], 'r2_temperature': ft['r2'],
             'r2_both': fj['r2'], 'r2_unique_gpp': fj['r2'] - ft['r2'], 'r2_unique_temperature': fj['r2'] - fg['r2']}
        put(r, 'gpp_alone', fg, 'gpp_cal')
        put(r, 'T_alone', ft, 'ta_pre')
        put(r, 'gpp_with_T', fj, 'gpp_cal')
        put(r, 'T_with_gpp', fj, 'ta_pre')
        put(r, 'gpp_with_weather', fa, 'gpp_cal')
        put(r, 'T_with_all', fa, 'ta_pre')
        put(r, 'SW_with_all', fa, 'sw_pre')
        put(r, 'water_with_all', fa, 'wb_cal')
        d_rows.append(r)

        # years in which temperature and GPP diverge
        dq = dd[dd.groupby(unit)[unit].transform('size') >= MIN_YEARS]
        an = dq[[TARGET, 'gpp_cal', 'ta_pre']] - dq.groupby(unit)[[TARGET, 'gpp_cal', 'ta_pre']].transform('mean')
        an['site_id'] = dq['site_id'].to_numpy()
        for name, st, sg, reading in QUADRANTS:
            q = an[(np.sign(an['ta_pre']) == st) & (np.sign(an['gpp_cal']) == sg)]
            if len(q) < 10 or q['site_id'].nunique() < MIN_SITES:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                f = sm.OLS(q[TARGET].to_numpy(float), np.ones(len(q))).fit(
                    cov_type='cluster', cov_kwds={'groups': pd.factorize(q['site_id'])[0]})
            q_rows.append({'group': group, 'eos_source': src, 'years': name, 'reading': reading, 'n_obs': len(q),
                           'n_sites': q['site_id'].nunique(), 'share_of_years': len(q) / len(an),
                           'mean_eos_anomaly_days': float(f.params[0]), 'se': float(f.bse[0]), 'p': float(f.pvalues[0])})
Dt, Dq = pd.DataFrame(d_rows), pd.DataFrame(q_rows)
Dt.to_csv(OUT_D, index=False)
Dq.to_csv(OUT_DQ, index=False)
print(f"\nD. Spring temperature vs GPP ({CAL_PRE} days before the solstice) -> '{OUT_D}'")
show = ['group', 'eos_source', 'n_obs', 'n_sites', 'r_gpp_temperature_within', 'gpp_alone_b', 'gpp_alone_p', 'T_alone_b',
        'T_alone_p', 'gpp_with_T_b', 'gpp_with_T_p', 'T_with_gpp_b', 'T_with_gpp_p', 'gpp_with_weather_b',
        'gpp_with_weather_p', 'r2_unique_gpp', 'r2_unique_temperature']
print(Dt[show].round(3).to_string(index=False))
print(f"\n   Years in which they diverge -> '{OUT_DQ}'")
print(Dq[Dq['group'] == 'all'].round(2).to_string(index=False))


# ---------------------------------------------------------------- figures
def forest(ax, d, specs, title):
    """specs: list of (prefix, label). One row of points per spec, one colour per EOS source."""
    srcs = [s for s in sources + ['stacked'] if s in set(d['eos_source'])]
    for si, src in enumerate(srcs):
        r = d[d['eos_source'] == src]
        if r.empty:
            continue
        r = r.iloc[0]
        y = np.arange(len(specs)) + (si - (len(srcs) - 1) / 2) * 0.15
        ax.errorbar([r.get(f'{p}_b', np.nan) for p, _ in specs], y,
                    xerr=[1.96 * r.get(f'{p}_se', np.nan) for p, _ in specs], fmt='o', ms=4, lw=1, capsize=2,
                    color=COLORS.get(src, 'gray'), label=f"{NAME[src]} (n={int(r['n_obs'])})")
    ax.axvline(0, color='k', lw=0.7)
    ax.set_yticks(range(len(specs)), [lab for _, lab in specs], fontsize=8)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6.5)
    ax.set_xlabel(f'{TARGET} shift (days per +1 within-site SD)')


if len(A):
    versions = [('sos', 'GPP leaf-out -> solstice'), ('cal', f'GPP {CAL_PRE} d before solstice'), ('rate', 'GPP rate leaf-out -> solstice')]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4), sharex=True)
    for ax, (gv, lab) in zip(axes, versions):
        forest(ax, A[(A['group'] == 'deciduous') & (A['gpp_version'] == gv)],
               [('gpp_alone', 'early GPP alone'), ('leafout_alone', 'leaf-out date alone'),
                ('gpp_joint', 'early GPP, leaf-out held fixed'), ('leafout_joint', 'leaf-out date, GPP held fixed')], lab)
    fig.suptitle("Deciduous forests: leaf-out date against early-season GPP (bars: 95% CI)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "leafout_vs_gpp.png", dpi=150)
    plt.close(fig)
if len(B):
    waters = [w for w in ('wb_pre', 'rew_pre') if WATER_VARS.get(w) in set(B['water'])]
    fig, axes = plt.subplots(len(waters), 2, figsize=(12, 4.4 * len(waters)), sharex=True, squeeze=False)
    for row, wv in zip(axes, waters):
        for ax, group in zip(row, ('dry-summer', 'summer-green')):
            forest(ax, B[(B['group'] == group) & (B['gpp_version'] == 'sos') & (B['water'] == WATER_VARS[wv])],
                   [('gpp_alone', 'pre-solstice GPP alone'), ('water_alone', 'water alone'),
                    ('gpp_with_water', 'GPP, water held fixed'), ('water_with_gpp', 'water, GPP held fixed'),
                    ('gpp_with_water_T', 'GPP, water and temperature held fixed')],
                   f'{group} sites - {WATER_VARS[wv]}')
    fig.suptitle(f"Pre-solstice GPP (leaf-out -> solstice) against water in the {WB_PRE} days before the solstice", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "water.png", dpi=150)
    plt.close(fig)
if len(C):
    ca = C[C['group'] == 'all']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, (pre, lab) in zip(axes, (('pre', f'{CAL_PRE} days before the solstice'), ('post', f'{CAL_POST} days from the solstice'))):
        for si, src in enumerate(sources):
            s = ca[ca['eos_source'] == src].set_index('variable').reindex(sink_vars)
            y = np.arange(len(sink_vars)) + (si - (len(sources) - 1) / 2) * 0.17
            ax.errorbar(s[f'{pre}_b'], y, xerr=1.96 * s[f'{pre}_se'], fmt='o', ms=4, lw=1, capsize=2,
                        color=COLORS.get(src, 'gray'), label=NAME[src])
        ax.axvline(0, color='k', lw=0.7)
        ax.set_yticks(range(len(sink_vars)), sink_vars)
        ax.set_title(lab, fontsize=9)
        ax.set_xlabel(f'{TARGET} shift (days per +1 within-site SD)')
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=7)
    fig.suptitle("Source (GPP) and sink variables in the same two-window model, all sites (bars: 95% CI)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "sink.png", dpi=150)
    plt.close(fig)
if len(Dt) and len(Dq):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    forest(axes[0], Dt[Dt['group'] == 'all'],
           [('gpp_alone', 'GPP alone'), ('T_alone', 'temperature alone'), ('gpp_with_T', 'GPP, temperature held fixed'),
            ('T_with_gpp', 'temperature, GPP held fixed'), ('gpp_with_weather', 'GPP, all weather held fixed')],
           f'all sites: slopes ({CAL_PRE} days before the solstice)')
    ax, qa = axes[1], Dq[Dq['group'] == 'all']
    srcs = [s for s in sources + ['stacked'] if s in set(qa['eos_source'])]
    names = [q[0] for q in QUADRANTS]
    for si, src in enumerate(srcs):
        s = qa[qa['eos_source'] == src].set_index('years').reindex(names)
        y = np.arange(len(names)) + (si - (len(srcs) - 1) / 2) * 0.15
        ax.errorbar(s['mean_eos_anomaly_days'], y, xerr=1.96 * s['se'], fmt='o', ms=4, lw=1, capsize=2,
                    color=COLORS.get(src, 'gray'), label=NAME[src])
    ax.axvline(0, color='k', lw=0.7)
    ax.set_yticks(range(len(names)), [f"{n}\n({r})" for n, _, _, r in QUADRANTS], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f'mean {TARGET} anomaly (days; negative = earlier)')
    ax.set_title('all sites: years classed by spring temperature and GPP anomaly', fontsize=9)
    ax.legend(fontsize=6.5)
    fig.suptitle("Spring temperature against pre-solstice GPP (bars: 95% CI)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "temperature.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- written summary
def c(r, prefix):
    b, p = r.get(f'{prefix}_b', np.nan), r.get(f'{prefix}_p', np.nan)
    return f"{b:+.1f}{ec.stars(p)}" if pd.notna(b) else '-'


lines = ["# Mechanism tests - why does high pre-solstice GPP go with an earlier EOS90?", "",
         "Generated by `Script/44_mechanism_tests.py`. Within-site models; days of EOS90 shift per +1 within-site SD. "
         "Stars: two-sided * p<0.05, ** p<0.01, *** p<0.001.", "",
         "## A. Leaf-out date against early GPP", "",
         "| sites | EOS source | early GPP | n (sites) | r(GPP, leaf-out) | GPP alone | leaf-out alone | GPP, leaf-out fixed | leaf-out, GPP fixed | R2 only GPP | R2 only leaf-out |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
for _, r in A.iterrows():
    lines.append(f"| {r['group']} | {NAME[r['eos_source']]} | {r['gpp_version']} | {int(r['n_obs'])} ({int(r['n_sites'])}) | "
                 f"{r['r_gpp_leafout_within']:+.2f} | {c(r, 'gpp_alone')} | {c(r, 'leafout_alone')} | {c(r, 'gpp_joint')} | "
                 f"{c(r, 'leafout_joint')} | {r['r2_unique_gpp']:.3f} | {r['r2_unique_leafout']:.3f} |")
lines += ["", f"## B. Water ({WB_PRE} days before the solstice) against pre-solstice GPP", "",
          "Plant-available soil water = measured soil water content relative to the site's own dry and wet ends.", "",
          "| sites | EOS source | water measure | GPP | n (sites) | r(GPP, water) | GPP alone | water alone | GPP, water fixed | water, GPP fixed | GPP, water + T fixed | R2 only GPP | R2 only water |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for _, r in B.iterrows():
    lines.append(f"| {r['group']} | {NAME[r['eos_source']]} | {r['water']} | {r['gpp_version']} | {int(r['n_obs'])} ({int(r['n_sites'])}) | "
                 f"{r['r_gpp_water_within']:+.2f} | {c(r, 'gpp_alone')} | {c(r, 'water_alone')} | {c(r, 'gpp_with_water')} | "
                 f"{c(r, 'water_with_gpp')} | {c(r, 'gpp_with_water_T')} | {r['r2_unique_gpp']:.3f} | {r['r2_unique_water']:.3f} |")
lines += ["", f"## C. Sink variables in place of GPP ({CAL_PRE} days before / {CAL_POST} days from the solstice)", "",
          "| sites | EOS source | variable | n (sites) | pre | post | R2 | R2 of GPP, same rows | r with GPP (pre) | variable next to GPP | GPP next to it |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
for _, r in C.iterrows():
    lines.append(f"| {r['group']} | {NAME[r['eos_source']]} | {r['variable']} | {int(r['n_obs'])} ({int(r['n_sites'])}) | {c(r, 'pre')} | "
                 f"{c(r, 'post')} | {r['r2']:.3f} | " + (f"{r['r2_gpp_same_rows']:.3f}" if pd.notna(r.get('r2_gpp_same_rows', np.nan)) else '-')
                 + " | " + (f"{r['r_with_gpp_pre_within']:+.2f}" if pd.notna(r.get('r_with_gpp_pre_within', np.nan)) else '-')
                 + f" | {c(r, 'pre_next_to_gpp')} | {c(r, 'gpp_pre_next_to_it')} |")
lines += ["", f"## D. Spring temperature against GPP ({CAL_PRE} days before the solstice)", "",
          "'all weather' = temperature, radiation and water balance of the same window.", "",
          "| sites | EOS source | n (sites) | r(GPP, T) | GPP alone | T alone | GPP, T fixed | T, GPP fixed | GPP, all weather fixed | T, all else fixed | R2 only GPP | R2 only T |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for _, r in Dt.iterrows():
    lines.append(f"| {r['group']} | {NAME[r['eos_source']]} | {int(r['n_obs'])} ({int(r['n_sites'])}) | "
                 f"{r['r_gpp_temperature_within']:+.2f} | {c(r, 'gpp_alone')} | {c(r, 'T_alone')} | {c(r, 'gpp_with_T')} | "
                 f"{c(r, 'T_with_gpp')} | {c(r, 'gpp_with_weather')} | {c(r, 'T_with_all')} | {r['r2_unique_gpp']:.3f} | "
                 f"{r['r2_unique_temperature']:.3f} |")
lines += ["", "**Mean EOS90 anomaly (days) by class of year** - the last two classes are the test", "",
          "| sites | EOS source | years | reading | n (sites) | share of years | mean EOS90 anomaly [95% CI] |",
          "|---|---|---|---|---|---|---|"]
for _, r in Dq.iterrows():
    m, se = r['mean_eos_anomaly_days'], r['se']
    lines.append(f"| {r['group']} | {NAME[r['eos_source']]} | {r['years']} | {r['reading']} | {int(r['n_obs'])} ({int(r['n_sites'])}) | "
                 f"{100 * r['share_of_years']:.0f}% | {m:+.1f}{ec.stars(r['p'])} [{m - 1.96 * se:+.1f}, {m + 1.96 * se:+.1f}] |")
OUT_MD.write_text("\n".join(lines) + "\n", encoding='utf-8')
print(f"\nSummary -> '{OUT_MD}'\nFigures -> '{FIG_DIR}'")
