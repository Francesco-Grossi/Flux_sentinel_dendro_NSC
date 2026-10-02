"""
PIPELINE STEP 35 - How much do carbon source / sink activities explain EOS
beyond environmental drivers? (Notion D.1-3a, Eq.1-3, after Zohner et al.
2023 Fig. 5). Extends the legacy step 08/09 model comparison to EOS10/EOS50.

Predictor blocks (fixed-anchor windows from step 27):
  ENV     Eq.1  TA, SW, P before and after the solstice
                ([SOS,SOL] and [SOL,EOS10-site-mean])
  SOS     leaf-out (SOS10)   <- the open "should SOS be a predictor?" question
  SOURCE  Eq.2  cumulative GPP  [SOS,SOL] and [SOL,EOS10-site-mean]
  SINK    Eq.3  cumulative NPP  [SOS,SOL] and [SOL,EOS10-site-mean]
                (only if data/npp_luo2025_daily.csv exists)

NOTE: unlike legacy step 09, ENV deliberately excludes temperature/photoperiod
"senescence rate" - those (legacy step 07) are computed over EOS90->EOS10 of the
same year, i.e. from the target itself, which is circular for an EOS10 target.

Models (all fit on identical complete-case rows so AIC and the nested tests
are comparable):
  M0 ENV | M1 ENV+SOS | M2 ENV+SOURCE | M3 ENV+SINK | M4 ENV+SOURCE+SINK
  | M5 ENV+SOS+SOURCE(+SINK)      (M3/M4 only when NPP is available)
All models are fitted WITHIN SITES (year minus site mean, sites with >= 3
years, SE clustered by site): the question is what explains a late or early
year at a site, not why sites differ. The random-intercept mixed models used
before also absorbed between-site differences and overstated every effect.
For each: AIC/BIC, within-site R2 (share of the year-to-year variance
explained), its gain over ENV, max VIF, leave-one-site-out CV RMSE/R2
(predicting the anomalies of a site the model has not seen), nested Wald
tests (cluster-robust), and the coefficients of the fullest model (days per
+1 within-site SD).

Input : data/eos_window_predictors_fixed_anchor.csv   (step 27)
Output: data/eos_env_vs_carbon_comparison.csv
        data/eos_env_vs_carbon_lrt.csv
        data/eos_env_vs_carbon_cv.csv
        data/eos_env_vs_carbon_coefficients.csv
"""
import numpy as np
import pandas as pd
import eos_common as ec

OUT_CMP = ec.DATA_DIR / "eos_env_vs_carbon_comparison.csv"
OUT_LRT = ec.DATA_DIR / "eos_env_vs_carbon_lrt.csv"
OUT_CV = ec.DATA_DIR / "eos_env_vs_carbon_cv.csv"
OUT_COEF = ec.DATA_DIR / "eos_env_vs_carbon_coefficients.csv"
ec.require(ec.WINDOW_FIXED_CSV, hint="Run 27_window_predictors.py first.")

w = pd.read_csv(ec.WINDOW_FIXED_CSV)
PRE, POST = 'SOS_to_SOL', 'SOL_to_EOS10'


def existing(cols):
    return [c for c in cols if c in w.columns]


ENV = existing([f'{v}__{win}' for win in (PRE, POST)
                for v in ('TA_mean', 'SW_mean', 'P_sum')])
SOS = existing(['SOS'])
SOURCE = existing([f'GPP_cum__{PRE}', f'GPP_cum__{POST}'])
SINK = existing([f'NPP_cum__{PRE}', f'NPP_cum__{POST}'])       # NPP from the annual CUE
SINK_D = existing([f'NPPd_cum__{PRE}', f'NPPd_cum__{POST}'])   # NPP from the seasonal (daily) CUE
print(f"ENV={ENV}\nSOS={SOS}\nSOURCE={SOURCE}\nSINK={SINK or 'n/a (no NPP file)'}")
if not ENV or not SOURCE:
    raise SystemExit("Missing ENV or SOURCE columns - check step 27 / merged flux columns.")

BLOCKS = {'M0_env': ENV, 'M1_env+SOS': ENV + SOS, 'M2_env+source': ENV + SOURCE}
if SINK:
    BLOCKS['M3_env+sink'] = ENV + SINK
    BLOCKS['M4_env+source+sink'] = ENV + SOURCE + SINK
if SINK_D:
    BLOCKS['M3d_env+sink_daily'] = ENV + SINK_D
    BLOCKS['M4d_env+source+sink_daily'] = ENV + SOURCE + SINK_D
BLOCKS['M5_env+SOS+source' + ('+sink' if SINK else '')] = ENV + SOS + SOURCE + SINK
FULL = list(BLOCKS)[-1]

cmp_rows, lrt_rows, cv_rows, coef_rows = [], [], [], []
for vi in sorted(w['vi_index'].unique()):
    for target in ec.PRIMARY_TARGETS:
        all_vars = sorted(set(sum(BLOCKS.values(), [])))
        d = w[w['vi_index'] == vi][[target, 'site_id'] + all_vars].dropna()
        d = d[d.groupby('site_id')['site_id'].transform('size') >= ec.MIN_YEARS_WITHIN].reset_index(drop=True)
        n_obs, n_sites = len(d), d['site_id'].nunique()
        tag = f"{vi} / {target}"
        if n_obs < ec.MIN_OBS or n_sites < ec.MIN_SITES:
            print(f"\n[{tag}] skipped: {n_obs} obs / {n_sites} sites")
            continue
        print(f"\n{'=' * 70}\n{tag}: n_obs={n_obs}, n_sites={n_sites}\n{'=' * 70}")

        fits_ml, r2 = {}, {}
        for name, xs in BLOCKS.items():
            f = ec.fit_within(d, target, xs)
            if f is None:
                continue
            fits_ml[name] = f
            r2[name] = float(f.rsquared)
            frame = ec.within_frame(d, target, xs)
            cmp_rows.append({'vi_index': vi, 'target': target, 'model': name, 'n_predictors': len(xs),
                             'n_obs': n_obs, 'n_sites': n_sites, 'aic': f.aic, 'bic': f.bic,
                             'loglik': f.llf, 'r2_within': float(f.rsquared),
                             'max_vif': ec.max_vif(frame, [f'v{i}' for i in range(len(xs))])})
            cvr = ec.loso_within(d, target, xs)
            cv_rows.append({'vi_index': vi, 'target': target, 'model': name, **cvr})
            if name == FULL:
                ct = ec.within_coef_table(f, xs)
                ct.insert(0, 'model', name); ct.insert(0, 'target', target); ct.insert(0, 'vi_index', vi)
                coef_rows.append(ct)

        base = r2.get('M0_env', np.nan)
        for r in cmp_rows:
            if r['vi_index'] == vi and r['target'] == target:
                r['delta_r2_within_vs_env'] = r['r2_within'] - base
        cur = pd.DataFrame([r for r in cmp_rows if r['vi_index'] == vi and r['target'] == target])
        cvd = pd.DataFrame([r for r in cv_rows if r['vi_index'] == vi and r['target'] == target])
        print(cur.merge(cvd[['model', 'loso_rmse_days', 'loso_r2']], on='model')
              [['model', 'aic', 'r2_within', 'delta_r2_within_vs_env', 'max_vif', 'loso_rmse_days', 'loso_r2']]
              .round(3).to_string(index=False))

        for red, xr in BLOCKS.items():
            for full, xf in BLOCKS.items():
                if set(xr) < set(xf) and red in fits_ml and full in fits_ml:
                    extra = [f'v{i}' for i, x in enumerate(xf) if x not in xr]   # the added terms are jointly zero?
                    try:
                        wt = fits_ml[full].wald_test(', '.join(f'{v} = 0' for v in extra), use_f=False, scalar=True)
                        stat, p = float(wt.statistic), float(wt.pvalue)
                    except Exception:
                        stat, p = np.nan, np.nan
                    lrt_rows.append({'vi_index': vi, 'target': target, 'reduced': red, 'full': full,
                                     'df_diff': len(extra), 'chi2': stat, 'p_value': p,
                                     'test': 'Wald, cluster-robust'})

pd.DataFrame(cmp_rows).to_csv(OUT_CMP, index=False)
pd.DataFrame(lrt_rows).to_csv(OUT_LRT, index=False)
pd.DataFrame(cv_rows).to_csv(OUT_CV, index=False)
if coef_rows:
    pd.concat(coef_rows).to_csv(OUT_COEF, index=False)
print(f"\nWritten: {OUT_CMP.name}, {OUT_LRT.name}, {OUT_CV.name}, {OUT_COEF.name}")
