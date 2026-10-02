"""
PIPELINE STEP 37 - Do growth RATE and accumulated carbon act the same way on
EOS? (Notion D.2: "GPPmean/NPPmean more directly link with EOS10 compared to
GPPcum/NPPcum?", "growth rate matters more", "effect differs among tree species".)

For a window, cumulative = rate x duration (cum = mean x len). So:
    R   : EOS ~ rate
    C   : EOS ~ cum
    RD  : EOS ~ rate + duration      <- separates the two ingredients of cum
All three are fitted WITHIN SITES (year minus site mean, sites with >= 3
years, SE clustered by site) and compared by AIC (identical rows), within-site
R2, leave-one-site-out CV R2, and the betas in RD (days per +1 within-site
SD). "Growth rate matters more" is supported if
|beta_rate| >> |beta_duration| in RD and R beats C out-of-sample.
`rate_negative` flags the Notion hypothesis "rate is always negative".

Interpretation caveat - anchors:
  * anchor 'fixed' (site-mean EOS as window end): windows that START at the
    solstice have constant length within a site, so cum and rate differ only
    by a site constant and cannot be told apart within-site. They are skipped
    in 'fixed' mode. SOS-anchored windows remain informative because
    duration varies with leaf-out.
  * anchor 'year' (this year's own EOS): reproduces the earlier analysis. There
    duration is tied to the target (longer season = later EOS), so a positive
    cum effect is expected mechanically; compare with 'fixed'.

Groups: 'ALL' plus every IGBP class with >= MIN_SITES sites. Optionally supply
data/site_leafout_strategy.csv (site_id, group) - e.g. deterministic vs
non-deterministic species after Baumgarten et al. 2026 - used as extra groups.

Input : data/eos_window_predictors_{fixed,year}_anchor.csv   (step 27)
        data/site_leafout_strategy.csv                        (optional)
Output: data/eos_rate_vs_cumulative.csv
"""
import numpy as np
import pandas as pd
import eos_common as ec

OUT_CSV = ec.DATA_DIR / "eos_rate_vs_cumulative.csv"
TARGETS_17 = ['EOS10']          # Notion: focus on EOS10 first; add 'EOS50' here if wanted
STRATEGY_CSV = ec.DATA_DIR / "site_leafout_strategy.csv"
WINDOWS = ['SOS_to_SOL', 'SOS_to_EOS90', 'SOS_to_EOS50', 'SOS_to_EOS10', 'SOL_to_EOS90', 'SOL_to_EOS10']
ec.require(ec.WINDOW_FIXED_CSV, ec.WINDOW_YEAR_CSV, hint="Run 27_window_predictors.py first.")

strategy = pd.read_csv(STRATEGY_CSV) if STRATEGY_CSV.exists() else None
if strategy is not None:
    print(f"Using site groups from '{STRATEGY_CSV.name}'.")

rows = []
for mode, path in (('fixed', ec.WINDOW_FIXED_CSV), ('year', ec.WINDOW_YEAR_CSV)):
    w = pd.read_csv(path)
    if strategy is not None:
        w = w.merge(strategy.rename(columns={'group': 'strategy_group'}), on='site_id', how='left')
    group_defs = [('ALL', 'ALL', w)]
    for col, label in (('igbp', 'IGBP'), ('strategy_group', 'strategy')):
        if col in w.columns:
            for g, sub in w.dropna(subset=[col]).groupby(col):
                group_defs.append((f'{label}:{g}', label, sub))
    carbon_vars = sorted({c.split('_cum__')[0] for c in w.columns if '_cum__' in c} & {'GPP', 'NPP', 'NPPd'})

    for gname, _, wg in group_defs:
        for vi in sorted(wg['vi_index'].unique()):
            for target in TARGETS_17:
                for cv in carbon_vars:
                    for win in WINDOWS:
                        if mode == 'fixed' and win.startswith('SOL_'):
                            continue
                        rate, cum, ln = f'{cv}_mean__{win}', f'{cv}_cum__{win}', f'len__{win}'
                        if rate not in wg.columns:
                            continue
                        d = wg[wg['vi_index'] == vi][[target, 'site_id', rate, cum, ln]].dropna()
                        d = d[d.groupby('site_id')['site_id'].transform('size') >= ec.MIN_YEARS_WITHIN].reset_index(drop=True)
                        if len(d) < ec.MIN_OBS or d['site_id'].nunique() < ec.MIN_SITES:
                            continue
                        models = {'R_rate': [rate], 'C_cum': [cum], 'RD_rate+duration': [rate, ln]}
                        rec = {'anchor': mode, 'group': gname, 'vi_index': vi, 'target': target, 'carbon': cv,
                               'window': win, 'n_obs': len(d), 'n_sites': d['site_id'].nunique()}
                        ok = True
                        for mname, xs in models.items():
                            f = ec.fit_within(d, target, xs)
                            if f is None:
                                ok = False
                                break
                            rec[f'aic_{mname}'] = f.aic
                            rec[f'r2w_{mname}'] = float(f.rsquared)
                            rec[f'loso_r2_{mname}'] = ec.loso_within(d, target, xs)['loso_r2']
                            ct = ec.within_coef_table(f, xs)
                            if mname == 'RD_rate+duration':
                                rec['beta_rate_in_RD'], rec['p_rate_in_RD'] = ct.iloc[0][['beta_days_per_sd', 'p_value']]
                                rec['beta_duration_in_RD'], rec['p_duration_in_RD'] = ct.iloc[1][['beta_days_per_sd', 'p_value']]
                            elif mname == 'R_rate':
                                rec['beta_rate'], rec['p_rate'] = ct.iloc[0][['beta_days_per_sd', 'p_value']]
                            else:
                                rec['beta_cum'], rec['p_cum'] = ct.iloc[0][['beta_days_per_sd', 'p_value']]
                        if not ok:
                            continue
                        rec['rate_negative'] = rec['beta_rate'] < 0
                        rec['rate_beats_cum_loso'] = (rec['loso_r2_R_rate'] > rec['loso_r2_C_cum']
                                                      if np.isfinite(rec['loso_r2_R_rate']) else np.nan)
                        rec['rate_over_duration_ratio'] = (abs(rec['beta_rate_in_RD']) /
                                                            max(abs(rec['beta_duration_in_RD']), 1e-9))
                        rows.append(rec)

res = pd.DataFrame(rows)
res.to_csv(OUT_CSV, index=False)
print(f"{len(res)} comparisons -> '{OUT_CSV}'")
if res.empty:
    raise SystemExit("No comparison had enough data.")

pd.set_option('display.width', 220)
show = res[(res['group'] == 'ALL') & (res['target'] == 'EOS10')]
print("\nALL sites, target EOS10 (within-site beta = days per +1 SD):")
print(show[['anchor', 'vi_index', 'carbon', 'window', 'n_obs', 'beta_rate', 'beta_cum',
            'loso_r2_R_rate', 'loso_r2_C_cum', 'rate_negative', 'rate_over_duration_ratio']]
      .round(3).to_string(index=False))
print("\nShare of comparisons where the rate coefficient is negative, by anchor mode:")
print(res.groupby('anchor')['rate_negative'].mean().round(2).to_string())
