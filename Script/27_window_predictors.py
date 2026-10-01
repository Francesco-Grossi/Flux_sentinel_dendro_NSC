"""
PIPELINE STEP 12 - Build carbon and environment WINDOW predictors for the
EOS <-> carbon source/sink analyses (steps 13-18).

For every site-year-index row of the step-5 phenology table, and for each
window below, compute from the daily FLUXNET data:
    <VAR>_mean__<window>   mean daily value ("rate", e.g. GPPmean / GPPrate)
    <VAR>_cum__<window>    cumulative value = mean x window length (cGPP)
    len__<window>          window length in days
for carbon variables VAR in {GPP, NEP, NPP*} (*NPP only if the optional
Luo et al. 2025 file data/npp_luo2025_daily.csv exists), and for the
environment (TA, SW, VPD mean; P sum) over the pre-/post-solstice windows.

Windows (anchors: SOS = leaf_out_10, SOL = summer solstice, EOSxx = end of season):
    SOS_to_SOL  SOL_to_EOS90  SOL_to_EOS10  SOS_to_EOS90  SOS_to_EOS50
    SOS_to_EOS10  EOS90_to_EOS50

TWO anchor modes are written (this is the key methodological point from the
Notion page - "EOS should be the fixed value in each site"):
  fixed : EOSxx anchors = each site's multi-year MEAN EOSxx (per vi_index).
          Window boundaries do not depend on the year's own EOS, so cumulative
          carbon is not mechanically tied to the target.
  year  : EOSxx anchors = that year's own EOSxx. This reproduces the earlier
          analysis (GPPcum over [SOS, EOSxx]); cumulative carbon there is
          partly a proxy for season length, which is what can make cGPP look
          "positively" related to EOS while GPPrate stays negative.

Input : data/phenology_double_logistic_by_site_year_index.csv   (step 5)
        data/fluxnet_landsat_merged.csv                          (step 4)
        data/npp_luo2025_daily.csv                               (optional)
Output: data/eos_window_predictors_fixed_anchor.csv
        data/eos_window_predictors_year_anchor.csv
"""
import numpy as np
import pandas as pd
import eos_common as ec

WINDOWS = {
    'SOS_to_SOL':     ('SOS', 'SOL'),
    'SOL_to_EOS90':   ('SOL', 'EOS90'),
    'SOL_to_EOS10':   ('SOL', 'EOS10'),
    'SOS_to_EOS90':   ('SOS', 'EOS90'),
    'SOS_to_EOS50':   ('SOS', 'EOS50'),
    'SOS_to_EOS10':   ('SOS', 'EOS10'),
    'EOS90_to_EOS50': ('EOS90', 'EOS50'),
}
ENV_WINDOWS = ['SOS_to_SOL', 'SOL_to_EOS10']     # "before / after summer solstice"
ENV_HOW = {'TA': 'mean', 'SW': 'mean', 'VPD': 'mean', 'P': 'sum'}
ANCHOR_STAT = 'mean'                              # Notion: "mean EOS10 across years"

print("Loading data ...")
pheno = ec.load_phenology()
flux = ec.load_flux_daily()
carbon_vars = [v for v in ec.CARBON_FLUX_VARS if v in flux.columns]
ratio_vars = [v for v in ec.RATIO_VARS if v in flux.columns]   # window mean only (e.g. seasonal CUE)
env_vars = [v for v in ENV_HOW if v in flux.columns]
lk = ec.FluxLookup(flux, carbon_vars + ratio_vars + env_vars)
print(f"Phenology rows: {len(pheno)} | vi_index: {sorted(pheno['vi_index'].unique())}")
print(f"Carbon variables: {carbon_vars} | ratios: {ratio_vars} | environment: {env_vars}")

meta = flux.groupby('site_id').agg(**{c: (c, 'first') for c in ('lat', 'igbp') if c in flux.columns})

for lvl in (10, 50, 90):
    col = f'EOS{lvl}'
    g = pheno.groupby(['site_id', 'vi_index'])[col]
    pheno[f'{col}_fixed'] = g.transform(ANCHOR_STAT).where(g.transform('count') >= ec.MIN_YEARS_ANCHOR)


def build(anchor_mode):
    recs = []
    for _, r in pheno.iterrows():
        site, year = r['site_id'], int(r['year'])
        anchors = {'SOS': r['leaf_out_10'], 'SOL': float(ec.solstice_doy(year))}
        for lvl in (10, 50, 90):
            anchors[f'EOS{lvl}'] = r[f'EOS{lvl}_fixed'] if anchor_mode == 'fixed' else r[f'EOS{lvl}']
        rec = {'site_id': site, 'year': year, 'vi_index': r['vi_index'], 'anchor_mode': anchor_mode,
               'SOS': r['leaf_out_10'], 'SOL': anchors['SOL'],
               'EOS10': r['EOS10'], 'EOS50': r['EOS50'], 'EOS90': r['EOS90'],
               'EOS10_anchor': anchors['EOS10'], 'EOS50_anchor': anchors['EOS50'],
               'EOS90_anchor': anchors['EOS90']}
        if site in meta.index:
            for c in meta.columns:
                rec[c] = meta.at[site, c]
        for wname, (a, b) in WINDOWS.items():
            s, e = anchors[a], anchors[b]
            n = ec.FluxLookup.n_days(s, e)
            rec[f'len__{wname}'] = n
            for v in carbon_vars:
                m = lk.window_mean(site, year, v, s, e)
                rec[f'{v}_mean__{wname}'] = m
                rec[f'{v}_cum__{wname}'] = m * n if np.isfinite(m) and np.isfinite(n) else np.nan
            for v in ratio_vars:
                rec[f'{v}_mean__{wname}'] = lk.window_mean(site, year, v, s, e)
            if wname in ENV_WINDOWS:
                for v in env_vars:
                    m = lk.window_mean(site, year, v, s, e)
                    if ENV_HOW[v] == 'sum':
                        rec[f'{v}_sum__{wname}'] = m * n if np.isfinite(m) and np.isfinite(n) else np.nan
                    else:
                        rec[f'{v}_mean__{wname}'] = m
        recs.append(rec)
    return pd.DataFrame(recs)


for mode, path in (('fixed', ec.WINDOW_FIXED_CSV), ('year', ec.WINDOW_YEAR_CSV)):
    df = build(mode)
    df.to_csv(path, index=False)
    print(f"\n[{mode}] -> '{path}' ({len(df)} rows). Non-null counts:")
    for wname in WINDOWS:
        parts = [f"{v}_cum={int(df[f'{v}_cum__{wname}'].notna().sum())}" for v in carbon_vars]
        print(f"   {wname:15s} " + "  ".join(parts))
