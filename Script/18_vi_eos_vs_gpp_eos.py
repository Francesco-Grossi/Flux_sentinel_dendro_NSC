"""
PIPELINE STEP 18 - Which vegetation index gives the best EOS for linking with
carbon fluxes? (Notion C.3, Wang et al. 2024 hypothesis: NIRv is probably a
better proxy for EOS10 than NDVI, since NDVI mostly carries canopy structure
while NIRv also carries physiology.)

Objective reference: EOS10 derived from the flux tower's own GPP curve, with
the same amplitude convention as the pipeline's EOS10 (first day after the
seasonal GPP peak on which the smoothed GPP has fallen to 10% of the seasonal
amplitude, i.e. near dormancy):
    smooth = 15-day centred rolling mean of daily GPP
    base   = 5th percentile of smooth over the year; amp = peak - base
    frac   = (smooth - base) / amp
Only EOS10 is derived from GPP. GPP-derived EOS90 / EOS50 were dropped: they
fall in mid-July / early September (end of peak photosynthesis, not leaf
senescence), 40-60 days before every canopy-based source, and do not track
PhenoCam or tower NDVI from year to year (within-site r <= 0.15, step 21).

For each vi_index the VI-derived EOS10 is compared with the
GPP-derived EOS10 of the same site-year: r, within-site r, bias (VI - GPP),
MAE, RMSE, n. Lower RMSE / higher r = the index tracks physiological
senescence better.

Input : data/phenology_double_logistic_by_site_year_index.csv (step 5)
        data/fluxnet_landsat_merged.csv                        (step 4)
Output: data/gpp_derived_eos_by_site_year.csv
        data/vi_eos_vs_gpp_eos.csv
        figure/vi_eos_vs_gpp_eos/<level>.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eos_common as ec

LEVELS = (10,)
SMOOTH_DAYS = 15
MIN_AMPLITUDE = 1.0        # gC m-2 d-1; skip site-years without a real seasonal cycle
MIN_VALID_DAYS = 240
OUT_GPP_EOS = ec.DATA_DIR / "gpp_derived_eos_by_site_year.csv"
OUT_CMP = ec.DATA_DIR / "vi_eos_vs_gpp_eos.csv"
FIG_DIR = ec.FIGURE_DIR / "vi_eos_vs_gpp_eos"
FIG_DIR.mkdir(parents=True, exist_ok=True)

pheno = ec.load_phenology()
flux = ec.load_flux_daily()
vis = sorted(pheno['vi_index'].unique())
print(f"vi_index in phenology table: {vis}")

recs = []
for (site, year), g in flux[['site_id', 'year', 'doy', 'GPP']].groupby(['site_id', 'year']):
    s = g.dropna(subset=['GPP']).drop_duplicates('doy').set_index('doy')['GPP'].reindex(np.arange(1, 367))
    if s.notna().sum() < MIN_VALID_DAYS:
        continue
    sm = s.interpolate(limit=7, limit_area='inside').rolling(SMOOTH_DAYS, center=True, min_periods=SMOOTH_DAYS // 2).mean()
    if sm.notna().sum() < MIN_VALID_DAYS:
        continue
    base, peak_doy = np.nanpercentile(sm, 5), int(sm.idxmax())
    amp = sm.max() - base
    if not np.isfinite(amp) or amp < MIN_AMPLITUDE:
        continue
    frac = (sm - base) / amp
    after = frac.loc[peak_doy:]
    rec = {'site_id': site, 'year': int(year), 'gpp_peak_doy': peak_doy, 'gpp_amplitude': float(amp)}
    for lvl in LEVELS:
        hit = after[after <= lvl / 100.0]
        rec[f'GPP_EOS{lvl}'] = float(hit.index[0]) if len(hit) else np.nan
    recs.append(rec)

gpp_eos = pd.DataFrame(recs)
gpp_eos.to_csv(OUT_GPP_EOS, index=False)
print(f"GPP-derived EOS for {len(gpp_eos)} site-years -> '{OUT_GPP_EOS}'")
if gpp_eos.empty:
    raise SystemExit("No site-year had a usable GPP seasonal cycle.")

m = pheno.merge(gpp_eos, on=['site_id', 'year'], how='inner')
rows = []
for vi in vis:
    for lvl in LEVELS:
        d = m[m['vi_index'] == vi][['site_id', f'EOS{lvl}', f'GPP_EOS{lvl}']].dropna()
        if len(d) < ec.MIN_OBS or d['site_id'].nunique() < ec.MIN_SITES:
            continue
        x, y = d[f'GPP_EOS{lvl}'], d[f'EOS{lvl}']
        dm = d[[f'GPP_EOS{lvl}', f'EOS{lvl}']] - d.groupby('site_id')[[f'GPP_EOS{lvl}', f'EOS{lvl}']].transform('mean')
        rows.append({'vi_index': vi, 'level': f'EOS{lvl}', 'n_obs': len(d), 'n_sites': d['site_id'].nunique(),
                     'r_pooled': float(np.corrcoef(x, y)[0, 1]),
                     'r_within_site': float(np.corrcoef(dm.iloc[:, 0], dm.iloc[:, 1])[0, 1]) if dm.std().min() > 0 else np.nan,
                     'bias_days': float((y - x).mean()), 'mae_days': float((y - x).abs().mean()),
                     'rmse_days': float(np.sqrt(((y - x) ** 2).mean()))})
cmp_df = pd.DataFrame(rows)
cmp_df.to_csv(OUT_CMP, index=False)
if cmp_df.empty:
    raise SystemExit("No VI x level had enough matched site-years.")
print(f"\nVI-derived EOS vs GPP-derived EOS -> '{OUT_CMP}'")
print(cmp_df.round(3).to_string(index=False))
for lvl in cmp_df['level'].unique():
    best = cmp_df[cmp_df['level'] == lvl].sort_values('rmse_days').iloc[0]
    print(f"  lowest RMSE for {lvl}: {best['vi_index']} ({best['rmse_days']:.1f} d)")

colors = dict(zip(vis, ('#2b6cb0', '#38a169', '#dd6b20', '#805ad5')))
for lvl in LEVELS:
    sub = cmp_df[cmp_df['level'] == f'EOS{lvl}']
    if sub.empty:
        continue
    fig, axes = plt.subplots(1, len(sub), figsize=(4.4 * len(sub), 4.2), squeeze=False)
    for ax, (_, r) in zip(axes[0], sub.iterrows()):
        d = m[m['vi_index'] == r['vi_index']][[f'EOS{lvl}', f'GPP_EOS{lvl}']].dropna()
        ax.scatter(d[f'GPP_EOS{lvl}'], d[f'EOS{lvl}'], s=10, alpha=0.5, color=colors.get(r['vi_index'], 'gray'))
        lim = [min(d.min()), max(d.max())]
        ax.plot(lim, lim, 'k--', lw=0.8)
        ax.set_xlabel(f'GPP-derived EOS{lvl} (DOY)')
        ax.set_ylabel(f'{r["vi_index"]}-derived EOS{lvl} (DOY)')
        ax.set_title(f'{r["vi_index"]}: r={r["r_pooled"]:.2f}, RMSE={r["rmse_days"]:.0f} d, n={int(r["n_obs"])}', fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f'EOS{lvl}.png', dpi=140)
    plt.close(fig)
print(f"Figures -> '{FIG_DIR}'")
