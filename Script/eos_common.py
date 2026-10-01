"""
Shared helpers for pipeline steps 27-39 (EOS <-> carbon source/sink analysis,
from the Notion page "Ideas in Sep, 2026"). NOT a pipeline step itself -
run_pipeline.sh does not call it; steps 27-39 import it as `eos_common`.

Conventions kept identical to legacy steps 06-09:
  * site-year-index rows, keyed (site_id, year, vi_index)
  * mixed-effects models with a site random intercept (statsmodels mixedlm)
  * predictors z-scored, so betas are "days of EOS shift per 1 SD"
  * leave-one-site-out CV predicts with FIXED EFFECTS ONLY

Test/CI override: set PHENO_DATA_DIR / PHENO_FIGURE_DIR to redirect I/O.
Set PHENO_VI="NDVI,GCC" to restrict which vi_index values are analysed
(default: every vi_index present in the phenology tables of steps 23, 24, 25).
"""
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("PHENO_DATA_DIR", ROOT / "data"))
FIGURE_DIR = Path(os.environ.get("PHENO_FIGURE_DIR", ROOT / "figure"))
DATA_DIR.mkdir(parents=True, exist_ok=True)

PHENOLOGY_CSV = DATA_DIR / "phenology_double_logistic_by_site_year_index.csv"   # step 23
# OPTIONAL extra EOS sources, same columns as the step-23 table; appended when present
# (vi_index = NDVI_tower / GCC), so steps 27-39 analyse them too.
EXTRA_PHENOLOGY_CSVS = [DATA_DIR / "phenology_tower_by_site_year_index.csv",       # step 24
                        DATA_DIR / "phenology_phenocam_by_site_year_index.csv"]    # step 25
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"                               # step 22
WINDOW_FIXED_CSV = DATA_DIR / "eos_window_predictors_fixed_anchor.csv"           # step 27
WINDOW_YEAR_CSV = DATA_DIR / "eos_window_predictors_year_anchor.csv"             # step 27
# OPTIONAL input from step 26 (Luo et al. 2025 CUE): site_id, date and any of
#   NPP        = annual CUE x GPP            -> 'NPP'   (gC m-2 d-1)
#   NPPd       = daily (seasonal) CUE x GPP  -> 'NPPd'  (gC m-2 d-1)
#   CUE_daily  = daily (seasonal) CUE        -> 'CUEd'  (ratio: only window MEANS are meaningful)
# Coarser cadence is linearly interpolated.
NPP_COLUMNS = {'NPP': 'NPP', 'NPPd': 'NPPd', 'CUE_daily': 'CUEd'}
CARBON_FLUX_VARS = ('GPP', 'NEP', 'NPP', 'NPPd')   # cumulative and mean windows
RATIO_VARS = ('CUEd',)                             # mean windows only
NPP_CSV = DATA_DIR / "npp_luo2025_daily.csv"

MIN_FIT_CORR = 0.8           # quality gate on the double-logistic fit (as legacy steps 6-9)
MIN_OBS, MIN_SITES = 20, 5   # as legacy steps 8/9
MIN_COVERAGE = 0.8           # min fraction of window days with valid flux data
MIN_WINDOW_DAYS = 5
MIN_YEARS_ANCHOR = 3         # min years to compute a per-site mean EOS anchor
TARGETS = ['EOS10', 'EOS50', 'EOS90']
PRIMARY_TARGETS = ['EOS10', 'EOS50']   # Notion: "may first focus on EOS10 at this stage"

# short name -> (FLUXNET column, plausible range on the raw column, sign applied after QC)
FLUX_VARS = {
    'GPP': ('GPP_NT_VUT_REF', (-5, 50), 1),
    'NEP': ('NEE_VUT_REF', (-50, 50), -1),      # NEP = -NEE
    'TA':  ('TA_F', (-60, 50), 1),
    'SW':  ('SW_IN_F', (0, 500), 1),
    'VPD': ('VPD_F', (0, 100), 1),
    'P':   ('P_F', (0, 300), 1),
}


# ---------------------------------------------------------------- I/O helpers
def require(*paths, hint=""):
    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing '{p}'. {hint}")


def load_phenology():
    require(PHENOLOGY_CSV, hint="Run 23_phenology_satellite.py first.")
    pheno = pd.concat([pd.read_csv(p) for p in [PHENOLOGY_CSV] + EXTRA_PHENOLOGY_CSVS if os.path.exists(p)],
                      ignore_index=True)
    pheno = pheno[(pheno['method'] == 'double_logistic') & (pheno['corr'] >= MIN_FIT_CORR)].copy()
    if 'qc_pass' in pheno.columns:  # step-23 QC flags (autumn coverage, R2, date order, ...)
        pheno = pheno[pheno['qc_pass'].astype(bool)].copy()
    only = os.environ.get("PHENO_VI")
    if only:
        pheno = pheno[pheno['vi_index'].isin([v.strip() for v in only.split(',')])]
    pheno['year'] = pheno['year'].astype(int)
    return pheno.reset_index(drop=True)


def solstice_doy(year):
    return pd.Timestamp(year=int(year), month=6, day=21).dayofyear


def _parse_dates(df, col):
    d = pd.to_datetime(df[col].astype(str), format='%Y%m%d', errors='coerce')
    return d.fillna(pd.to_datetime(df[col].astype(str), errors='coerce'))


def load_flux_daily():
    """Daily flux table with short-named columns: GPP, NEP, TA, SW, VPD, P
    (+ NPP, NPPd, CUEd if NPP_CSV exists), plus site_id, year, doy, lat, igbp."""
    require(FLUX_CSV, hint="Run 22_merge_fluxnet_hls.py first.")
    raw = pd.read_csv(FLUX_CSV, low_memory=False)
    raw['date'] = _parse_dates(raw, 'TIMESTAMP')
    raw = raw.dropna(subset=['date']).copy()
    out = pd.DataFrame({'site_id': raw['site_id'].to_numpy(),
                        'year': raw['date'].dt.year.to_numpy(),
                        'doy': raw['date'].dt.dayofyear.to_numpy()})
    for meta in ('lat', 'igbp'):
        if meta in raw.columns:
            out[meta] = raw[meta].to_numpy()
    for short, (col, (lo, hi), sign) in FLUX_VARS.items():
        if col not in raw.columns:
            print(f"  [note] '{col}' not in merged table -> '{short}' skipped.")
            continue
        v = pd.to_numeric(raw[col], errors='coerce')
        v = v.where((v >= lo) & (v <= hi))
        out[short] = (sign * v).to_numpy()

    if os.path.exists(NPP_CSV):
        npp = pd.read_csv(NPP_CSV)
        npp = npp.rename(columns=NPP_COLUMNS)
        vals = [c for c in NPP_COLUMNS.values() if c in npp.columns]
        npp['date'] = pd.to_datetime(npp['date'], errors='coerce')
        npp = npp.dropna(subset=['date']).dropna(subset=vals, how='all').copy()
        npp['year'], npp['doy'] = npp['date'].dt.year, npp['date'].dt.dayofyear
        step = npp.sort_values(['site_id', 'date']).groupby(['site_id', 'year'])['doy'].diff().median()
        if pd.notna(step) and step > 1:
            print(f"  [note] NPP cadence ~{step:.0f} d -> linearly interpolated to daily.")
            parts = []
            for (s, y), g in npp.groupby(['site_id', 'year']):
                g = g.sort_values('doy').drop_duplicates('doy')
                idx = np.arange(int(g['doy'].min()), int(g['doy'].max()) + 1)
                part = {'site_id': s, 'year': y, 'doy': idx}
                for c in vals:
                    ok = g[c].notna()
                    part[c] = np.interp(idx, g.loc[ok, 'doy'], g.loc[ok, c]) if ok.sum() > 1 else np.nan
                parts.append(pd.DataFrame(part))
            npp = pd.concat(parts, ignore_index=True)
        out = out.merge(npp[['site_id', 'year', 'doy'] + vals].drop_duplicates(['site_id', 'year', 'doy']),
                        on=['site_id', 'year', 'doy'], how='left')
        print(f"  [note] {vals} loaded from '{NPP_CSV.name}'.")
    else:
        print(f"  [note] '{NPP_CSV.name}' not found -> NPP (sink) analyses skipped, GPP/NEP only.")
    return out


class FluxLookup:
    """Fast per-(site, year) window statistics on the daily table."""

    def __init__(self, flux, cols):
        self.cols = [c for c in cols if c in flux.columns]
        self._d = {}
        for key, g in flux.groupby(['site_id', 'year'], sort=False):
            g = g.sort_values('doy')
            self._d[(key[0], int(key[1]))] = (
                g['doy'].to_numpy(), {c: g[c].to_numpy(dtype=float) for c in self.cols})

    @staticmethod
    def n_days(start, end):
        if not (np.isfinite(start) and np.isfinite(end)):
            return np.nan
        n = int(np.floor(end)) - int(np.ceil(start)) + 1
        return float(n) if n >= MIN_WINDOW_DAYS else np.nan

    def window_mean(self, site, year, col, start, end, min_cov=MIN_COVERAGE):
        """Mean over [start, end] (DOY) of valid days; NaN if coverage < min_cov."""
        n_exp = self.n_days(start, end)
        if np.isnan(n_exp):
            return np.nan
        rec = self._d.get((site, int(year)))
        if rec is None or col not in rec[1]:
            return np.nan
        doy, arrs = rec
        s, e = int(np.ceil(start)), int(np.floor(end))
        v = arrs[col][(doy >= s) & (doy <= e)]
        v = v[np.isfinite(v)]
        if len(v) < min_cov * n_exp:
            return np.nan
        return float(v.mean())


# --------------------------------------------------------------- statistics
def zscore(s):
    return (s - s.mean()) / s.std()


def fit_lme(d, y, xs, reml=True):
    """Site-random-intercept LME of y on z-scored xs. Predictors are renamed
    v0..vk internally (robust to odd column names); returns fit or None."""
    d = d[[y, 'site_id'] + list(xs)].dropna()
    if len(d) < 5 or (d[xs].std() == 0).any() or np.isnan(d[xs].std()).any():
        return None
    X = d[xs].apply(zscore)
    X.columns = [f'v{i}' for i in range(len(xs))]
    dd = pd.concat([d[[y, 'site_id']].reset_index(drop=True), X.reset_index(drop=True)], axis=1)
    formula = f"{y} ~ " + " + ".join(X.columns)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return smf.mixedlm(formula, dd, groups=dd['site_id']).fit(reml=reml)
    except Exception:
        return None


def _bse(fit, k):
    """Fixed-effect SE; NaN (not a warning) when the fit's covariance is not positive definite."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return float(fit.bse_fe[k])


def coef_table(fit, xs):
    """Fixed-effect table keyed by the original predictor names."""
    rows = []
    for i, name in enumerate(xs):
        k = f'v{i}'
        rows.append({'predictor': name, 'beta_days_per_sd': float(fit.fe_params[k]),
                     'std_err': _bse(fit, k), 'p_value': float(fit.pvalues[k])})
    return pd.DataFrame(rows)


def lme_slope(d, y, x):
    """Single-predictor LME. Returns dict or None if data are insufficient."""
    d = d[[y, x, 'site_id']].dropna()
    n, ns = len(d), d['site_id'].nunique()
    if n < MIN_OBS or ns < MIN_SITES:
        return None
    fit = fit_lme(d, y, [x])
    if fit is None:
        return None
    r_pool = pearsonr(d[x], d[y])[0]
    dm = d[[x, y]] - d.groupby('site_id')[[x, y]].transform('mean')
    r_within = pearsonr(dm[x], dm[y])[0] if dm[x].std() > 0 and dm[y].std() > 0 else np.nan
    return {'beta_days_per_sd': float(fit.fe_params['v0']), 'std_err': _bse(fit, 'v0'),
            'p_value': float(fit.pvalues['v0']), 'r_pooled': float(r_pool),
            'r_within_site': float(r_within), 'n_obs': n, 'n_sites': ns}


def r2_nakagawa(fit):
    """(marginal, conditional) R2 of a random-intercept LME."""
    fe = fit.fe_params.to_numpy()
    var_f = float(np.var(fit.model.exog @ fe))
    var_r = float(np.asarray(fit.cov_re).reshape(-1)[0])
    var_e = float(fit.scale)
    tot = var_f + var_r + var_e
    return var_f / tot, (var_f + var_r) / tot


def max_vif(d, xs):
    if len(xs) < 2:
        return np.nan
    try:
        return float(np.max(np.diag(np.linalg.inv(np.corrcoef(d[xs].to_numpy().T)))))
    except np.linalg.LinAlgError:
        return np.inf


def loso_cv(d, y, xs):
    """Leave-one-site-out CV, fixed effects only (same protocol as legacy steps 8/9)."""
    d = d.reset_index(drop=True)
    preds = np.full(len(d), np.nan)
    for site in d['site_id'].unique():
        te = (d['site_id'] == site).to_numpy()
        train, test = d[~te], d[te]
        mu, sd = train[xs].mean(), train[xs].std()
        if (sd == 0).any() or sd.isna().any() or train['site_id'].nunique() < 3:
            continue
        Xtr = ((train[xs] - mu) / sd)
        Xtr.columns = [f'v{i}' for i in range(len(xs))]
        tr = pd.concat([train[[y, 'site_id']].reset_index(drop=True), Xtr.reset_index(drop=True)], axis=1)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fit = smf.mixedlm(f"{y} ~ " + " + ".join(Xtr.columns), tr, groups=tr['site_id']).fit(reml=True)
        except Exception:
            continue
        Xte = np.column_stack([np.ones(len(test))] + [((test[c] - mu[c]) / sd[c]).to_numpy() for c in xs])
        preds[te] = Xte @ fit.fe_params.to_numpy()
    ok = ~np.isnan(preds)
    if ok.sum() < 5:
        return {'n_predicted': int(ok.sum()), 'loso_rmse_days': np.nan, 'loso_r2': np.nan}
    yt = d[y].to_numpy()[ok]
    res = yt - preds[ok]
    ss_tot = float(np.sum((yt - yt.mean()) ** 2))
    return {'n_predicted': int(ok.sum()), 'loso_rmse_days': float(np.sqrt(np.mean(res ** 2))),
            'loso_r2': float(1 - np.sum(res ** 2) / ss_tot) if ss_tot > 0 else np.nan}


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
