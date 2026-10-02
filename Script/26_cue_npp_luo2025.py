"""
PIPELINE STEP 26 - Carbon use efficiency (CUE) from eddy-covariance data -
annual, by temperature, and seasonal / daily - and daily NPP = CUE x GPP for
the carbon-SINK analyses (steps 27 and 33-38).

Python port of the MATLAB code "CUE_fluxnet" v1 by Xiangzhong (Remi) Luo
(main_CUE_v1.m, readfluxnet_data.m, estimate_site_CUE.m, est_CUE_ss2.m,
est_CUE_fixed_ss2.m, and the DRAM sampler of mcmcstat/mcmcrun.m by M. Laine).
Reference: Luo et al., https://www.researchsquare.com/article/rs-3989566/v1

Method (as in the MATLAB code), per site, with GPP and Reco the mean of the
night-time and day-time partitioned fluxes:
  Round 1 - per year, 10-day windows every 5 days over the growing season
            (GPP > 20% of its range and Tair > 5 C). All pairs of days in a
            window give dReco, dGPP, dcumGPP. MCMC (DRAM, 200 samples) fits
                dReco = gR*CUE*dGPP + mR*CUE*dcumGPP*(1-tau)        (weight 0.05)
                dReco = (1-CUE)*dGPP                                 (weight 0.95)
            for gR, mR, CUE, tau. Their site means become fixed values /
            starting points for round 2.
  Round 2 - per year, 1-degree air-temperature bins (0-35 C, rain < 2 mm).
            Same pairs; gR and tau fixed, maintenance respiration scaled by
            an Arrhenius term exp(-Ea/k * (1/T - 1/T0)); MCMC fits CUE and Ea.
  Annual  - bootstrap (1000 draws from each bin's mean and MCMC sd) gives the
            year's gR, CUE, Ea, tau and their sd.

Seasonal / daily CUE (NOT in the MATLAB code; an extension of its round 2):
  Round 3 - per year, SEASON_WINDOW_DAYS-day windows every SEASON_STEP_DAYS
            days. Pairs of days inside a window are kept only if their air
            temperatures differ by <= SEASON_MAX_DT C (rain < 2 mm), so that,
            as in round 2, the pairs compare days of similar temperature and
            the seasonal cooling of soil respiration does not leak into CUE.
            gR, tau, mR0, T0 are fixed as in round 2 and Ea is fixed to the
            year's value; MCMC fits CUE only.
  Daily   - the window values (3-window running median) are interpolated
            linearly to every day of the year; before the first / after the
            last window the nearest value is held. A site-year needs >=
            SEASON_MIN_WINDOWS converged windows, else its daily CUE is NaN.
  This is noisier than the annual value and has not been validated against
  the paper; read it as a seasonal shape, not as exact daily numbers.

INDEXING_MODE
  'corrected'  (default) fixes the two points below.
  'matlab'     reproduces the MATLAB code exactly as written, including:
               (a) round 1 applies the growing-season window mask (length =
                   number of growing-season days) to the FULL-YEAR arrays, so
                   it selects the first days of the calendar year, not the
                   growing-season days;
               (b) round 2, for a bin whose chain did not move, copies the
                   previous bin's CUE and Ea VALUES into the uncertainty
                   columns.
               'corrected' (a) uses the growing-season days the window refers
               to, and (b) copies the previous bin's uncertainties.
Random numbers differ from MATLAB's, so results match only statistically.
The year-to-year pattern is stable across random seeds (r ~ 0.9+ at US-Ha1),
the absolute level less so (+/- ~0.1): use within-site differences.

Flux partitioning: MATLAB needs both NT and DT columns. Where a file has only
one of them, that one is used and flux_partition records which ('NT+DT',
'NT' or 'DT').

Input : data_raw/*_<site>_FLUXNET_*.zip (daily DD file), for the sites in
        data/fluxnet_landsat_merged.csv
Respiration terms (an extension; my reading of the model in the MATLAB cost
function, which the MATLAB code itself does not output):
    Rg = gR x CUE x GPP                         growth respiration
    Rm = mR0 x exp(-Ea/k (1/T - 1/T0)) x CUE x cumGPP x (1 - tau)
                                                maintenance respiration of the
         biomass built since 1 January, at the day's temperature. It does not
         contain the day's GPP, unlike (1 - CUE) x GPP.
         The model is fitted on differences between days, so this expression
         fixes the seasonal COURSE of Rm, not its level (taken literally it
         exceeds GPP). Rm is therefore scaled per site-year so that
         Rg + Rm = (1 - CUE) x GPP over the year: the year's autotrophic
         respiration is distributed over the days by temperature and
         accumulated biomass instead of by the day's GPP.
mR0 and T0 are site constants from round 1 and are saved with the yearly table.

Saved per site: each site's result is written to data_raw/_cue_cache/<site>.pkl
as soon as it is finished. A rerun reads those back and computes only the
sites that are missing, or whose input archive or method settings changed.
An interrupted run therefore loses only the sites it was working on.
Set CUE_RESTART=1 to recompute everything.

Output: data/cue_luo2025_site_year.csv          site_id, year, gR, CUE, Ea, tau, *_sd, mR0, T0, ...
        data/cue_luo2025_temperature_bins.csv   CUE and Ea per 1-degree bin (round 2)
        data/cue_luo2025_seasonal.csv           CUE per sliding window (round 3)
        data/npp_luo2025_daily.csv              site_id, date, CUE_annual, CUE_daily,
                                                NPP  = CUE_annual x GPP_NT_VUT_REF
                                                NPPd = CUE_daily  x GPP_NT_VUT_REF
                                                Rg, Rm (see above)
Usage : python 26_cue_npp_luo2025.py                     all pipeline sites
        python 26_cue_npp_luo2025.py path/to/FLX_..._DD_....csv [more.csv]   only these files
        CUE_INDEXING=matlab python 26_cue_npp_luo2025.py      (writes *_matlab.csv)
"""
import glob
import io
import os
import re
import sys
import zipfile
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
RAW_DIR = ROOT / "data_raw"
FLUX_CSV = DATA_DIR / "fluxnet_landsat_merged.csv"

INDEXING_MODE = os.environ.get("CUE_INDEXING", "corrected")  # 'corrected' | 'matlab'
SUFFIX = "" if INDEXING_MODE == "corrected" else f"_{INDEXING_MODE}"
OUT_CUE = DATA_DIR / f"cue_luo2025_site_year{SUFFIX}.csv"
OUT_BINS = DATA_DIR / f"cue_luo2025_temperature_bins{SUFFIX}.csv"
OUT_SEASON = DATA_DIR / f"cue_luo2025_seasonal{SUFFIX}.csv"
OUT_NPP = DATA_DIR / f"npp_luo2025_daily{SUFFIX}.csv"
N_WORKERS = max(1, (os.cpu_count() or 2) - 1)
CACHE_DIR = RAW_DIR / f"_cue_cache{SUFFIX}"      # one file per finished site
CACHE_VERSION = 2                                 # raise when the method changes, to recompute every site
RESTART = os.environ.get("CUE_RESTART", "0") == "1"

NSIMU = 200                 # options.nsimu
ADAPT_INT = 100             # mcmcrun default
DR_SCALE = 5.0              # mcmcrun default, second-stage proposal shrink
N_BOOT = 1000
# round 3 (seasonal CUE)
SEASON_WINDOW_DAYS = 30
SEASON_STEP_DAYS = 5
SEASON_MAX_DT = 1.5         # max air-temperature difference within a pair of days (C)
SEASON_MIN_DAYS = 10
SEASON_MIN_PAIRS = 20
SEASON_MIN_WINDOWS = 4      # converged windows needed to build a daily series for a site-year
CUE_BOUNDS = (0.1, 1.0)
BOLTZMANN_EV = 8.617333262e-5
NEEDED = ['TIMESTAMP', 'TA_F', 'SW_IN_F', 'VPD_F', 'P_F', 'GPP_NT_VUT_REF', 'GPP_DT_VUT_REF',
          'RECO_NT_VUT_REF', 'RECO_DT_VUT_REF']
PARAMS = ['gR', 'CUE', 'Ea', 'tau']


# ------------------------------------------------------------------ mcmcrun (DRAM)
def dram(ssfun, p0, low, upp, rng, nsimu=NSIMU, adaptint=ADAPT_INT):
    """Delayed-rejection adaptive Metropolis as mcmcstat's mcmcrun with its
    defaults: sigma2 = 1, flat priors inside [low, upp], 2 tries, initial
    proposal sd = 5% of the start value, adaptation every `adaptint` samples.
    Returns (chain mean, chain)."""
    p0, low, upp = (np.asarray(v, float) for v in (p0, low, upp))
    npar = len(p0)
    qvar = (np.abs(p0) * 0.05) ** 2
    qvar[qvar == 0] = 1.0
    R = np.diag(np.sqrt(qvar))                      # upper Cholesky factor of the proposal covariance
    invR = np.linalg.inv(R)
    R2 = R / DR_SCALE
    qcov_scale = 2.4 / np.sqrt(npar)

    chain = np.empty((nsimu, npar))
    old, ss = p0.copy(), ssfun(p0)
    chain[0] = old
    for i in range(1, nsimu):
        new = old + rng.standard_normal(npar) @ R
        if np.any(new < low) or np.any(new > upp):
            accept, tst, ss1 = False, 0.0, np.inf
        else:
            ss1 = ssfun(new)
            with np.errstate(over='ignore', invalid='ignore'):
                tst = np.exp(-0.5 * (ss1 - ss))
            accept = bool(tst >= 1) or bool(tst > 0 and tst > rng.random())   # NaN -> reject, as in MATLAB
        if not accept:                              # delayed rejection, second try
            z = old + rng.standard_normal(npar) @ R2
            if not (np.any(z < low) or np.any(z > upp)):
                ssz = ssfun(z)
                with np.errstate(over='ignore', invalid='ignore'):
                    a1 = 1.0 - min(1.0, tst)                               # 1 - alpha(x, y)
                    a2 = 1.0 - min(1.0, np.exp(-0.5 * (ss1 - ssz)))        # 1 - alpha(z, y)
                    if a2 == 0:
                        alpha = 0.0
                    else:
                        q = -0.5 * (np.sum(((new - z) @ invR) ** 2) - np.sum(((new - old) @ invR) ** 2))
                        alpha = min(1.0, np.exp(-0.5 * (ssz - ss) + q) * a2 / a1)
                if alpha >= 1 or rng.random() < alpha:
                    accept, new, ss1 = True, z, ssz
        if accept:
            old, ss = new, ss1
        chain[i] = old

        if (i + 1) % adaptint == 0:                 # adapt the proposal to the chain so far
            cov = np.atleast_2d(np.cov(chain[:i + 1], rowvar=False))
            for jitter in (0.0, 1e-8):
                try:
                    R = np.linalg.cholesky(cov + np.eye(npar) * jitter).T * qcov_scale
                    invR, R2 = np.linalg.inv(R), R / DR_SCALE
                    break
                except np.linalg.LinAlgError:
                    continue
    return chain.mean(axis=0), chain


def pair_diffs(v):
    """v[j] - v[i] for all pairs i < j (MATLAB combntns(v, 2), second minus first)."""
    i, j = np.triu_indices(len(v), k=1)
    return v[j] - v[i]


# ------------------------------------------------------------------ estimate_site_CUE
def year_arrays(y):
    gpp = (0.5 * (y['GPP_DT'] + y['GPP_NT'])).to_numpy(float).copy()
    reco = (0.5 * (y['RECO_DT'] + y['RECO_NT'])).to_numpy(float)
    gpp[gpp < 0.01] = np.nan
    cumgpp = np.cumsum(np.nan_to_num(gpp))
    cli = y[['TA', 'SW', 'VPD', 'P']].to_numpy(float)
    return gpp, reco, cumgpp, cli, y['doy'].to_numpy(float)


def nanmean_rows(a):
    with np.errstate(invalid='ignore'):
        return np.nanmean(a, axis=0) if len(a) else np.full(a.shape[1], np.nan)


def estimate_site_cue(df, rng, mode=INDEXING_MODE):
    """df: daily rows with year, doy, TA, SW, VPD, P, GPP_NT, GPP_DT, RECO_NT, RECO_DT.
    Returns (annual, bins, seasonal): one row per year (gR, CUE, Ea, tau and
    their sd), per year x temperature bin, and per year x sliding window."""
    years = np.sort(df['year'].unique())
    by_year = {yr: year_arrays(df[df['year'] == yr]) for yr in years}

    # ---- round 1: gR, mR, CUE, tau per 10-day growing-season window
    out1 = []
    for yr in years:
        gpp, reco, cumgpp, cli, doy = by_year[yr]
        if np.all(np.isnan(gpp)):
            continue
        threshold = (np.nanmax(gpp) - np.nanmin(gpp)) * 0.2 + np.nanmin(gpp)
        with np.errstate(invalid='ignore'):
            ind_gs = (gpp > threshold) & (cli[:, 0] > 5)
        gs_doy = doy[ind_gs]
        if len(gs_doy) == 0:
            continue
        gs_pos = np.flatnonzero(ind_gs)
        for t_c in np.arange(gs_doy.min(), gs_doy.max() + 1e-9, 5):
            t_index = (gs_doy > t_c) & (gs_doy < t_c + 10)
            if t_index.sum() <= 5:
                continue
            # 'matlab': the mask over growing-season days indexes the full-year arrays from their start
            sel = np.flatnonzero(t_index) if mode == 'matlab' else gs_pos[t_index]
            dy, dx1, dx2 = pair_diffs(reco[sel]), pair_diffs(gpp[sel]), pair_diffs(cumgpp[sel])

            def ss1(k):
                model = k[0] * k[2] * dx1 + k[1] * k[2] * dx2 * (1 - k[3])
                return np.sum((dy - model) ** 2) * 0.05 + np.sum((dy - (1 - k[2]) * dx1) ** 2 * 0.95)

            mean, chain = dram(ss1, [0.20, 0.1, 0.5, 1.0], [0.15, 0.0, 0.1, 0.1], [0.30, 0.20, 1.0, 1.0], rng)
            out1.append([yr, t_c, *mean, *nanmean_rows(cli[sel]), *chain.std(axis=0)])

    cols = ['year'] + PARAMS + [f'{p}_sd' for p in PARAMS] + ['n_temperature_bins', 'n_bins_not_converged']
    bin_cols = ['year', 'tair_bin', 'CUE', 'Ea', 'CUE_sd', 'Ea_sd', 'converged']
    season_cols = ['year', 'doy_center', 'CUE', 'CUE_sd', 'n_days', 'n_pairs', 'tair_mean', 'converged']
    if not out1:
        return (pd.DataFrame({'year': years}).reindex(columns=cols), pd.DataFrame(columns=bin_cols),
                pd.DataFrame(columns=season_cols))
    out1 = np.array(out1, float)   # year, t_c, gR, mR, CUE, tau, Tair, SW, VPD, P, sd(gR, mR, CUE, tau)

    with np.errstate(invalid='ignore'):
        gR_fixed, mR_fixed, CUE_fixed, tau_fixed = np.nanmean(out1[:, 2:6], axis=0)
        gR_sd, mR_sd, CUE_sd, tau_sd = np.nanmean(out1[:, 10:14], axis=0)
    mR0 = np.nanpercentile(out1[:, 3], 10, method='hazen')           # MATLAB prctile
    T0 = out1[int(np.nanargmin(np.abs(mR0 - out1[:, 3]))), 6]

    # ---- round 2: CUE and Ea per 1-degree air-temperature bin
    out2, stuck = [], []           # rows: year, t_c, gR, CUE, Ea, tau, sd(gR, CUE, Ea, tau)
    for yr in years:
        gpp, reco, cumgpp, cli, doy = by_year[yr]
        tair, rain = cli[:, 0], cli[:, 3]
        for t_c in range(0, 36):
            with np.errstate(invalid='ignore'):
                sel = np.flatnonzero((tair > t_c) & (tair < t_c + 1) & (rain < 2))
            if len(sel) <= 5:
                continue
            dy, dx1, dx2 = pair_diffs(reco[sel]), pair_diffs(gpp[sel]), pair_diffs(cumgpp[sel])
            arrh_arg = 1.0 / (t_c + 0.5 + 273.15) - 1.0 / (T0 + 273.15)

            def ss2(kb):
                model = gR_fixed * kb[0] * dx1 + np.exp((-kb[1] / BOLTZMANN_EV) * arrh_arg) \
                    * mR0 * kb[0] * dx2 * (1 - tau_fixed)
                return np.sum((dy - model) ** 2) * 0.05 + np.sum((dy - (1 - kb[0]) * dx1) ** 2 * 0.95)

            mean, chain = dram(ss2, [CUE_fixed, 0.65], [0.1, 0.0], [1.0, 1.2], rng)
            moved = abs(mean[0] - CUE_fixed) > 0.001 and abs(mean[1] - 0.65) > 0.01
            if moved or not out2:
                cue_ea, cue_ea_sd = mean, chain.std(axis=0)
            else:                  # chain did not move: carry the previous bin forward
                cue_ea = out2[-1][3:5]
                cue_ea_sd = out2[-1][3:5] if mode == 'matlab' else out2[-1][7:9]
            out2.append([yr, t_c, gR_fixed, cue_ea[0], cue_ea[1], tau_fixed, gR_sd, cue_ea_sd[0], cue_ea_sd[1], tau_sd])
            stuck.append(not moved)
    out2, stuck = np.array(out2, float).reshape(-1, 10), np.array(stuck, bool)

    # ---- annual values: bootstrap over the year's bins
    rows = []
    for yr in years:
        y_ind = out2[:, 0] == yr
        row = {'year': yr, 'n_temperature_bins': int(y_ind.sum()), 'n_bins_not_converged': int(stuck[y_ind].sum())}
        if y_ind.any():
            mu, sd = out2[y_ind, 2:6], out2[y_ind, 6:10]
            sd = np.where(sd < 0, np.nan, sd)
            draws = rng.normal(mu, sd, size=(N_BOOT, *mu.shape))
            with np.errstate(invalid='ignore'):
                for p, name in enumerate(PARAMS):
                    row[name] = np.nanmean(draws[:, :, p])
                    row[f'{name}_sd'] = np.nanstd(draws[:, :, p])
        rows.append(row)
    annual = pd.DataFrame(rows).reindex(columns=cols)
    # site constants of the respiration model (the same for every year of the site)
    annual['mR0'], annual['T0'], annual['tau_site'], annual['gR_site'] = mR0, T0, tau_fixed, gR_fixed
    bins =pd.DataFrame({'year': out2[:, 0].astype(int), 'tair_bin': out2[:, 1] + 0.5, 'CUE': out2[:, 3],
                         'Ea': out2[:, 4], 'CUE_sd': out2[:, 7], 'Ea_sd': out2[:, 8], 'converged': ~stuck})

    # ---- round 3: seasonal CUE in sliding time windows, pairs of similar air temperature
    season = []
    half = SEASON_WINDOW_DAYS / 2.0
    ann = annual.set_index('year')
    for yr in years:
        gpp, reco, cumgpp, cli, doy = by_year[yr]
        tair, rain = cli[:, 0], cli[:, 3]
        cue0, ea_y = ann.at[yr, 'CUE'], ann.at[yr, 'Ea']
        cue0 = float(np.clip(cue0 if np.isfinite(cue0) else CUE_fixed, *CUE_BOUNDS))
        ea_y = ea_y if np.isfinite(ea_y) else 0.65
        with np.errstate(invalid='ignore'):
            valid = np.isfinite(gpp) & np.isfinite(reco) & np.isfinite(tair) & (rain < 2)
        for c in np.arange(half, 366 - half + 1e-9, SEASON_STEP_DAYS):
            sel = np.flatnonzero(valid & (doy >= c - half) & (doy < c + half))
            if len(sel) < SEASON_MIN_DAYS:
                continue
            i, j = np.triu_indices(len(sel), k=1)
            keep = np.abs(tair[sel][j] - tair[sel][i]) <= SEASON_MAX_DT
            if keep.sum() < SEASON_MIN_PAIRS:
                continue
            i, j = sel[i[keep]], sel[j[keep]]
            dy, dx1, dx2 = reco[j] - reco[i], gpp[j] - gpp[i], cumgpp[j] - cumgpp[i]
            arrh = np.exp((-ea_y / BOLTZMANN_EV) * (1.0 / (0.5 * (tair[i] + tair[j]) + 273.15) - 1.0 / (T0 + 273.15)))
            maint = arrh * mR0 * dx2 * (1 - tau_fixed)

            def ss3(k):
                model = gR_fixed * k[0] * dx1 + k[0] * maint
                return np.sum((dy - model) ** 2) * 0.05 + np.sum((dy - (1 - k[0]) * dx1) ** 2 * 0.95)

            mean, chain = dram(ss3, [cue0], [CUE_BOUNDS[0]], [CUE_BOUNDS[1]], rng)
            season.append([yr, c, mean[0], chain.std(), len(sel), int(keep.sum()), float(np.mean(tair[sel])),
                           abs(mean[0] - cue0) > 0.001])
    seasonal = pd.DataFrame(season, columns=season_cols)
    return annual, bins, seasonal


def daily_cue_from_windows(seasonal):
    """Window CUE -> one value per day of year, per year. Returns year, doy, CUE_daily."""
    parts = []
    for yr, g in seasonal[seasonal['converged'].astype(bool)].groupby('year'):
        g = g.sort_values('doy_center')
        if len(g) < SEASON_MIN_WINDOWS:
            continue
        smooth = g['CUE'].rolling(3, center=True, min_periods=1).median().to_numpy()
        days = np.arange(1, 367)
        parts.append(pd.DataFrame({'year': int(yr), 'doy': days,
                                   'CUE_daily': np.interp(days, g['doy_center'].to_numpy(), smooth)}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=['year', 'doy', 'CUE_daily'])


# ------------------------------------------------------------------ readfluxnet_data
def read_fluxnet_dd(source):
    """Daily FLUXNET2015 / ONEFlux table -> the columns estimate_site_cue needs.
    `source`: path to a DD csv, or to a FLUXNET zip containing one."""
    source = str(source)
    if source.lower().endswith('.zip'):
        with zipfile.ZipFile(source) as z:
            dd = [n for n in z.namelist() if re.search(r'_DD_', n) and n.endswith('.csv') and 'ERA5' not in n
                  and 'BIF' not in n]
            if not dd:
                return None, None
            raw_bytes = io.BytesIO(z.read(sorted(dd, key=lambda n: 'FULLSET' not in n.upper())[0]))
        raw = pd.read_csv(raw_bytes, usecols=lambda c: c in NEEDED)
    else:
        raw = pd.read_csv(source, usecols=lambda c: c in NEEDED)
    raw = raw.reindex(columns=NEEDED).apply(pd.to_numeric, errors='coerce')
    raw = raw.where(raw >= -1000)                                    # -9999 -> NaN
    date = pd.to_datetime(raw['TIMESTAMP'].astype('Int64').astype(str), format='%Y%m%d', errors='coerce')
    d = pd.DataFrame({'date': date, 'year': date.dt.year, 'doy': date.dt.dayofyear,
                      'TA': raw['TA_F'], 'SW': raw['SW_IN_F'], 'VPD': raw['VPD_F'] / 10.0, 'P': raw['P_F'],
                      'GPP_NT': raw['GPP_NT_VUT_REF'], 'GPP_DT': raw['GPP_DT_VUT_REF'],
                      'RECO_NT': raw['RECO_NT_VUT_REF'], 'RECO_DT': raw['RECO_DT_VUT_REF']}).dropna(subset=['date'])
    has_nt = d['GPP_NT'].notna().any() and d['RECO_NT'].notna().any()
    has_dt = d['GPP_DT'].notna().any() and d['RECO_DT'].notna().any()
    if has_nt and has_dt:
        partition = 'NT+DT'
    elif has_nt:
        partition, d['GPP_DT'], d['RECO_DT'] = 'NT', d['GPP_NT'], d['RECO_NT']
    elif has_dt:
        partition, d['GPP_NT'], d['RECO_NT'] = 'DT', d['GPP_DT'], d['RECO_DT']
    else:
        return None, None
    d['year'] = d['year'].astype(int)
    return d, partition


def cache_path(site_id):
    return CACHE_DIR / f"{site_id}.pkl"


def cache_key(source):
    """Changes when the input archive or any setting of the method changes."""
    st = os.stat(source)
    return (os.path.basename(str(source)), st.st_size, int(st.st_mtime), INDEXING_MODE, NSIMU, ADAPT_INT, N_BOOT,
            SEASON_WINDOW_DAYS, SEASON_STEP_DAYS, SEASON_MAX_DT, SEASON_MIN_DAYS, SEASON_MIN_PAIRS, CACHE_VERSION)


def load_cached(site_id, source):
    p = cache_path(site_id)
    if RESTART or not p.exists():
        return None
    try:
        saved = pd.read_pickle(p)
    except Exception:
        return None
    return saved['tables'] if saved.get('key') == cache_key(source) else None


def run_site(args):
    """One site, saved to the cache as soon as it is done - an interrupted run loses only the sites in progress."""
    site_id, source = args
    site_id, out, err = compute_site(site_id, source)
    if out is not None:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        tmp = cache_path(site_id).with_suffix('.tmp')
        pd.to_pickle({'key': cache_key(source), 'tables': out}, tmp)
        os.replace(tmp, cache_path(site_id))
    return site_id, out, err


def compute_site(site_id, source):
    try:
        d, partition = read_fluxnet_dd(source)
    except zipfile.BadZipFile:
        return site_id, None, 'corrupt zip (re-download with step 11)'
    if d is None:
        return site_id, None, 'no daily GPP / Reco'
    rng = np.random.default_rng(zlib.crc32(site_id.encode()))        # reproducible per site
    annual, bins, seasonal = estimate_site_cue(d, rng)
    daily = daily_cue_from_windows(seasonal)
    for t in (annual, bins, seasonal, daily):
        t.insert(0, 'site_id', site_id)
    annual['flux_partition'], annual['indexing_mode'] = partition, INDEXING_MODE
    return site_id, (annual, bins, seasonal, daily), None


if __name__ == '__main__':
    files = sys.argv[1:]
    if files:                       # explicit files, e.g. the US-Ha1 sample shipped with the MATLAB code
        jobs = []
        for f in files:
            m = re.search(r'[A-Z]{2}-[A-Za-z0-9]{3}', os.path.basename(f))
            jobs.append((m.group(0) if m else Path(f).stem, f))
    else:
        if not os.path.exists(FLUX_CSV):
            raise FileNotFoundError(f"Missing '{FLUX_CSV}'. Run steps 11-22 first.")
        sites = sorted(pd.read_csv(FLUX_CSV, usecols=['site_id'])['site_id'].unique())
        zips = glob.glob(str(RAW_DIR / "*.zip"))
        jobs = []
        for s in sites:
            match = [z for z in zips if re.search(rf'(^|[_-]){re.escape(s)}([_-]|$)', os.path.basename(z))]
            if match:
                jobs.append((s, match[0]))
    # sites finished in an earlier run (same input file, same settings) are read back, not recomputed
    by_site, todo = {}, []
    for site_id, source in jobs:
        saved = None if files else load_cached(site_id, source)
        if saved is None:
            todo.append((site_id, source))
        else:
            by_site[site_id] = saved
    print(f"Estimating CUE for {len(jobs)} sites (indexing mode: '{INDEXING_MODE}', {N_WORKERS} workers): "
          f"{len(by_site)} already done in '{CACHE_DIR.name}', {len(todo)} to compute "
          f"(set CUE_RESTART=1 to recompute all)...", flush=True)

    skipped = {}
    if todo:
        with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(todo))) as pool:
            pending = [pool.submit(compute_site if files else run_site, *((j,) if not files else j)) for j in todo]
            for k, fut in enumerate(as_completed(pending), start=1):
                site_id, out, err = fut.result()
                if err:
                    skipped[site_id] = err
                    continue
                by_site[site_id] = out
                annual, _, seasonal, daily = out
                ok = annual['CUE'].notna()
                print(f"  [{k}/{len(todo)}] {site_id}: {int(ok.sum())}/{len(annual)} years, "
                      f"mean CUE {annual.loc[ok, 'CUE'].mean():.3f}; seasonal windows "
                      f"{int(seasonal['converged'].sum()) if len(seasonal) else 0}, "
                      f"daily CUE for {daily['year'].nunique()} years - saved", flush=True)
    results = [by_site[s] for s, _ in jobs if s in by_site]          # fixed order, whatever finished first
    if not results:
        raise SystemExit(f"No site produced a CUE estimate. Skipped: {skipped}")
    cue, bins, seasonal, daily = (pd.concat([r[n] for r in results if len(r[n])], ignore_index=True)
                                  for n in range(4))

    if files:                       # sample run: print, do not touch the pipeline outputs
        pd.set_option('display.width', 200)
        print(cue.drop(columns=['indexing_mode']).round(3).to_string(index=False))
        conv = seasonal[seasonal['converged'].astype(bool)]
        print("\nSeasonal CUE (median over years of the converged windows), by window centre:")
        print(conv.groupby('doy_center')['CUE'].agg(['median', 'count']).round(3).T.to_string())
        sys.exit(0)

    cue.to_csv(OUT_CUE, index=False)
    bins.to_csv(OUT_BINS, index=False)
    seasonal.to_csv(OUT_SEASON, index=False)
    print(f"\nSite-year CUE -> '{OUT_CUE}' ({int(cue['CUE'].notna().sum())} site-years with CUE, "
          f"{cue.loc[cue['CUE'].notna(), 'site_id'].nunique()} sites)")
    print(f"CUE by temperature bin -> '{OUT_BINS}' ({len(bins):,} bins)")
    print(f"Seasonal CUE -> '{OUT_SEASON}' ({int(seasonal['converged'].sum()):,} converged windows of "
          f"{len(seasonal):,})")
    if skipped:
        print(f"Skipped {len(skipped)} sites: {skipped}")

    # daily NPP on the pipeline's own GPP (QC-passing site-years), with the annual and the daily CUE
    flux = pd.read_csv(FLUX_CSV, usecols=['site_id', 'TIMESTAMP', 'GPP_NT_VUT_REF', 'TA_F'])
    flux['date'] = pd.to_datetime(flux['TIMESTAMP'].astype(str), format='%Y%m%d', errors='coerce')
    flux = flux.dropna(subset=['date']).sort_values(['site_id', 'date'])
    flux['year'], flux['doy'] = flux['date'].dt.year, flux['date'].dt.dayofyear
    gpp = flux['GPP_NT_VUT_REF'].where(flux['GPP_NT_VUT_REF'].between(-5, 50))      # -9999 = missing
    ta = flux['TA_F'].where(flux['TA_F'].between(-60, 50))
    flux['cum_gpp'] = gpp.clip(lower=0).fillna(0).groupby([flux['site_id'], flux['year']]).cumsum()
    par = cue[['site_id', 'year', 'CUE', 'gR', 'Ea', 'tau', 'mR0', 'T0']].rename(columns={'CUE': 'CUE_annual'})
    npp = flux.assign(gpp=gpp, ta=ta).merge(par, on=['site_id', 'year'], how='inner')
    npp = npp.merge(daily, on=['site_id', 'year', 'doy'], how='left')
    npp['NPP'] = npp['CUE_annual'] * npp['gpp']
    npp['NPPd'] = npp['CUE_daily'] * npp['gpp']
    # respiration terms of the fitted model. Rm does not use the day's GPP: it is the maintenance cost of
    # the biomass built since 1 January (CUE x cumulative GPP, less turnover), at the day's temperature.
    arrh = np.exp((-npp['Ea'] / BOLTZMANN_EV) * (1.0 / (npp['ta'] + 273.15) - 1.0 / (npp['T0'] + 273.15)))
    npp['Rg'] = npp['gR'] * npp['CUE_annual'] * npp['gpp']
    shape = arrh * npp['mR0'] * npp['CUE_annual'] * npp['cum_gpp'] * (1 - npp['tau'])
    # The model is fitted on day-to-day DIFFERENCES, so this expression gives the seasonal course of Rm but
    # not its level (taken literally it exceeds GPP at most sites). It is therefore scaled, per site-year,
    # so that Rg + Rm adds up to the year's autotrophic respiration (1 - CUE) x GPP.
    key = [npp['site_id'], npp['year']]
    ra_year = ((1 - npp['CUE_annual']) * npp['gpp']).groupby(key).transform('sum')
    rm_year = (ra_year - npp['Rg'].groupby(key).transform('sum')).clip(lower=0)
    npp['Rm'] = shape / shape.groupby(key).transform('sum').where(lambda s: s > 0) * rm_year
    npp = npp.dropna(subset=['NPP', 'NPPd'], how='all')[['site_id', 'date', 'CUE_annual', 'CUE_daily', 'NPP', 'NPPd',
                                                         'Rg', 'Rm']]
    npp.to_csv(OUT_NPP, index=False)
    print(f"Daily NPP -> '{OUT_NPP}' ({len(npp):,} days, {npp['site_id'].nunique()} sites; "
          f"NPP on {int(npp['NPP'].notna().sum()):,} days, NPPd on {int(npp['NPPd'].notna().sum()):,})")
