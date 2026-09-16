"""
PIPELINE STEP 11 - Observed vs. predicted EOS90, broken out by model.

Two figures are produced for each of the two model sets already fit in
steps 08/09:

  1. IN-SAMPLE fit: the model sees all the data it's evaluated on (fixed
     effects + each site's own fitted random intercept). This shows how
     well each model's *shape* matches the data, but will always look
     better than real predictive skill because of the random intercept.
  2. OUT-OF-SAMPLE (leave-one-site-out CV): the same models, but each
     point is predicted from a fit that never saw that site, using fixed
     effects only (no access to that site's random intercept). This is
     the honest version - compare the two figures side by side and the
     gap between them IS the overfitting/site-memorization the random
     intercept was providing.

Two model sets, one figure pair each:
  A) baseline_no_gpp / total_gpp / split_gpp_pre_post          (step 08 models)
  B) A_environment_only / B_environment_plus_leafout /
     C_environment_plus_splitgpp / D_environment_plus_leafout_plus_splitgpp  (step 09 models)

Each panel is annotated with RMSE (days) and R^2, and shows the 1:1 line
(perfect prediction) in addition to the point cloud.

Inputs: data/phenology_flux_predictors_by_site_year_index.csv  (step 6)
        data/rate_of_change_predictors_by_site_year_index.csv  (step 7)

Output: figure/eos90_obs_vs_pred/<model_set>_<in_sample|loso>.png
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIGURE_DIR = Path(__file__).resolve().parent.parent / "figure" / "eos90_obs_vs_pred"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

PREDICTORS_CSV = DATA_DIR / "phenology_flux_predictors_by_site_year_index.csv"  # step 6
RATE_CSV = DATA_DIR / "rate_of_change_predictors_by_site_year_index.csv"        # step 7

VI_INDICES = ['NDVI', 'NIRv']
MIN_OBS, MIN_SITES = 20, 5
TARGET = 'EOS90'

BASE_PREDICTORS = ['leaf_out_90', 'temperature_senescence_rate', 'photoperiod_senescence_rate']
MODEL_SET_A = {
    'baseline_no_gpp': BASE_PREDICTORS,
    'total_gpp': BASE_PREDICTORS + ['gpp_total_raw'],
    'split_gpp_pre_post': BASE_PREDICTORS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90'],
}
ENV_VARS = ['temperature_senescence_rate', 'photoperiod_senescence_rate',
            'mean_radiation_growing_season', 'mean_temperature_growing_season',
            'total_reco_growing_season']
MODEL_SET_B = {
    'A_environment_only': ENV_VARS,
    'B_environment_plus_leafout': ENV_VARS + ['leaf_out_90'],
    'C_environment_plus_splitgpp': ENV_VARS + ['gpp_sos10_to_solstice', 'gpp_solstice_to_eos90'],
    'D_environment_plus_leafout_plus_splitgpp': ENV_VARS + ['leaf_out_90', 'gpp_sos10_to_solstice',
                                                              'gpp_solstice_to_eos90'],
}

for p in (PREDICTORS_CSV, RATE_CSV):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing '{p}'. Run steps 06 and 07 first.")


def zscore(s):
    return (s - s.mean()) / s.std()


pred = pd.read_csv(PREDICTORS_CSV)
rate = pd.read_csv(RATE_CSV)[['site_id', 'year', 'vi_index',
                               'temperature_senescence_rate', 'photoperiod_senescence_rate']]
merged = pred.merge(rate, on=['site_id', 'year', 'vi_index'], how='inner')
merged['gpp_total_raw'] = merged['gpp_sos10_to_solstice'] + merged['gpp_solstice_to_eos90']


def get_common_data(vi, model_set):
    """A single subset per VI index, with NaN dropped on the UNION of every
    variable used across the whole model set - so every model in the
    comparison is evaluated on exactly the same site-years (matches the
    consistent-sample approach in steps 08/09; without this, a smaller
    model with fewer required columns would be scored on an easier,
    larger sample than a bigger model and the RMSE/R^2 comparison would
    not be apples-to-apples)."""
    all_vars = sorted(set(v for raw_vars in model_set.values() for v in raw_vars))
    needed = [TARGET, 'site_id'] + all_vars
    return merged[merged['vi_index'] == vi][needed].dropna().reset_index(drop=True)


def in_sample_fit(d, raw_vars):
    dz = d.copy()
    for v in raw_vars:
        dz[v + '_z'] = zscore(dz[v])
    formula = f"{TARGET} ~ " + " + ".join(v + '_z' for v in raw_vars)
    fit = smf.mixedlm(formula, dz, groups=dz['site_id']).fit(reml=True)
    return dz[TARGET].to_numpy(), fit.fittedvalues.to_numpy()


def loso_predict(d, raw_vars):
    """Leave-one-site-out, fixed-effects-only prediction (see step 08/09)."""
    sites = d['site_id'].unique()
    yhat = np.full(len(d), np.nan)
    for held_out_site in sites:
        train_mask = d['site_id'] != held_out_site
        test_mask = ~train_mask
        if test_mask.sum() == 0 or train_mask.sum() < MIN_OBS:
            continue
        train, test = d.loc[train_mask].copy(), d.loc[test_mask].copy()
        means = {v: train[v].mean() for v in raw_vars}
        stds = {v: train[v].std() for v in raw_vars}
        for v in raw_vars:
            train[v + '_zc'] = (train[v] - means[v]) / stds[v]
            test[v + '_zc'] = (test[v] - means[v]) / stds[v]
        formula = f"{TARGET} ~ " + " + ".join(v + '_zc' for v in raw_vars)
        try:
            fit = smf.mixedlm(formula, train, groups=train['site_id']).fit(reml=True)
        except Exception:
            continue
        coef_names = ['Intercept'] + [v + '_zc' for v in raw_vars]
        X_test = np.column_stack([np.ones(len(test))] + [test[v + '_zc'].to_numpy() for v in raw_vars])
        coefs = np.array([fit.params.get(c, 0.0) for c in coef_names])
        yhat[test.index.to_numpy()] = X_test @ coefs
    return d[TARGET].to_numpy(), yhat


def rmse_r2(y, yhat):
    mask = ~np.isnan(yhat)
    resid = y[mask] - yhat[mask]
    rmse = float(np.sqrt(np.mean(resid ** 2)))
    ss_res, ss_tot = float(np.sum(resid ** 2)), float(np.sum((y[mask] - y[mask].mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return rmse, r2, mask


def plot_grid(model_set, mode, out_name):
    n_models = len(model_set)
    fig, axes = plt.subplots(len(VI_INDICES), n_models, figsize=(4.2 * n_models, 8), squeeze=False)

    for row, vi in enumerate(VI_INDICES):
        d_common = get_common_data(vi, model_set)
        for col, (model_name, raw_vars) in enumerate(model_set.items()):
            ax = axes[row][col]
            d = d_common
            if len(d) < MIN_OBS or d['site_id'].nunique() < MIN_SITES:
                ax.text(0.5, 0.5, "not enough data", ha='center', va='center',
                        fontsize=9, color='gray', transform=ax.transAxes)
                ax.set_title(f"{vi} | {model_name}", fontsize=9)
                continue

            if mode == 'in_sample':
                y, yhat = in_sample_fit(d, raw_vars)
            else:
                y, yhat = loso_predict(d, raw_vars)

            rmse, r2, mask = rmse_r2(y, yhat)
            ax.scatter(yhat[mask], y[mask], alpha=0.6, s=22, edgecolor='white', linewidth=0.4)
            lims = [min(y[mask].min(), yhat[mask].min()) - 5, max(y[mask].max(), yhat[mask].max()) + 5]
            ax.plot(lims, lims, color='#c53030', linewidth=1.5, linestyle='--', label='1:1')
            ax.set_xlim(lims)
            ax.set_ylim(lims)
            ax.set_xlabel("Predicted EOS90 (DOY)")
            ax.set_ylabel("Observed EOS90 (DOY)")
            ax.set_title(f"{vi} | {model_name}\nRMSE={rmse:.1f}d, R2={r2:.2f}, n={int(mask.sum())}", fontsize=9)
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.legend(loc='upper left', fontsize=7)

    mode_label = "in-sample fit" if mode == 'in_sample' else "leave-one-site-out (out-of-sample)"
    fig.suptitle(f"Observed vs. predicted EOS90 - {mode_label}", fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = FIGURE_DIR / out_name
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved '{out_path}'.")


print("Model set A (step 08 models): baseline / total_gpp / split_gpp")
plot_grid(MODEL_SET_A, 'in_sample', 'model_set_A_step08_in_sample.png')
plot_grid(MODEL_SET_A, 'loso', 'model_set_A_step08_loso.png')

print("\nModel set B (step 09 models): environment / +leafout / +splitgpp / all")
plot_grid(MODEL_SET_B, 'in_sample', 'model_set_B_step09_in_sample.png')
plot_grid(MODEL_SET_B, 'loso', 'model_set_B_step09_loso.png')

print(f"\nDone. Figures saved under '{FIGURE_DIR}'.")
