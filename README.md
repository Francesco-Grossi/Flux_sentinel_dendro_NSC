# Phenology - carbon pipeline (FLUXNET, HLS, PhenoCam)

Does carbon uptake during the growing season (source: GPP; sink: NPP, CUE)
shift the end of season (EOS) of vegetation at eddy-covariance sites north of
30 N? The pipeline downloads the data, derives EOS from four independent
sources, estimates carbon use efficiency, and tests the relationships.

## Setup (once per machine)

1. Python 3.10 or newer and git.
2. Create the environment and install the packages:

   ```bash
   python -m venv .venv
   source .venv/Scripts/activate      # Windows Git Bash; Linux/macOS: source .venv/bin/activate
   pip install -r requirements.txt
   ```

3. Google Earth Engine (needed by step 12 only): a Google account registered
   for Earth Engine and a cloud project. Sign in once and name the project:

   ```bash
   earthengine authenticate
   export GEE_PROJECT=<your-gee-project-id>
   ```

4. Nothing else needs an account: FLUXNET data come through `fluxnet-shuttle`
   and PhenoCam data through its public API.

If HTTPS requests fail with `CERTIFICATE_VERIFY_FAILED`, an antivirus or
proxy is re-signing connections; add its root certificate to the bundle
shown by `python -m certifi`, or turn off its HTTPS scanning.

## Running

```bash
./run_pipeline.sh              # everything, from scratch
./run_pipeline.sh --list       # show the steps
./run_pipeline.sh 21           # resume from step 21 to the end
./run_pipeline.sh 30 46        # only the analysis
```

Each step can also be run on its own (`python Script/<step>.py`); it stops
with a message naming the step to run first if an input is missing. Logs go
to `logs/<step>.log`.

Steps 11 and 12 are slow (hours: about 30 GB of FLUXNET archives, and Earth
Engine requests for every site). Both resume: step 11 skips sites already
downloaded and re-fetches unreadable archives; step 12 skips sites already
extracted.

## Steps

Numbered by stage. Every step reads only files written by earlier steps.

### 1x - raw downloads

| Step | Script | Writes |
|---|---|---|
| 11 | `11_fluxnet_download.py` | `data_raw/*.zip`, `data/fluxnet_daily_all_vars.csv` (sites > 30 N, natural land cover, >= 10 consecutive years) |
| 12 | `12_hls_extraction.py` | HLS Landsat + Sentinel-2 reflectance within 1 km of each tower (own vegetation type only): `data/fluxnet_all_highlat_landsat_raw_bands.csv`, `..._indices.csv` |
| 13 | `13_phenocam_download.py` | `data/phenocam_site_matches.csv`, `data_raw/phenocam/*_1day.csv` |
| 14 | `14_dendrometer_download.py` | open dendrometer datasets at flux sites: `data_raw/dendro/` (see `Output/dendrometer_datasets.md`) |

### 2x - filtering and fitting

| Step | Script | Writes |
|---|---|---|
| 21 | `21_growing_season_qc_filter.py` | flux site-years with good growing-season data quality |
| 22 | `22_merge_fluxnet_hls.py` | `data/fluxnet_landsat_merged.csv` (daily flux + satellite indices) |
| 23 | `23_phenology_satellite.py` | leaf-out / EOS from satellite NDVI and NIRv |
| 24 | `24_phenology_tower_ndvi.py` | leaf-out / EOS from tower broadband NDVI (radiation sensors) |
| 25 | `25_phenology_phenocam.py` | leaf-out / EOS from PhenoCam greenness (GCC) |
| 26 | `26_cue_npp_luo2025.py` | CUE (annual, by temperature, seasonal, daily) after Luo et al.; daily NPP |
| 27 | `27_window_predictors.py` | carbon and climate window predictors per site-year, for every EOS source |
| 28 | `28_autumn_phenology_predictors.py` | every EOS source with split-GPP windows (year-anchored and fixed), for the 1:1 plots |
| 29 | `29_dendrometer_growth.py` | stem growth per site-year from dendrometers: amount, rate before / after the solstice, onset and cessation dates |

Steps 23-25 share one fitting routine (`pheno_fit.py`: double logistic,
dormant-season background, climatology gap-fill, QC flags), so EOS90 / EOS50 /
EOS10 mean the same thing for every source.

### 3x - analysis

| Step | Script | Question |
|---|---|---|
| 30 | `30_eos_satellite_vs_gpp.py` | GPP-derived EOS10, and each EOS source against it |
| 31 | `31_compare_eos_sources.py` | Do the EOS sources agree from year to year? |
| 32 | `32_cue_seasonal_summary.py` | Seasonal course of CUE; annual vs seasonal vs daily |
| 33 | `33_eos_window_scan.py` | Which carbon window relates to which EOS, and with what sign? |
| 34 | `34_solstice_sliding_window_scan.py` | When around the solstice is the relation strongest? |
| 35 | `35_eos_env_vs_carbon_models.py` | Does carbon explain EOS beyond climate? |
| 36 | `36_eos_path_analysis.py` | Climate -> carbon -> EOS pathways |
| 37 | `37_eos_rate_vs_cumulative.py` | Uptake rate vs cumulative uptake |
| 38 | `38_eos_results_by_leaf_habit.py` | One results table by plant type, with multiple-testing correction |
| 39 | `39_anomaly_timing_drought.py` | Timing of positive / negative anomalies; drought vs non-drought years |
| 40 | `40_plot_predictor_correlations.py` | 1:1 scatter of every EOS x predictor pair (split-GPP effect) |
| 41 | `41_split_gpp_cancellation_test.py` | Are Zani and Lu both right? Pre-solstice GPP negative, post-solstice positive, whole-season none (within sites) |
| 42 | `42_compare_phenocam_satellite.py` | Do PhenoCam and satellite EOS give the same answer on the same site-years? |
| 43 | `43_presolstice_robustness.py` | Is the pre-solstice effect solid? Leave one site out, other QC thresholds, statistical power |
| 44 | `44_mechanism_tests.py` | Why? Leaf-out date, water balance and sink variables (respiration, sink index) against GPP |
| 45 | `45_growth_vs_senescence.py` | Is measured stem growth coupled to GPP, and does it predict the onset of senescence? (few sites) |
| 46 | `46_results_report.py` | Builds `Output/results_report.md`: all results with their figures |

### Not pipeline steps

- `eos_common.py`, `pheno_fit.py` - shared code imported by the steps.
- `repair_fluxnet_downloads.py` - tests every FLUXNET archive and downloads the unreadable ones again, with retries.
- `Script/legacy/` - earlier hypothesis tests (07, 08, 09, 11, 12). They read
  the table of step 28 and are kept for reference; their post-solstice GPP windows end
  at the same year's EOS, which builds a positive relation in by construction;
  steps 33-38 replace them.

## Folders

| Folder | Content |
|---|---|
| `Script/` | the pipeline |
| `data_raw/` | downloads (FLUXNET archives, PhenoCam files); not in git |
| `data/` | tables written by the steps |
| `figure/` | figures written by the analysis steps |
| `Output/` | `results_report.md` (all results and figures) and per-analysis summaries |
| `logs/` | one log per step |

## Method notes

- EOS90 / EOS50 / EOS10: day of year when the fitted greenness curve has
  fallen to 90 / 50 / 10 % of its seasonal amplitude after the peak.
- Carbon windows end at the site's mean EOS ("fixed anchors"), not at the
  same year's EOS, to avoid the window-length effect described above.
- Effects are reported as days of EOS shift per +1 within-site SD of the
  predictor, from within-site models (year minus site mean, standard errors
  clustered by site). Mixed models with a site random intercept also pick up
  differences between sites and give larger effects; step 33 keeps them for
  comparison only.
- CUE: Python port of the MATLAB code of Luo et al.
  (https://www.researchsquare.com/article/rs-3989566/v1); see the docstring
  of step 26 for the two points where it departs from the MATLAB code.
