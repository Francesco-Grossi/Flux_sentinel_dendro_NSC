#!/usr/bin/env bash
#
# Runs the FLUXNET / HLS / PhenoCam phenology - carbon pipeline in order.
# Stops immediately if any step fails, so a later script never silently runs
# on stale or missing data from a broken earlier step.
#
# Steps are numbered by stage:
#   1x  raw downloads            (network; slow)
#   2x  filtering and fitting    (QC, merge, phenology, CUE / NPP, predictors)
#   3x-42  analysis              (results, figures, summaries)
#
# Usage:
#   ./run_pipeline.sh              run everything, from scratch
#   ./run_pipeline.sh 21           resume from step 21 to the end
#   ./run_pipeline.sh 21 27        run steps 21 through 27
#   ./run_pipeline.sh 30 42        only the analysis
#   ./run_pipeline.sh --list       show the steps and exit
#
# Logs: each step's stdout+stderr goes to logs/<step>.log AND the console at
# the same time (tee), so a long-running step can be watched live.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="${ROOT}/Script"
LOG_DIR="${ROOT}/logs"
mkdir -p "$LOG_DIR"

STEPS=(
  # ---- 1x raw downloads
  "11_fluxnet_download.py"              # FLUXNET daily files -> data_raw/, data/fluxnet_daily_all_vars.csv
  "12_hls_extraction.py"                # HLS Landsat + Sentinel-2 reflectance (Google Earth Engine)
  "13_phenocam_download.py"             # PhenoCam GCC for the matched flux sites
  # ---- 2x filtering and fitting
  "21_growing_season_qc_filter.py"      # drop site-years with poor growing-season flux quality
  "22_merge_fluxnet_hls.py"             # daily flux + satellite indices
  "23_phenology_satellite.py"           # leaf-out / EOS from NDVI, NIRv
  "24_phenology_tower_ndvi.py"          # leaf-out / EOS from tower broadband NDVI
  "25_phenology_phenocam.py"            # leaf-out / EOS from PhenoCam GCC
  "26_cue_npp_luo2025.py"               # CUE (annual, seasonal, daily) and daily NPP
  "27_window_predictors.py"             # carbon / climate window predictors per site-year
  "28_autumn_phenology_predictors.py"   # satellite EOS + split-GPP windows, for the 1:1 plots
  # ---- 3x analysis
  "30_eos_satellite_vs_gpp.py"          # GPP-derived EOS10; each EOS source against it
  "31_compare_eos_sources.py"           # agreement between EOS sources
  "32_cue_seasonal_summary.py"          # seasonal course of CUE; annual vs seasonal
  "33_eos_window_scan.py"               # which carbon window relates to which EOS
  "34_solstice_sliding_window_scan.py"  # when around the solstice the relation is strongest
  "35_eos_env_vs_carbon_models.py"      # does carbon add to climate in explaining EOS
  "36_eos_path_analysis.py"             # climate -> carbon -> EOS pathways
  "37_eos_rate_vs_cumulative.py"        # rate vs cumulative uptake
  "38_eos_results_by_leaf_habit.py"     # one results table, split by plant type
  "39_anomaly_timing_drought.py"        # timing of anomalies; drought vs non-drought years
  "40_plot_predictor_correlations.py"   # 1:1 scatter of every EOS x predictor pair (split-GPP effect)
  "41_split_gpp_cancellation_test.py"   # formal test: pre negative, post positive, total none
  "42_results_report.py"                # one Markdown report with all results and figures
)

step_number() { echo "${1%%_*}"; }

if [[ "${1:-}" == "--list" || "${1:-}" == "-l" ]]; then
    for s in "${STEPS[@]}"; do echo "  $(step_number "$s")  ${s}"; done
    exit 0
fi

START="${1:-$(step_number "${STEPS[0]}")}"
END="${2:-$(step_number "${STEPS[${#STEPS[@]}-1]}")}"

# Find a real Python interpreter. An activated virtual environment comes first
# on PATH. On Windows, `python3` is often just the Microsoft Store stub (prints
# a Store-redirect message and exits nonzero), so each candidate is checked for
# actually working, not just for existing on PATH.
# The command is kept as an array so that a project path containing spaces
# (and the two-word "py -3") both work.
PYTHON=()
for candidate in "${ROOT}/.venv/Scripts/python" "${ROOT}/.venv/bin/python" python3 python "py -3"; do
    if [[ "$candidate" == "py -3" ]]; then cmd=(py -3); else cmd=("$candidate"); fi
    if "${cmd[@]}" --version >/dev/null 2>&1; then
        PYTHON=("${cmd[@]}")
        break
    fi
done
if [ ${#PYTHON[@]} -eq 0 ]; then
    echo "ERROR: no working Python interpreter found (tried: .venv, python3, python, py -3)."
    echo "Create the environment first - see README.md, section 'Setup'."
    exit 1
fi
echo "Using Python interpreter: ${PYTHON[*]} ($("${PYTHON[@]}" --version 2>&1))"

echo "=================================================================="
echo "Phenology - carbon pipeline: steps ${START} to ${END}"
echo "Logs: ${LOG_DIR}"
echo "=================================================================="

PIPELINE_START=$(date +%s)
ran=0

for script in "${STEPS[@]}"; do
    num=$(step_number "$script")
    if (( 10#$num < 10#$START || 10#$num > 10#$END )); then
        continue
    fi
    log_file="${LOG_DIR}/${script%.py}.log"

    echo
    echo "------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step ${num}: ${script}"
    echo "------------------------------------------------------------------"

    step_start=$(date +%s)
    # -u = unbuffered stdout. Without it, Python block-buffers output when
    # it's piped (as it is here, through tee), so print() statements don't
    # show up in the log/console until the buffer fills or the script exits.
    if "${PYTHON[@]}" -u "${SCRIPT_DIR}/${script}" 2>&1 | tee "${log_file}"; then
        step_end=$(date +%s)
        echo "[OK] ${script} finished in $((step_end - step_start))s (log: ${log_file})"
        ran=$((ran + 1))
    else
        echo
        echo "[FAILED] ${script} - see ${log_file} for details."
        echo "Fix the issue and resume with: $0 ${num}"
        exit 1
    fi
done

PIPELINE_END=$(date +%s)
echo
echo "=================================================================="
echo "${ran} step(s) complete in $((PIPELINE_END - PIPELINE_START))s."
echo "=================================================================="
