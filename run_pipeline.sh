#!/usr/bin/env bash
#
# Runs the full FLUXNET/Landsat phenology-GPP pipeline (Script/01-09) in
# order. Stops immediately if any step fails, so a downstream script never
# silently runs on stale/missing data from a broken upstream step.
#
# Usage:
#   ./run_pipeline.sh              # run all 9 steps
#   ./run_pipeline.sh 5            # resume from step 5 (skip 01-04)
#   ./run_pipeline.sh 5 7          # run only steps 5 through 7
#
# Logs: each step's stdout+stderr goes to logs/<NN>_<name>.log AND the
# console at the same time (tee), so a long-running step (01, 02) can be
# monitored live without losing the saved log.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/Script"
LOG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/logs"
mkdir -p "$LOG_DIR"

# Find a real Python interpreter. On Windows, `python3` is often just the
# Microsoft Store stub (prints a Store-redirect message and exits nonzero)
# even when a real Python is installed as `python` - so python3 is checked
# for actually working (--version succeeds), not just for existing on PATH.
PYTHON=""
for candidate in python3 python "py -3"; do
    if $candidate --version >/dev/null 2>&1; then
        PYTHON="$candidate"
        break
    fi
done
if [ -z "$PYTHON" ]; then
    echo "ERROR: no working Python interpreter found (tried: python3, python, py -3)."
    echo "Install Python from https://www.python.org/downloads/ (check 'Add to PATH'"
    echo "during install), or if it's already installed, open Windows Settings ->"
    echo "Apps -> Advanced app settings -> App execution aliases, and turn OFF the"
    echo "'python.exe'/'python3.exe' entries (those are the Store stub intercepting the command)."
    exit 1
fi
echo "Using Python interpreter: $PYTHON ($($PYTHON --version 2>&1))"

STEPS=(
  "01_fluxnet_download.py"
  "02_landsat_hls_extraction.py"
  "03_growing_season_qc_filter.py"
  "04_merge_fluxnet_landsat.py"
  "05_double_logistic_phenology.py"
  "06_autumn_phenology_predictors.py"
  "07_rate_of_change_predictors.py"
  "08_eos90_split_gpp_hypothesis_test.py"
  "09_eos90_four_model_comparison.py"
  "10_plot_predictor_correlations.py"
  "11_plot_eos90_obs_vs_pred.py"
)

START="${1:-1}"
END="${2:-9}"

echo "=================================================================="
echo "FLUXNET / Landsat phenology-GPP pipeline"
echo "Steps ${START} to ${END} of ${#STEPS[@]}"
echo "Logs: ${LOG_DIR}"
echo "=================================================================="

PIPELINE_START=$(date +%s)

for i in "${!STEPS[@]}"; do
    step_num=$((i + 1))
    if (( step_num < START || step_num > END )); then
        continue
    fi

    script="${STEPS[$i]}"
    log_file="${LOG_DIR}/${script%.py}.log"

    echo
    echo "------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step ${step_num}/${#STEPS[@]}: ${script}"
    echo "------------------------------------------------------------------"

    step_start=$(date +%s)
    # -u = unbuffered stdout. Without it, Python block-buffers output when
    # it's piped (as it is here, through tee) instead of going to a real
    # terminal, so print() statements sit in a buffer and don't show up in
    # the log/console until the buffer fills or the script exits - making a
    # long-running step (like 01/02) look silent even while it's working.
    if $PYTHON -u "${SCRIPT_DIR}/${script}" 2>&1 | tee "${log_file}"; then
        step_end=$(date +%s)
        echo "[OK] ${script} finished in $((step_end - step_start))s (log: ${log_file})"
    else
        echo
        echo "[FAILED] ${script} - see ${log_file} for details."
        echo "Fix the issue and resume with: $0 ${step_num}"
        exit 1
    fi
done

PIPELINE_END=$(date +%s)
echo
echo "=================================================================="
echo "Pipeline complete in $((PIPELINE_END - PIPELINE_START))s."
echo "=================================================================="