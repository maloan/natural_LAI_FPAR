#!/bin/bash
# =============================================================================
# run_all_taus.sh — Run pipeline for all tau thresholds sequentially
# =============================================================================
# Runs the complete pipeline for tau_0.05, tau_0.1, and tau_0.2
# Usage: bash run_all_taus.sh
# =============================================================================

set -e  # Exit on first error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SNU_LAI_FPAR_ROOT="$SCRIPT_DIR"
cd "$SCRIPT_DIR/R"

TAUS=("0.05" "0.1" "0.2")
TIMESTAMP=$(date "+%Y-%m-%d_%H-%M-%S")
LOG_DIR="$SCRIPT_DIR/logs"

mkdir -p "$LOG_DIR"

echo "==============================================================="
echo "Pipeline run started: $TIMESTAMP"
echo "Targets: tau_${TAUS[@]}"
echo "Log directory: $LOG_DIR"
echo "==============================================================="

echo ""

for TAU in "${TAUS[@]}"; do
  RUN_TAG="tau_${TAU}"
  SETUP_LOG="$LOG_DIR/setup_${RUN_TAG}_${TIMESTAMP}.log"
  PIPELINE_LOG="$LOG_DIR/pipeline_${RUN_TAG}_${TIMESTAMP}.log"
  ANALYSIS_LOG="$LOG_DIR/analysis_${RUN_TAG}_${TIMESTAMP}.log"

  echo ">>> Running setup (script 00) for $RUN_TAG..."
  if make setup RUN_TAG="$RUN_TAG" TAU_CCI="$TAU" 2>&1 | tee "$SETUP_LOG"; then
    echo ">>> Setup completed: $(date '+%Y-%m-%d %H:%M:%S')"
  else
    echo ">>> Setup failed for $RUN_TAG!"
    exit 1
  fi

  echo ""
  echo ">>> Running pipeline for $RUN_TAG (log: $PIPELINE_LOG)"
  echo ">>> Started: $(date '+%Y-%m-%d %H:%M:%S')"

  if make pipeline RUN_TAG="$RUN_TAG" TAU_CCI="$TAU" 2>&1 | tee "$PIPELINE_LOG"; then
    echo ">>> Pipeline completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ">>> Status: SUCCESS"
  else
    echo ">>> Pipeline completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ">>> Status: FAILED"
    echo "ERROR: Pipeline failed for $RUN_TAG"
    exit 1
  fi

  echo ""
  echo ">>> Running analysis for $RUN_TAG (log: $ANALYSIS_LOG)"
  echo ">>> Started: $(date '+%Y-%m-%d %H:%M:%S')"

  if make analysis RUN_TAG="$RUN_TAG" TAU_CCI="$TAU" 2>&1 | tee "$ANALYSIS_LOG"; then
    echo ">>> Analysis completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ">>> Status: SUCCESS"
  else
    echo ">>> Analysis completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ">>> Status: FAILED"
    echo "WARNING: Analysis failed for $RUN_TAG (continuing)"
  fi

done

echo ""
echo "==============================================================="
echo "All pipelines and analyses completed!"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "==============================================================="
