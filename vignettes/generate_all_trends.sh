#!/bin/bash
# =============================================================================
# generate_all_trends.sh — Generate trend files for all tau/variable/mask combinations
# =============================================================================
# Prerequisite: run_all_taus.sh must have completed successfully
# This script calls trends/build_trends_masked_0p25.sh for each combination
#
# Usage: bash generate_all_trends.sh
# =============================================================================

set -e  # Exit on first error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SNU_LAI_FPAR_ROOT="$SCRIPT_DIR"

TAUS=("0.05" "0.1" "0.2")
VARS=("LAI" "FPAR")
MASKS=("CCI" "GLC")
TIMESTAMP=$(date "+%Y-%m-%d_%H-%M-%S")
LOG_DIR="$SCRIPT_DIR/logs"

mkdir -p "$LOG_DIR"

echo "==============================================================="
echo "Trend generation started: $TIMESTAMP"
echo "Root: $SNU_LAI_FPAR_ROOT"
echo "Log directory: $LOG_DIR"
echo "==============================================================="

# Check if trend script exists
if [[ ! -f "$SCRIPT_DIR/trends/build_trends_masked_0p25.sh" ]]; then
  echo "ERROR: Trend generation script not found: $SCRIPT_DIR/trends/build_trends_masked_0p25.sh"
  exit 1
fi

total_jobs=0
failed_jobs=0

for TAU in "${TAUS[@]}"; do
  RUN_TAG="tau_${TAU}"

  for VAR in "${VARS[@]}"; do
    for MASK in "${MASKS[@]}"; do
      total_jobs=$((total_jobs + 1))
      LOG_FILE="$LOG_DIR/trend_${RUN_TAG}_${VAR}_${MASK}_${TIMESTAMP}.log"

      echo ""
      echo ">>> Generating trends for $RUN_TAG / $VAR / $MASK"
      echo ">>> Log: $LOG_FILE"
      echo ">>> Started: $(date '+%Y-%m-%d %H:%M:%S')"

      # Check if input directory exists
      IN_DIR="$SCRIPT_DIR/output/${RUN_TAG}/masked_0p25/${VAR}/masked_${VAR}_${MASK}"
      if [[ ! -d "$IN_DIR" ]]; then
        echo "WARNING: Input directory not found: $IN_DIR"
        echo "Skipping $RUN_TAG / $VAR / $MASK"
        failed_jobs=$((failed_jobs + 1))
        continue
      fi

      if bash "$SCRIPT_DIR/trends/build_trends_masked_0p25.sh" "$RUN_TAG" "$VAR" "$MASK" 2>&1 | tee "$LOG_FILE"; then
        echo ">>> Completed: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ">>> Status: SUCCESS"
      else
        echo ">>> Completed: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ">>> Status: FAILED"
        echo "ERROR: Trend generation failed for $RUN_TAG / $VAR / $MASK"
        failed_jobs=$((failed_jobs + 1))
      fi
    done
  done
done

echo ""
echo "==============================================================="
echo "Trend generation completed!"
echo "Total jobs: $total_jobs"
echo "Failed jobs: $failed_jobs"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "==============================================================="

if [[ $failed_jobs -gt 0 ]]; then
  echo "WARNING: $failed_jobs job(s) failed"
  exit 1
fi

exit 0
