#!/usr/bin/env bash
# ==============================================================================
# 02_batch_build_trends_masked.sh
# Batch wrapper for build_trends_masked_0p25.sh
# Runs the three CCI thresholds and the single GLC branch used in the paper.
# ==============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_SCRIPT="$SCRIPT_DIR/build_trends_masked_0p25.sh"
[[ -x "$BUILD_SCRIPT" ]] || {
  echo "ERROR: Missing executable script: $BUILD_SCRIPT" >&2
  exit 1
}

VARS_ARR=(LAI FPAR)
RUNS=(
  "alpha_0.05 CCI"
  "alpha_0.1 CCI"
  "alpha_0.2 CCI"
  "alpha_0.1 GLC"
)

total_jobs=0
failed_jobs=0
FAILED_LIST=()

for run in "${RUNS[@]}"; do
  read -r RUN_TAG MASK <<< "$run"
  for VAR in "${VARS_ARR[@]}"; do
    total_jobs=$((total_jobs + 1))
    IN_DIR="$ROOT/output/${RUN_TAG}/masked_0p25/${VAR}/masked_${VAR}_${MASK}"
    if [[ ! -d "$IN_DIR" ]]; then
      echo "WARNING: Missing input directory: $IN_DIR"
      failed_jobs=$((failed_jobs + 1))
      FAILED_LIST+=("$RUN_TAG/$VAR/$MASK")
      continue
    fi
    echo ">>> $RUN_TAG / $VAR / $MASK"
    if bash "$BUILD_SCRIPT" "$RUN_TAG" "$VAR" "$MASK"; then
      echo ">>> SUCCESS"
    else
      echo ">>> FAILED"
      failed_jobs=$((failed_jobs + 1))
      FAILED_LIST+=("$RUN_TAG/$VAR/$MASK")
    fi
    echo ""
  done
done

echo "============================================================"
echo "Total jobs:  $total_jobs"
echo "Failed jobs: $failed_jobs"
if [[ "$failed_jobs" -gt 0 ]]; then
  echo ""
  echo "Failed combinations:"
  printf '  %s\n' "${FAILED_LIST[@]}"
  exit 1
fi
