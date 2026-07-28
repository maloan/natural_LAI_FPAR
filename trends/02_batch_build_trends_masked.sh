#!/usr/bin/env bash
# ==============================================================================
# 02_batch_build_trends_masked.sh
# Batch wrapper for build_trends_masked_0p25.sh
# Generates masked trend products for multiple (ALPHA, VAR, MASKTAG) combinations
#
# Usage:
#   ALPHAS="0.05 0.1 0.2" VARS="LAI FPAR" MASKS="CCI GLC" bash ./02_batch_build_trends_masked.sh
#
# Environment variables (optional)
#   ALPHAS:  space-separated list of alpha values (default: "0.05 0.1 0.2")
#   VARS:  space-separated list of variables (default: "LAI FPAR")
#   MASKS: space-separated list of mask sources (default: "CCI GLC")
# ==============================================================================

set -euo pipefail
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
SCRIPT_DIR="$ROOT/trends"
BUILD_SCRIPT="$SCRIPT_DIR/build_trends_masked_0p25.sh"
[[ -x "$BUILD_SCRIPT" ]] || {
  echo "ERROR: Missing executable script: $BUILD_SCRIPT" >&2
  exit 1
}

read -r -a ALPHAS_ARR <<< "${ALPHAS:-0.05 0.1 0.2}"
read -r -a VARS_ARR <<< "${VARS:-LAI FPAR}"
read -r -a MASKS_ARR <<< "${MASKS:-CCI GLC}"

total_jobs=0
failed_jobs=0
FAILED_LIST=()

for ALPHA in "${ALPHAS_ARR[@]}"; do
  [[ "$ALPHA" =~ ^[0-9]*\.?[0-9]+$ ]] || {
    echo "WARNING: Invalid ALPHA value: $ALPHA"
    continue
  }
  RUN_TAG="alpha_${ALPHA}"
  for VAR in "${VARS_ARR[@]}"; do
    for MASK in "${MASKS_ARR[@]}"; do
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
