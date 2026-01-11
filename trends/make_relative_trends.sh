#!/usr/bin/env bash
# ==============================================================================
# make_relative_trends.sh
#
# Purpose
#   Create relative (relative-to-local-mean) trend rasters from existing
#   masked 0.25° products in:
#     output/<tau>/eval/trend_{LAI,FPAR}_{CCI,GLC}/
#
# Definition (per pixel)
#   relative_trend [yr^-1] = slope_peryear / mean_level
#
# Notes
#   - If the "mean level" file contains a time dimension, we first collapse it
#     to a 2-D climatological mean (timmean) to avoid time-resolved outputs.
# ==============================================================================

set -euo pipefail

TAU="${1:-tau_0.05}"

ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
EVAL_DIR="${ROOT}/output/${TAU}/eval"

cd "${EVAL_DIR}"

echo "============================================================"
echo "Building relative trends in: ${EVAL_DIR}"
echo "TAU = ${TAU}"
echo "============================================================"

# thresholds to avoid division by ~0 (tune if you want)
EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"

METRICS=(yearmean yearmax yearmin yearamp)
DIRS=(trend_FPAR_CCI trend_FPAR_GLC trend_LAI_CCI trend_LAI_GLC)

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need_cmd cdo

for d in "${DIRS[@]}"; do
  [[ -d "$d" ]] || { echo "[WARN] missing dir: $d (skip)"; continue; }
  echo
  echo "=== ${d} ==="
  cd "$d"

  VAR="$(echo "$d" | cut -d_ -f2)"   # LAI or FPAR
  MASK="$(echo "$d" | cut -d_ -f3)"  # CCI or GLC

  if [[ "$VAR" == "LAI" ]]; then
    EPS="$EPS_LAI"
  else
    EPS="$EPS_FPAR"
  fi

  echo "VAR  = ${VAR}"
  echo "MASK = ${MASK}"
  echo "EPS  = ${EPS}"

  for m in "${METRICS[@]}"; do
    sl="${VAR}_${m}_trend_slope_peryear_0p25.nc"
    mu="${VAR}_${m}_0p25.nc"
    mu_mean="${VAR}_${m}_meanlevel_0p25.nc"
    out="${VAR}_${m}_trend_relative_peryear_0p25.nc"

    if [[ ! -f "$sl" || ! -f "$mu" ]]; then
      echo "  [skip] ${m}: missing $( [[ -f "$sl" ]] || echo "$sl" ) $( [[ -f "$mu" ]] || echo "$mu" )"
      continue
    fi

    echo "  meanlevel: timmean(${mu}) -> ${mu_mean}"
    cdo -O timmean "$mu" "$mu_mean"

    echo "  -> ${out}"
    # relative trend = slope / meanlevel, masked where meanlevel < EPS
    cdo -O setname,slope \
      -ifthen -gec,"${EPS}" "$mu_mean" \
      -div "$sl" "$mu_mean" \
      "$out"

    rm -f "$mu_mean"
  done

  cd ..
done

echo
echo "Done. relative trend files written under:"
echo "  ${EVAL_DIR}/trend_*"
