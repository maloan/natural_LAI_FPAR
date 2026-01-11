#!/usr/bin/env bash
# ==============================================================================
# make_relative_trends_georef.sh
#
# Purpose
#   Create relative (relative-to-local-mean) trend rasters from existing
#   UNMASKED / GEOREF 0.25° NetCDF products.
#
# Definition (per pixel)
#   relative_trend [yr^-1] = slope_peryear / mean_level
#
# Usage
#   VARS="LAI FPAR" METRICS="yearmean yearmax yearamp" \
#   bash make_relative_trends_georef.sh
# ==============================================================================

set -euo pipefail

ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
DIR025="${ROOT}/analysis/unmasked/0p25"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need_cmd cdo

# Defaults (override via environment if you want)
VARS_ARR=(${VARS:-"LAI FPAR"})
METRICS_ARR=(${METRICS:-"yearmean yearmax yearmin yearamp"})

EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"

cd "$DIR025"

echo "============================================================"
echo "relative trends from existing GEOREF 0.25° files"
echo "DIR  = ${DIR025}"
echo "VARS = ${VARS_ARR[*]}"
echo "METS = ${METRICS_ARR[*]}"
echo "EPS  = LAI:${EPS_LAI}  FPAR:${EPS_FPAR}"
echo "============================================================"

for VAR in "${VARS_ARR[@]}"; do
  if [[ "$VAR" == "LAI" ]]; then
    EPS="$EPS_LAI"
  else
    EPS="$EPS_FPAR"
  fi

  echo
  echo "=== VAR: $VAR (EPS=$EPS) ==="
for MET in "${METRICS_ARR[@]}"; do
  mu="${VAR}_georef_${MET}_0p25.nc"
  sl="${VAR}_georef_${MET}_trend_slope_peryear_0p25.nc"
  out="${VAR}_georef_${MET}_trend_relative_peryear_0p25.nc"
  mu_mean="${VAR}_georef_${MET}_meanlevel_0p25.nc"

  if [[ ! -f "$mu" || ! -f "$sl" ]]; then
    echo "  [skip] ${MET}: missing $( [[ -f "$mu" ]] || echo "$mu" ) $( [[ -f "$sl" ]] || echo "$sl" )"
    continue
  fi

  echo "  meanlevel: timmean(${mu}) -> ${mu_mean}"
  cdo -O timmean "$mu" "$mu_mean"

  echo "  relative: ${sl} / ${mu_mean} -> ${out}"
  cdo -O setname,slope \
    -ifthen -gec,"${EPS}" "$mu_mean" \
    -div "$sl" "$mu_mean" \
    "$out"

  # Optional: remove the intermediate file if you don't need it
rm -f "$mu_mean"
done

done

echo
echo "Done. Wrote:"
echo "  *_georef_*_trend_relative_peryear_0p25.nc"
echo "in:"
echo "  ${DIR025}"
