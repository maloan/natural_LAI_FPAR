#!/usr/bin/env bash
# ==============================================================================
# make_relative_trends_georef.sh
#   Create relative (baseline-normalized) trend rasters from existing
#   UNMASKED / GEOREF 0.25° NetCDF products.
#
# Usage
#   VARS="LAI FPAR" METRICS="yearmean yearmax yearamp" bash make_relative_trends_georef.sh
# ==============================================================================

set -euo pipefail

ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
DIR025="${ROOT}/analysis/unmasked/0p25"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need_cmd cdo

# Defaults
VARS_ARR=(${VARS:-"LAI FPAR"})
METRICS_ARR=(${METRICS:-"yearmean yearmax yearmin yearamp"})

EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"

cd "$DIR025"

echo "============================================================"
echo "Relative trends from GEOREF 0.25° files"
echo "DIR   = ${DIR025}"
echo "VARS  = ${VARS_ARR[*]}"
echo "METS  = ${METRICS_ARR[*]}"
echo "EPS   = LAI:${EPS_LAI}  FPAR:${EPS_FPAR}"
echo "Rule  = reltrend = slope_peryear / meanlevel ; only if meanlevel >= EPS"
echo "============================================================"

for VAR in "${VARS_ARR[@]}"; do
  case "$VAR" in
    LAI)  EPS="$EPS_LAI" ;;
    FPAR) EPS="$EPS_FPAR" ;;
    *)    echo "Unknown VAR: $VAR (expected LAI/FPAR)"; exit 1 ;;
  esac

  echo
  echo "=== VAR: $VAR (EPS=$EPS) ==="

  for MET in "${METRICS_ARR[@]}"; do
    ts="${VAR}_georef_${MET}_0p25.nc"
    sl="${VAR}_georef_${MET}_trend_slope_peryear_0p25.nc"

    meanlevel="${VAR}_georef_${MET}_meanlevel_0p25.nc"
    mask="${VAR}_georef_${MET}_meanlevel_ge${EPS}_mask_0p25.nc"
    out="${VAR}_georef_${MET}_trend_relative_peryear_0p25.nc"

    if [[ ! -f "$ts" || ! -f "$sl" ]]; then
      echo "  [skip] ${MET}: missing file(s): $( [[ -f "$ts" ]] || echo "$ts" ) $( [[ -f "$sl" ]] || echo "$sl" )"
      continue
    fi

    echo "  meanlevel: timmean(${ts}) -> ${meanlevel}"
    cdo -O timmean "$ts" "$meanlevel"

    echo "  mask: meanlevel >= ${EPS} -> ${mask}"
    cdo -O gec,"${EPS}" "$meanlevel" "$mask"

    echo "  reltrend: (slope / meanlevel) masked -> ${out}"
    # ratio first, then mask; set output variable name to reflect meaning
    cdo -O setname,trend_relative \
      -ifthen "$mask" \
      -div "$sl" "$meanlevel" \
      "$out"
  done
done

