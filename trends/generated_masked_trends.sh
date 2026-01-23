#!/usr/bin/env bash
# ==============================================================================
# build_trends_masked_0p25.sh — Build masked 0.25° monthly NetCDF, annual metrics,
#                               OLS trends, and relative trends (slope/meanlevel)
#
# Purpose
#   From masked monthly GeoTIFFs (e.g., VAR_masked_YYYYMM_0p25.tif):
#     1) GeoTIFF -> dated NetCDF (monthly) and mergetime to one monthly NetCDF
#     2) Annual metrics: yearmean, yearmax, yearmin, yearamp (=max-min)
#     3) OLS trends (cdo trend): intercept + slope_peryear on each annual metric
#     4) Relative trends: slope_peryear / meanlevel, masked where meanlevel < EPS
#     5) Move outputs to ../../../eval/trend_${VAR}_${MASKTAG}
#
# Usage
#   ./build_trends_masked_0p25.sh VAR MASKTAG [RES]
#
# Examples
#   ./build_trends_masked_0p25.sh FPAR GLC
#   ./build_trends_masked_0p25.sh FPAR CCI
#   ./build_trends_masked_0p25.sh LAI  GLC
#   ./build_trends_masked_0p25.sh LAI  CCI 0p25
#
# Notes
#   - Expects to be run inside a directory containing:
#       ${VAR}_masked_*.tif  (e.g., FPAR_masked_198201_0p25.tif)
#   - Safe cleanup: deletes only files produced by this script.
#   - Relative trend definition matches Methods:
#       r = b / meanlevel, with meanlevel = timmean(annual_metric),
#       and r undefined where meanlevel < EPS.
# ==============================================================================
set -euo pipefail

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need_cmd gdal_translate
need_cmd cdo

# ------------------------------------------------------------------------------
# Args
# ------------------------------------------------------------------------------
VAR=${1:?Provide VAR: LAI or FPAR}
MASKTAG=${2:?Provide MASKTAG: e.g. CCI or GLC}
RES=${3:-0p25}

# ------------------------------------------------------------------------------
# Thresholds for relative trends (paper definition)
# ------------------------------------------------------------------------------
EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"

case "$VAR" in
  LAI)  EPS="$EPS_LAI" ;;
  FPAR) EPS="$EPS_FPAR" ;;
  *) echo "VAR must be LAI or FPAR, got: $VAR" >&2; exit 1 ;;
esac

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------
PATTERN="${VAR}_masked_*.tif"          # expects: VAR_masked_YYYYMM_0p25.tif
MONTHLY="${VAR}_masked_monthly_${RES}.nc"

# ------------------------------------------------------------------------------
# Outputs we generate (for safe cleanup + move)
# ------------------------------------------------------------------------------
ANNUAL_FILES=(
  "${VAR}_yearmean_${RES}.nc"
  "${VAR}_yearmax_${RES}.nc"
  "${VAR}_yearmin_${RES}.nc"
  "${VAR}_yearamp_${RES}.nc"
)

TREND_FILES=(
  "${VAR}_yearmean_trend_intercept_${RES}.nc" "${VAR}_yearmean_trend_slope_peryear_${RES}.nc"
  "${VAR}_yearmax_trend_intercept_${RES}.nc"  "${VAR}_yearmax_trend_slope_peryear_${RES}.nc"
  "${VAR}_yearmin_trend_intercept_${RES}.nc"  "${VAR}_yearmin_trend_slope_peryear_${RES}.nc"
  "${VAR}_yearamp_trend_intercept_${RES}.nc"  "${VAR}_yearamp_trend_slope_peryear_${RES}.nc"
)

REL_FILES=(
  "${VAR}_yearmean_trend_relative_peryear_${RES}.nc"
  "${VAR}_yearmax_trend_relative_peryear_${RES}.nc"
  "${VAR}_yearmin_trend_relative_peryear_${RES}.nc"
  "${VAR}_yearamp_trend_relative_peryear_${RES}.nc"
)

# ------------------------------------------------------------------------------
# Safety cleanup (only what we create)
# ------------------------------------------------------------------------------
rm -f "$MONTHLY" "${ANNUAL_FILES[@]}" "${TREND_FILES[@]}" "${REL_FILES[@]}"
rm -rf tmp_nc
mkdir -p tmp_nc

# ------------------------------------------------------------------------------
# 1) GeoTIFF -> dated NetCDF (monthly) + merge time series
# ------------------------------------------------------------------------------
n_tif=$(ls -1 ${PATTERN} 2>/dev/null | wc -l | awk '{print $1}')
if [[ "$n_tif" -eq 0 ]]; then
  echo "No files match: ${PATTERN}" >&2
  exit 1
fi

echo "============================================================"
echo "VAR     = ${VAR}"
echo "MASKTAG = ${MASKTAG}"
echo "RES     = ${RES}"
echo "EPS     = ${EPS}"
echo "N_TIF   = ${n_tif}"
echo "PWD     = $(pwd)"
echo "============================================================"

echo "Step 1: Converting GeoTIFFs -> dated NetCDF (tmp_nc/)..."
for tif in ${PATTERN}; do
  base=$(basename "$tif")
  yyyymm=${base#${VAR}_masked_}
  yyyymm=${yyyymm%_${RES}.tif}   # e.g. 198201

  year=${yyyymm:0:4}
  mon=${yyyymm:4:2}

  raw_nc="tmp_nc/${VAR}_${yyyymm}.nc"
  dated_nc="tmp_nc/${VAR}_${yyyymm}_dated.nc"

  echo "  ${tif} -> ${dated_nc}"
  gdal_translate -of NETCDF "$tif" "$raw_nc" >/dev/null
  cdo -O setdate,${year}-${mon}-01 "$raw_nc" "$dated_nc" >/dev/null
done

echo "Step 2: Merging monthly time series -> ${MONTHLY}"
cdo -O mergetime tmp_nc/${VAR}_*_dated.nc "${MONTHLY}" >/dev/null
rm -rf tmp_nc

# ------------------------------------------------------------------------------
# 2) Annual metrics (mean/max/min/amp)
# ------------------------------------------------------------------------------
echo "Step 3: Computing annual metrics..."
cdo -O yearmean "${MONTHLY}" "${VAR}_yearmean_${RES}.nc" >/dev/null
cdo -O yearmax  "${MONTHLY}" "${VAR}_yearmax_${RES}.nc"  >/dev/null
cdo -O yearmin  "${MONTHLY}" "${VAR}_yearmin_${RES}.nc"  >/dev/null
cdo -O sub "${VAR}_yearmax_${RES}.nc" "${VAR}_yearmin_${RES}.nc" "${VAR}_yearamp_${RES}.nc" >/dev/null

# ------------------------------------------------------------------------------
# 3) OLS trends on annual metrics (slope per year + intercept)
# ------------------------------------------------------------------------------
echo "Step 4: Computing OLS trends (cdo trend)..."
for MET in yearmean yearmax yearmin yearamp; do
  cdo -O trend "${VAR}_${MET}_${RES}.nc" \
    "${VAR}_${MET}_trend_intercept_${RES}.nc" \
    "${VAR}_${MET}_trend_slope_peryear_${RES}.nc" >/dev/null
done

# ------------------------------------------------------------------------------
# 4) Relative trends: r = slope / meanlevel with EPS threshold on meanlevel
#    meanlevel = timmean(annual metric time series)
# ------------------------------------------------------------------------------
echo "Step 5: Computing relative trends (r = slope/meanlevel; meanlevel >= EPS)..."
for MET in yearmean yearmax yearmin yearamp; do
  ts="${VAR}_${MET}_${RES}.nc"
  sl="${VAR}_${MET}_trend_slope_peryear_${RES}.nc"
  mu="tmp_${VAR}_${MET}_meanlevel_${RES}.nc"
  out="${VAR}_${MET}_trend_relative_peryear_${RES}.nc"

  cdo -O timmean "$ts" "$mu" >/dev/null

  # Keep values only where meanlevel >= EPS (else set to missing)
  cdo -O setname,trend_relative \
    -ifthen -gec,"${EPS}" "$mu" \
    -div "$sl" "$mu" \
    "$out" >/dev/null

  rm -f "$mu"
done

# ------------------------------------------------------------------------------
# 5) Move outputs
# ------------------------------------------------------------------------------
OUT_EVAL="../../../eval/trend_${VAR}_${MASKTAG}"
mkdir -p "$OUT_EVAL"

mv "$MONTHLY" "${ANNUAL_FILES[@]}" "${TREND_FILES[@]}" "${REL_FILES[@]}" "$OUT_EVAL/"

echo "============================================================"
echo "Done."
echo "Wrote outputs to: ${OUT_EVAL}"
echo "Files:"
ls -1 "$OUT_EVAL" | sed 's/^/  - /'
echo "============================================================"
