#!/usr/bin/env bash
# ==============================================================================
# build_trends_masked_0p25.sh — Run masked 0.25° trend workflow from repo root
#
# Examples
#   ./build_trends_masked_0p25.sh alpha_0.2 FPAR GLC
#   ./build_trends_masked_0p25.sh alpha_0.1 LAI  CCI
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
cleanup() {
  rm -rf tmp_nc
  rm -f tmp_*.nc
}
trap cleanup EXIT
trap 'echo "ERROR at line $LINENO" >&2; exit 1' ERR

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Missing dependency: $1" >&2
    exit 1
  }
}
need_cmd gdal_translate
need_cmd cdo
need_cmd Rscript
if [[ $# -lt 3 ]]; then
  cat >&2 << 'EOF'
Usage: ./build_trends_masked_0p25.sh ALPHA VAR MASKTAG
Example:
  ./build_trends_masked_0p25.sh alpha_0.1 LAI CCI
EOF
  exit 1
fi
ALPHA="$1"
VAR="$2"
MASKTAG="$3"
RES="0p25"

# ------------------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------------------
[[ "$VAR" =~ ^(LAI|FPAR)$ ]] || {
  echo "ERROR: VAR must be LAI or FPAR" >&2
  exit 1
}
[[ "$MASKTAG" =~ ^(CCI|GLC)$ ]] || {
  echo "ERROR: MASKTAG must be CCI or GLC" >&2
  exit 1
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
IN_DIR="${ROOT}/output/${ALPHA}/masked_0p25/${VAR}/masked_${VAR}_${MASKTAG}"
OUT_EVAL="${ROOT}/output/${ALPHA}/eval/trend_${VAR}_${MASKTAG}"
mkdir -p "$OUT_EVAL"
LOG_FILE="${OUT_EVAL}/build_trends.log"
: > "$LOG_FILE"

# ------------------------------------------------------------------------------
# Thresholds
# ------------------------------------------------------------------------------
EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"
if [[ "$VAR" == "LAI" ]]; then
  EPS="$EPS_LAI"
else
  EPS="$EPS_FPAR"
fi

# ------------------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------------------
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}
log "========== START =========="
log "ALPHA=$ALPHA VAR=$VAR MASKTAG=$MASKTAG"

# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------
[[ -d "$IN_DIR" ]] || {
  log "ERROR: Missing input directory: $IN_DIR"
  exit 1
}
cd "$IN_DIR"
shopt -s nullglob
tifs=( "${VAR}_masked_"*.tif )
[[ ${#tifs[@]} -gt 0 ]] || {
  log "ERROR: No GeoTIFF files found"
  exit 1
}

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
MONTHLY="${VAR}_masked_monthly_${RES}.nc"
rm -f \
  "$MONTHLY" \
  "${VAR}_yearmean_${RES}.nc" \
  "${VAR}_yearmax_${RES}.nc" \
  "${VAR}_yearmin_${RES}.nc" \
  "${VAR}_yearamp_${RES}.nc" \
  "${VAR}"_*trend*_"${RES}.nc"
mkdir -p tmp_nc

# ------------------------------------------------------------------------------
# GeoTIFF → NetCDF
# ------------------------------------------------------------------------------
count=0
failed=0
for tif in "${tifs[@]}"; do
  count=$((count + 1))
  base=$(basename "$tif" ".tif")
  yyyymm=$(echo "$base" | grep -oE '[0-9]{6}' || true)
  if [[ ! "$yyyymm" =~ ^[0-9]{6}$ ]]; then
    log "WARNING: Invalid filename: $tif"
    failed=$((failed + 1))
    continue
  fi
  year="${yyyymm:0:4}"
  mon="${yyyymm:4:2}"
  raw_nc="tmp_nc/${VAR}_${yyyymm}.nc"
  dated_nc="tmp_nc/${VAR}_${yyyymm}_dated.nc"
  if ! gdal_translate -of NETCDF \
    "$tif" "$raw_nc" >>"$LOG_FILE" 2>&1; then
    log "WARNING: gdal_translate failed for $tif"
    failed=$((failed + 1))
    continue
  fi
  if ! cdo -O setdate,"${year}-${mon}-01" \
    "$raw_nc" "$dated_nc" >>"$LOG_FILE" 2>&1; then
    log "WARNING: cdo setdate failed for $tif"
    rm -f "$raw_nc"
    failed=$((failed + 1))
    continue
  fi
  rm -f "$raw_nc"
done
success=$((count - failed))
log "Converted $success/$count files"
[[ "$success" -gt 0 ]] || {
  log "ERROR: No files successfully converted"
  exit 1
}

# ------------------------------------------------------------------------------
# Merge monthly time series
# ------------------------------------------------------------------------------
mapfile -t dated_files < <(
  printf '%s\n' tmp_nc/"${VAR}"_*_dated.nc | sort
)
[[ ${#dated_files[@]} -gt 0 ]] || {
  log "ERROR: No dated NetCDF files found"
  exit 1
}
if ! cdo -O mergetime \
  "${dated_files[@]}" \
  "$MONTHLY" >>"$LOG_FILE" 2>&1; then
  log "ERROR: cdo mergetime failed"
  exit 1
fi

# ------------------------------------------------------------------------------
# Validate time axis
# ------------------------------------------------------------------------------
n_time=$(cdo -s ntime "$MONTHLY")
START_YEAR=1982
END_YEAR=2024
expected_months=$(( (END_YEAR - START_YEAR + 1) * 12 ))
if [[ "$n_time" -ne "$expected_months" ]]; then
  log "ERROR: Expected $expected_months timesteps, got $n_time"
  exit 1
fi
log "Validated monthly time axis ($n_time months)"

# ------------------------------------------------------------------------------
# Annual metrics
# ------------------------------------------------------------------------------
for met in yearmean yearmax yearmin; do
  cdo -O "$met" \
    "$MONTHLY" \
    "${VAR}_${met}_${RES}.nc" \
    >>"$LOG_FILE" 2>&1
done
cdo -O sub \
  "${VAR}_yearmax_${RES}.nc" \
  "${VAR}_yearmin_${RES}.nc" \
  "${VAR}_yearamp_${RES}.nc" \
  >>"$LOG_FILE" 2>&1

# ------------------------------------------------------------------------------
# OLS trends
# ------------------------------------------------------------------------------
for met in yearmean yearmax yearmin yearamp; do
  cdo -O trend \
    "${VAR}_${met}_${RES}.nc" \
    "${VAR}_${met}_trend_intercept_${RES}.nc" \
    "${VAR}_${met}_trend_slope_peryear_${RES}.nc" \
    >>"$LOG_FILE" 2>&1
done

# ------------------------------------------------------------------------------
# Relative trends (% yr^-1)
# ------------------------------------------------------------------------------
for met in yearmean yearmax yearmin yearamp; do
  ts="${VAR}_${met}_${RES}.nc"
  sl="${VAR}_${met}_trend_slope_peryear_${RES}.nc"
  mu="tmp_${VAR}_${met}_mean.nc"
  mask="tmp_${VAR}_${met}_mask.nc"
  out="${VAR}_${met}_trend_relative_percent_peryear_${RES}.nc"
  if [[ "$met" == "yearamp" ]]; then
    EPS_MET=0.01
  else
    EPS_MET="$EPS"
  fi
  cdo -O timmean "$ts" "$mu" >>"$LOG_FILE" 2>&1
  cdo -O gec,"${EPS_MET}" "$mu" "$mask" >>"$LOG_FILE" 2>&1
  cdo -O setname,trend_relative_percent \
    -ifthen "$mask" \
    -mulc,100 -div "$sl" "$mu" \
    "$out" >>"$LOG_FILE" 2>&1
  rm -f "$mu" "$mask"
done

# ------------------------------------------------------------------------------
# Move outputs
# ------------------------------------------------------------------------------
mv \
  "$MONTHLY" \
  "${VAR}_yearmean_${RES}.nc" \
  "${VAR}_yearmax_${RES}.nc" \
  "${VAR}_yearmin_${RES}.nc" \
  "${VAR}_yearamp_${RES}.nc" \
  "${VAR}_yearmean_trend_intercept_${RES}.nc" \
  "${VAR}_yearmean_trend_slope_peryear_${RES}.nc" \
  "${VAR}_yearmax_trend_intercept_${RES}.nc" \
  "${VAR}_yearmax_trend_slope_peryear_${RES}.nc" \
  "${VAR}_yearmin_trend_intercept_${RES}.nc" \
  "${VAR}_yearmin_trend_slope_peryear_${RES}.nc" \
  "${VAR}_yearamp_trend_intercept_${RES}.nc" \
  "${VAR}_yearamp_trend_slope_peryear_${RES}.nc" \
  "${VAR}_yearmean_trend_relative_percent_peryear_${RES}.nc" \
  "${VAR}_yearmax_trend_relative_percent_peryear_${RES}.nc" \
  "${VAR}_yearmin_trend_relative_percent_peryear_${RES}.nc" \
  "${VAR}_yearamp_trend_relative_percent_peryear_${RES}.nc" \
  "$OUT_EVAL/"

# ------------------------------------------------------------------------------
# Mann-Kendall
# ------------------------------------------------------------------------------
pids=()
for met in yearmean yearmax yearmin yearamp; do
  (
    cd "$ROOT"
    RUN_MODE=masked \
    RUN_TAG="$ALPHA" \
    MASK="$MASKTAG" \
    VAR="$VAR" \
    METRIC="$met" \
    Rscript "trends/compute_mk_pval.R"
  ) >>"$LOG_FILE" 2>&1 &
  pids+=($!)
done
mk_failed=0
for pid in "${pids[@]}"; do
  wait "$pid" || mk_failed=$((mk_failed + 1))
done
if [[ "$mk_failed" -gt 0 ]]; then
  log "WARNING: MK failed for $mk_failed metric(s)"
else
  log "Mann-Kendall complete"
fi
log "========== COMPLETE =========="
