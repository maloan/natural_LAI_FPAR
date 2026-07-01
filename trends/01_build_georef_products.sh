#!/usr/bin/env bash
# ==============================================================================
# 01_build_georef_products.sh - Full LAI/FPAR georeferenced workflow
# (GeoTIFF → trends → MK p-values → relative trends)
# Usage:
# ./01_build_georef_products.sh LAI 0p05
# ./01_build_georef_products.sh FPAR 0p05
# Or run both with VARS="LAI FPAR" bash ./01_build_georef_products.sh 0p05
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
cleanup() {
  rm -rf tmp_nc_*
  rm -f "${GRID_025:-}"
}

trap cleanup EXIT
trap 'echo "ERROR at line $LINENO" >&2; exit 1' ERR

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing dependency: $1" >&2
    exit 1
  }
}

need_cmd gdal_translate
need_cmd cdo
need_cmd Rscript

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------
if [[ $# -eq 2 ]]; then
  VAR="$1"
  RES005="$2"
  VARS_ARR=("$VAR")
elif [[ $# -eq 1 ]]; then
  RES005="$1"
  VARS_ARR=(${VARS:-"LAI FPAR"})
else
  cat >&2 << EOF
Usage (single variable):
  ./01_build_georef_products.sh LAI 0p05
Usage (both variables):
  VARS="LAI FPAR" bash ./01_build_georef_products.sh 0p05
EOF
  exit 1
fi

for VAR in "${VARS_ARR[@]}"; do
  case "$VAR" in
    LAI|FPAR) ;;
    *)
      echo "VAR must be LAI or FPAR, got: $VAR" >&2
      exit 1
      ;;
  esac
done

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"
START_YEAR=1982
END_YEAR=2024
N_YEARS=$((END_YEAR - START_YEAR + 1))

# ------------------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------------------
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}
log_error() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

# ------------------------------------------------------------------------------
# Grid definition
# ------------------------------------------------------------------------------
GRID_025=$(mktemp)
cat > "$GRID_025" << 'EOF'
gridtype = lonlat
gridsize = 1036800
xsize    = 1440
ysize    = 720
xname    = longitude
xlongname = longitude
xunits   = degrees_east
yname    = latitude
ylongname = latitude
yunits   = degrees_north
xfirst   = -179.875
xinc     = 0.25
yfirst   = -89.875
yinc     = 0.25
EOF

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
remap_to_025() {
  local in_nc="$1"
  local out_nc="$2"
  cdo -O remapcon,"$GRID_025" "$in_nc" "$out_nc"
}

set_taxis_if_ntime_matches() {
  local f="$1"
  local expected="$2"
  local start_date="$3"
  local step="$4"
  local ntime
  ntime=$(cdo -s ntime "$f")
  if [[ "$ntime" -eq "$expected" ]]; then
    cdo -O settaxis,"$start_date",00:00:00,"$step" \
      "$f" "${f%.nc}_time.nc"
    mv "${f%.nc}_time.nc" "$f"
  else
    echo "Skipping time-axis reset for $(basename "$f"): expected $expected, got $ntime"
  fi
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------
for VAR in "${VARS_ARR[@]}"; do
  echo ""
  echo "=========================================================="
  echo "Processing: $VAR at $RES005"
  echo "=========================================================="
  case "$VAR" in
    LAI)  EPS="$EPS_LAI" ;;
    FPAR) EPS="$EPS_FPAR" ;;
  esac

  # ---------------------------------------------------------------------------
  # Paths
  # ---------------------------------------------------------------------------
  IN_DIR="$ROOT/output/nonvegetated_only_0p05/${VAR}"
  OUT005="$ROOT/analysis/unmasked/${RES005}"
  OUT025="$ROOT/analysis/unmasked/0p25"

  mkdir -p "$OUT005" "$OUT025"
  cd "$IN_DIR" || {
    log_error "Cannot change to input directory: $IN_DIR"
    exit 1
  }

  # ------------------------------------------------------------------------------
  # GeoTIFF → CF-compliant NetCDF (robust CDO pipeline)
  # ------------------------------------------------------------------------------

  TMP_DIR="tmp_nc_${VAR}_${RES005}"
  rm -rf "$TMP_DIR"
  mkdir -p "$TMP_DIR"

  shopt -s nullglob
  tifs=( "${VAR}"_*_"${RES005}_masked_nonvegetated.tif" )

  [[ ${#tifs[@]} -gt 0 ]] || {
    echo "No input GeoTIFFs found" >&2
    exit 1
  }

  for tif in "${tifs[@]}"; do
    base=$(basename "$tif")
    conversion_count=$((conversion_count + 1))

    yyyymm=$(echo "$base" | grep -oE '[0-9]{6}' | head -1 || true)
    if [[ ! "$yyyymm" =~ ^[0-9]{6}$ ]]; then
      log_error "Bad filename format: $base (expected YYYYMM in name)"
      conversion_failed=$((conversion_failed + 1))
      continue
    fi

    year=${yyyymm:0:4}
    mon=${yyyymm:4:2}
    if ! [[ "$year" =~ ^[0-9]{4}$ ]] || ! [[ "$mon" =~ ^[0-9]{2}$ ]] || (( mon < 1 || mon > 12 )); then
      log_error "Invalid year/month from $base: $year-$mon"
      conversion_failed=$((conversion_failed + 1))
      continue
    fi

    raw_nc="${TMP_DIR}/${VAR}_${yyyymm}_raw.nc"
    out_nc="${TMP_DIR}/${VAR}_${yyyymm}_dated.nc"

    # ---------------------------------------------------------------------------
    # Step 1: GeoTIFF → NetCDF via GDAL
    # ---------------------------------------------------------------------------
    if ! gdal_translate -of NetCDF "$tif" "$raw_nc" >/dev/null 2>&1; then
      log_error "gdal_translate failed for $tif"
      conversion_failed=$((conversion_failed + 1))
      continue
    fi

    # ---------------------------------------------------------------------------
    # Step 2: assign CF-compliant time axis
    # ---------------------------------------------------------------------------
    if ! cdo -O settaxis,${year}-${mon}-01,00:00:00,1month \
      "$raw_nc" "$out_nc" >/dev/null 2>&1; then
      log_error "cdo settaxis failed for $tif"
      conversion_failed=$((conversion_failed + 1))
      rm -f "$raw_nc"
      continue
    fi

    rm -f "$raw_nc"
  done
  
  conversion_count=${#tifs[@]}
  conversion_success=$((conversion_count - conversion_failed))
  log "GeoTIFF conversion: $conversion_success/$conversion_count files succeeded"
  if [[ "$conversion_success" -eq 0 ]]; then
    log_error "No files successfully converted"
    exit 1
  fi
  # ---------------------------------------------------------------------------
  # Merge monthly time series (CF-safe)
  # ---------------------------------------------------------------------------

  MONTHLY_005="${VAR}_georef_monthly_${RES005}.nc"

  mapfile -t dated_files < <(ls "${TMP_DIR}"/*_dated.nc | sort)

  [[ ${#dated_files[@]} -gt 0 ]] || {
    echo "No dated NetCDF files found" >&2
    exit 1
  }

  TMP_MERGE="${TMP_DIR}/tmp_merge.nc"

  if ! cdo -O mergetime "${dated_files[@]}" "$TMP_MERGE" >/dev/null 2>&1; then
    log_error "cdo mergetime failed"
    exit 1
  fi

  if ! cdo -O settaxis,${START_YEAR}-01-01,00:00:00,1month \
    "$TMP_MERGE" "$MONTHLY_005" >/dev/null 2>&1; then
    log_error "cdo settaxis on merged file failed"
    exit 1
  fi

  rm -f "$TMP_MERGE"
  
  n_time=$(cdo -s ntime "$MONTHLY_005" 2>/dev/null || echo 0)
  expected_months=$(( (END_YEAR - START_YEAR + 1) * 12 ))
  if [[ "$n_time" -ne "$expected_months" ]]; then
    log_error "Expected $expected_months monthly timesteps, got $n_time"
    exit 1
  fi
  log "Validated monthly time series: $n_time timesteps"


  # ---------------------------------------------------------------------------
  # Annual metrics
  # ---------------------------------------------------------------------------
  YEARMEAN_005="${VAR}_georef_yearmean_${RES005}.nc"
  YEARMAX_005="${VAR}_georef_yearmax_${RES005}.nc"
  YEARMIN_005="${VAR}_georef_yearmin_${RES005}.nc"
  YEARAMP_005="${VAR}_georef_yearamp_${RES005}.nc"
  
  if ! cdo -O yearmean "$MONTHLY_005" "$YEARMEAN_005" >/dev/null 2>&1; then
    log_error "cdo yearmean failed"
    exit 1
  fi
  if ! cdo -O yearmax "$MONTHLY_005" "$YEARMAX_005" >/dev/null 2>&1; then
    log_error "cdo yearmax failed"
    exit 1
  fi
  if ! cdo -O yearmin "$MONTHLY_005" "$YEARMIN_005" >/dev/null 2>&1; then
    log_error "cdo yearmin failed"
    exit 1
  fi
  if ! cdo -O sub "$YEARMAX_005" "$YEARMIN_005" "$YEARAMP_005" >/dev/null 2>&1; then
    log_error "cdo sub (yearamp) failed"
    exit 1
  fi
  
  log "Annual metrics computed: yearmean, yearmax, yearmin, yearamp"

  # ---------------------------------------------------------------------------
  # OLS trends
  # ---------------------------------------------------------------------------
  for metric_nc in \
    "$YEARMEAN_005" \
    "$YEARMAX_005" \
    "$YEARMIN_005" \
    "$YEARAMP_005"; do
    metric=$(basename "$metric_nc")
    metric=${metric#${VAR}_georef_}
    metric=${metric%_${RES005}.nc}
    intercept_file="${VAR}_georef_${metric}_trend_intercept_${RES005}.nc"
    slope_file="${VAR}_georef_${metric}_trend_slope_peryear_${RES005}.nc"
    
    if ! cdo -O trend "$metric_nc" "$intercept_file" "$slope_file" >/dev/null 2>&1; then
      log_error "cdo trend failed for metric $metric"
      exit 1
    fi
  done
  log "OLS trends computed for all metrics"

  # ---------------------------------------------------------------------------
  # Move 0.05° outputs
  # ---------------------------------------------------------------------------
  mv \
    "$MONTHLY_005" \
    "$YEARMEAN_005" \
    "$YEARMAX_005" \
    "$YEARMIN_005" \
    "$YEARAMP_005" \
    "${VAR}_georef_yearmean_trend_intercept_${RES005}.nc" \
    "${VAR}_georef_yearmean_trend_slope_peryear_${RES005}.nc" \
    "${VAR}_georef_yearmax_trend_intercept_${RES005}.nc" \
    "${VAR}_georef_yearmax_trend_slope_peryear_${RES005}.nc" \
    "${VAR}_georef_yearmin_trend_intercept_${RES005}.nc" \
    "${VAR}_georef_yearmin_trend_slope_peryear_${RES005}.nc" \
    "${VAR}_georef_yearamp_trend_intercept_${RES005}.nc" \
    "${VAR}_georef_yearamp_trend_slope_peryear_${RES005}.nc" \
    "$OUT005/"
  rm -rf "$TMP_DIR"

  # ---------------------------------------------------------------------------
  # Remap to 0.25°
  # ---------------------------------------------------------------------------
  cd "$OUT005"
  for metric in yearmean yearmax yearmin yearamp; do
    in_nc="${VAR}_georef_${metric}_${RES005}.nc"
    out_nc="${OUT025}/${VAR}_georef_${metric}_0p25.nc"
    remap_to_025 "$in_nc" "$out_nc"
    set_taxis_if_ntime_matches \
      "$out_nc" \
      "$N_YEARS" \
      "${START_YEAR}-01-01" \
      "1years"
  done

  for metric in yearmean yearmax yearmin yearamp; do
    for kind in slope_peryear intercept; do
      in_nc="${VAR}_georef_${metric}_trend_${kind}_${RES005}.nc"
      [[ -f "$in_nc" ]] || continue
      out_nc="${OUT025}/${VAR}_georef_${metric}_trend_${kind}_0p25.nc"
      remap_to_025 "$in_nc" "$out_nc"
    done
  done

  # ---------------------------------------------------------------------------
  # Mann-Kendall p-values (parallel)
  # ---------------------------------------------------------------------------
  mk_pids=()
  for metric in yearmean yearmax yearmin yearamp; do
    (
      cd "$ROOT"
      RUN_MODE=unmasked \
      VAR="$VAR" \
      METRIC="$metric" \
      Rscript "trends/compute_mk_pval.R"
    ) &
    mk_pids+=("$!")
  done
  
  mk_failed=0
  for pid in "${mk_pids[@]}"; do
    if ! wait "$pid"; then
      mk_failed=$((mk_failed + 1))
    fi
  done
  
  if [[ "$mk_failed" -gt 0 ]]; then
    log_error "Mann-Kendall test failed for $mk_failed metric(s)"
    exit 1
  fi
  log "Mann-Kendall p-values computed for all metrics"

  # ---------------------------------------------------------------------------
  # Relative trends (% yr^-1)
  # ---------------------------------------------------------------------------
  cd "$OUT025"
  for MET in yearmean yearmax yearmin yearamp; do
    ts="${VAR}_georef_${MET}_0p25.nc"
    sl="${VAR}_georef_${MET}_trend_slope_peryear_0p25.nc"
    mu="tmp_${VAR}_${MET}_mean.nc"
    mask="tmp_${VAR}_${MET}_mask.nc"
    out="${VAR}_georef_${MET}_trend_relative_percent_peryear_0p25.nc"
    [[ -f "$ts" && -f "$sl" ]] || continue
    if [[ "$MET" == "yearamp" ]]; then
      EPS_MET=0.01
    else
      EPS_MET="$EPS"
    fi
    cdo -O timmean "$ts" "$mu"
    cdo -O gec,"${EPS_MET}" "$mu" "$mask"
    cdo -O setname,trend_relative_percent \
      -ifthen "$mask" \
      -mulc,100 -div "$sl" "$mu" \
      "$out"
    rm -f "$mu" "$mask"
  done
  echo "Completed: $VAR"
done
echo ""
echo "All variables processed successfully"
