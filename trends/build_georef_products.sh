## =============================================================================
# build_georef_products.sh — Full LAI/FPAR georeferenced workflow (deduplicated)
#
# Usage
#   ./build_georef_products.sh LAI 0p05
#   ./build_georef_products.sh FPAR 0p05
## =============================================================================
set -euo pipefail
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------
VAR=${1:?Provide variable: LAI or FPAR}
RES005=${2:?Provide input resolution: e.g. 0p05}

case "$VAR" in
  LAI|FPAR) ;;
  *) echo "VAR must be LAI or FPAR, got: $VAR" >&2; exit 1 ;;
esac

echo "=========================================================="
echo "Building products for: $VAR at $RES005"
echo "=========================================================="

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
IN_DIR="$ROOT/output/nonvegetated_only_0p05/${VAR}"
OUT005="$ROOT/analysis/unmasked/${RES005}"
OUT025="$ROOT/analysis/unmasked/0p25"
REF025="$ROOT/src/ref_0p25.nc"
[[ -f "$REF025" ]] || { echo "Missing 0.25 reference grid: $REF025" >&2; exit 1; }

mkdir -p "$IN_DIR" "$OUT005" "$OUT025"

# Create CDO grid text files for remapping
# (CDO's remapcon doesn't read CRS metadata from NetCDF, so explicit grid files are needed)
GRID_005=$(mktemp)
GRID_025=$(mktemp)
cat > "$GRID_005" << 'GRIDEOF'
gridtype = lonlat
gridsize = 25920000
xsize    = 7200
ysize    = 3600
xname    = lon
xlongname = longitude
xunits   = degrees_east
yname    = lat
ylongname = latitude
yunits   = degrees_north
xfirst   = -179.975
xinc     = 0.05
yfirst   = -89.975
yinc     = 0.05
GRIDEOF

cat > "$GRID_025" << 'GRIDEOF'
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
GRIDEOF

echo "Input directory:       $IN_DIR"
echo "Output (0.05°):        $OUT005"
echo "Output (0.25°):        $OUT025"
echo "Reference 0.25° grid:  $REF025"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
need_cmd gdal_translate
need_cmd cdo

remap_to_025() {
  # remap_to_025 <in_nc> <out_nc>
  local in_nc="$1"
  local out_nc="$2"
  cdo -O remapcon,"$GRID_025" "$in_nc" "$out_nc"
}

set_taxis_if_ntime_matches() {
  # set_taxis_if_ntime_matches <file> <expected_ntime> <start_date> <step>
  # step examples: 1month | 1years
  local f="$1"
  local expected="$2"
  local start_date="$3"
  local step="$4"

  local ntime
  ntime=$(cdo ntime "$f" 2>/dev/null | awk '{print $NF}')
  if [[ "$ntime" -eq "$expected" ]]; then
    cdo -O settaxis,"$start_date",00:00:00,"$step" "$f" "${f%.nc}_time.nc"
    mv "${f%.nc}_time.nc" "$f"
  else
    echo "    Skipping time-axis reset for $(basename "$f"): expected $expected, got $ntime"
  fi
}

# ------------------------------------------------------------------------------
# 1) GeoTIFF to dated monthly NetCDF (tmp)
# ------------------------------------------------------------------------------
cd "$IN_DIR"

echo "Step 1: Converting GeoTIFFs → dated NetCDF (monthly)"

TMP_DIR="tmp_nc_${VAR}_${RES005}"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

shopt -s nullglob
tifs=( "${VAR}"_*_"${RES005}_masked_nonvegetated.tif" )
if [[ "${#tifs[@]}" -eq 0 ]]; then
  echo "No input GeoTIFFs found matching: ${VAR}_*_${RES005}_masked_nonvegetated.tif" >&2
  exit 1
fi

for tif in "${tifs[@]}"; do
  base=$(basename "$tif")
  # Expect: VAR_YYYYMM_RES005_masked_nonvegetated.tif
  yyyymm=${base#${VAR}_}
  yyyymm=${yyyymm%_${RES005}_masked_nonvegetated.tif}

  if [[ ! "$yyyymm" =~ ^[0-9]{6}$ ]]; then
    echo "Unexpected filename format: $base" >&2
    exit 1
  fi

  year=${yyyymm:0:4}
  mon=${yyyymm:4:2}

  raw_nc="${TMP_DIR}/${VAR}_${yyyymm}.nc"
  dated_nc="${TMP_DIR}/${VAR}_${yyyymm}_dated.nc"

  echo "  $base → $(basename "$dated_nc")"

  gdal_translate -of NETCDF "$tif" "$raw_nc" >/dev/null
  cdo -O setdate,"${year}-${mon}-01" "$raw_nc" "$dated_nc"
done

# ------------------------------------------------------------------------------
# 2) Merge to monthly time series
# ------------------------------------------------------------------------------
echo "Step 2: Building monthly stack"

MONTHLY_005="${VAR}_georef_monthly_${RES005}.nc"
cdo -O mergetime "${TMP_DIR}/${VAR}_"*_dated.nc "$MONTHLY_005"

echo "Monthly NetCDF created: $MONTHLY_005"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 3) Annual aggregates
# ------------------------------------------------------------------------------
echo "Step 3: Computing yearly aggregates"

YEARMEAN_005="${VAR}_georef_yearmean_${RES005}.nc"
YEARMAX_005="${VAR}_georef_yearmax_${RES005}.nc"
YEARMIN_005="${VAR}_georef_yearmin_${RES005}.nc"
YEARAMP_005="${VAR}_georef_yearamp_${RES005}.nc"

cdo -O yearmean "$MONTHLY_005" "$YEARMEAN_005"
cdo -O yearmax  "$MONTHLY_005" "$YEARMAX_005"
cdo -O yearmin  "$MONTHLY_005" "$YEARMIN_005"
cdo -O sub "$YEARMAX_005" "$YEARMIN_005" "$YEARAMP_005"

echo "Yearmean / yearmax / yearmin / yearamp created."
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 4) Trends
# ------------------------------------------------------------------------------
echo "Step 4: Computing trends (annual metrics)"

for metric_nc in "$YEARMEAN_005" "$YEARMAX_005" "$YEARMIN_005" "$YEARAMP_005"; do
  metric=$(basename "$metric_nc")
  metric=${metric#${VAR}_georef_}
  metric=${metric%_${RES005}.nc}

  echo "  Trend for: $metric"

  cdo -O trend "$metric_nc" \
    "${VAR}_georef_${metric}_trend_intercept_${RES005}.nc" \
    "${VAR}_georef_${metric}_trend_slope_peryear_${RES005}.nc"
done

echo "Trends computed."
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 5) Move 0.05° outputs to analysis folder
# ------------------------------------------------------------------------------
echo "Step 5: Moving 0.05° NetCDF outputs → $OUT005"

for f in \
  "$MONTHLY_005" \
  "$YEARMEAN_005" "$YEARMAX_005" "$YEARMIN_005" "$YEARAMP_005" \
  "${VAR}_georef_yearmean_trend_intercept_${RES005}.nc" \
  "${VAR}_georef_yearmean_trend_slope_peryear_${RES005}.nc" \
  "${VAR}_georef_yearmax_trend_intercept_${RES005}.nc" \
  "${VAR}_georef_yearmax_trend_slope_peryear_${RES005}.nc" \
  "${VAR}_georef_yearmin_trend_intercept_${RES005}.nc" \
  "${VAR}_georef_yearmin_trend_slope_peryear_${RES005}.nc" \
  "${VAR}_georef_yearamp_trend_intercept_${RES005}.nc" \
  "${VAR}_georef_yearamp_trend_slope_peryear_${RES005}.nc"; do
  [[ -f "$f" ]] && mv "$f" "$OUT005/"
done
rm -rf "$TMP_DIR"

# ------------------------------------------------------------------------------
# 6) Remap to 0.25° (annual metrics + slope fields)
# ------------------------------------------------------------------------------
echo "Step 6: Remapping to 0.25°"
cd "$OUT005"

START_YEAR=1982
END_YEAR=2024
N_YEARS=$((END_YEAR - START_YEAR + 1))

# (6b) Annual aggregates to 0.25° and enforce yearly axis
for metric in yearmean yearmax yearmin yearamp; do
  in_nc="${VAR}_georef_${metric}_${RES005}.nc"
  out_nc="${OUT025}/${VAR}_georef_${metric}_0p25.nc"
  remap_to_025 "$in_nc" "$out_nc"
  set_taxis_if_ntime_matches "$out_nc" "$N_YEARS" "${START_YEAR}-01-01" "1years"
done

# (6c) Trend slope to 0.25°
for metric in yearmean yearmax yearmin yearamp; do
  for kind in slope_peryear intercept; do
    in_nc="${VAR}_georef_${metric}_trend_${kind}_${RES005}.nc"
    [[ -f "$in_nc" ]] || continue
    out_nc="${OUT025}/${VAR}_georef_${metric}_trend_${kind}_0p25.nc"
    remap_to_025 "$in_nc" "$out_nc"
  done
done

# Clean up temporary grid files
rm -f "$GRID_005" "$GRID_025"
