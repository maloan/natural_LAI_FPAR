#!/bin/bash
## =============================================================================
# build_georef_products.sh — Full LAI/FPAR georeferenced workflow
#
# Purpose
#   Convert monthly GeoTIFFs (0.05°) into dated NetCDF, merge into a proper
#   time series, compute yearly aggregates, compute trends, and remap to 0.25°.
#
# Usage
#   ./build_georef_products.sh LAI 0p05
#   ./build_georef_products.sh FPAR 0p05
#
# Arguments
#   VAR     = LAI | FPAR
#   RES005  = 0p05  (input resolution)
#
# Requirements
#   gdal_translate, cdo
#
## =============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------
VAR=${1:?Provide variable: LAI or FPAR}
RES005=${2:?Provide input resolution: e.g. 0p05}

echo "=========================================================="
echo "Building products for: $VAR at $RES005"
echo "=========================================================="

# ------------------------------------------------------------------------------
# Define directories
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
IN_DIR="$ROOT/output/georef_biotic/georef_biotic_${VAR}_${RES005}"
OUT005="$ROOT/analysis/unmasked/${RES005}"
OUT025="$ROOT/analysis/unmasked/0p25"
REF025="$ROOT/src/ref_0p25.nc"

mkdir -p "$OUT005" "$OUT025"
mkdir -p "$IN_DIR"
echo "Input directory:       $IN_DIR"
echo "Output (0.05°):        $OUT005"
echo "Output (0.25°):        $OUT025"
echo "Reference 0.25° grid:  $REF025"
echo "----------------------------------------------------------"

cd "$IN_DIR"

# ------------------------------------------------------------------------------
# 1. Convert GeoTIFF → NetCDF with timestamps
# ------------------------------------------------------------------------------
echo "Step 1: Converting GeoTIFFs → NetCDF with time axis"

rm -rf tmp_nc
mkdir -p tmp_nc

for tif in ${VAR}_*.tif; do
    base=$(basename "$tif")
    yyyymm=${base#${VAR}_}
    yyyymm=${yyyymm%_${RES005}.tif}

    year=${yyyymm:0:4}
    mon=${yyyymm:4:2}

    raw_nc="tmp_nc/${VAR}_${yyyymm}.nc"
    dated_nc="tmp_nc/${VAR}_${yyyymm}_dated.nc"

    echo "  $tif → $dated_nc"

    # GeoTIFF → raw NetCDF
    gdal_translate -of NETCDF "$tif" "$raw_nc"

    # Add timestamp
    cdo -O setdate,${year}-${mon}-01 "$raw_nc" "$dated_nc"
done

echo "Done converting GeoTIFFs."
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 2. Merge timestamps into a monthly time series stack
# ------------------------------------------------------------------------------
echo "Step 2: Building monthly stack"
cdo -O mergetime tmp_nc/${VAR}_*_dated.nc ${VAR}_georef_monthly_${RES005}.nc
echo "Monthly NetCDF created: ${VAR}_georef_monthly_${RES005}.nc"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 3. Yearly aggregates: mean, max, min, amplitude
# ------------------------------------------------------------------------------
echo "Step 3: Computing yearly aggregates"

cdo -O yearmean ${VAR}_georef_monthly_${RES005}.nc ${VAR}_georef_yearmean_${RES005}.nc
cdo -O yearmax  ${VAR}_georef_monthly_${RES005}.nc ${VAR}_georef_yearmax_${RES005}.nc
cdo -O yearmin  ${VAR}_georef_monthly_${RES005}.nc ${VAR}_georef_yearmin_${RES005}.nc

cdo -O sub ${VAR}_georef_yearmax_${RES005}.nc \
           ${VAR}_georef_yearmin_${RES005}.nc \
           ${VAR}_georef_yearamp_${RES005}.nc

echo "Yearmean / yearmax / yearmin / yearamp created."
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 4. Compute trends (per year)
# ------------------------------------------------------------------------------
echo "Step 4: Computing trends"

for metric in yearmean yearmax yearmin yearamp monthly; do

    echo "  Trend for: $metric"

    cdo -O trend ${VAR}_georef_${metric}_${RES005}.nc \
        ${VAR}_georef_${metric}_trend_intercept_${RES005}.nc \
        ${VAR}_georef_${metric}_trend_slope_peryear_${RES005}.nc

done

echo "Trends computed."
echo "----------------------------------------------------------"


# ------------------------------------------------------------------------------
# 5. Move all 0p05 files to analysis folder
# ------------------------------------------------------------------------------
echo "Step 5: Moving 0.05° products → $OUT005"

mv *.nc "$OUT005"/
rm -rf tmp_nc

echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 5b. Remap monthly dataset to 0.25° and restore monthly time axis
# ------------------------------------------------------------------------------

cd "$OUT005"

echo "Step 5b: Remapping MONTHLY product to 0.25°"

monthly_in="${VAR}_georef_monthly_${RES005}.nc"
monthly_out="${VAR}_georef_monthly_0p25.nc"

if [[ -f "$monthly_in" ]]; then

    echo "  Remapping ${monthly_in} → ${OUT025}/${monthly_out}"

    # Conservative remap to 0.25°
    cdo -O remapcon,$REF025 "$monthly_in" "${OUT025}/${monthly_out}"

    # Restore monthly time axis
    echo "  Restoring monthly time axis"

    cdo -settaxis,1982-01-01,00:00:00,1month \
        "${OUT025}/${monthly_out}" \
        "${OUT025}/${monthly_out%.nc}_time.nc"

    mv "${OUT025}/${monthly_out%.nc}_time.nc" "${OUT025}/${monthly_out}"

else
    echo "  No monthly file found: $monthly_in"
fi

echo "Monthly 0.25° file created: ${OUT025}/${monthly_out}"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 6. Remap selected fields to 0.25° and restore time axis
# ------------------------------------------------------------------------------
cd "$OUT005"

echo "Step 6: Remapping to 0.25° resolution and restoring time metadata"

declare -a FILES_TO_REMAP=(
    "${VAR}_georef_yearmean_${RES005}.nc"
    "${VAR}_georef_yearmax_${RES005}.nc"
    "${VAR}_georef_yearmin_${RES005}.nc"
    "${VAR}_georef_yearamp_${RES005}.nc"
)

START_YEAR=1982
END_YEAR=2024
EXPECTED_LAYERS=$((END_YEAR - START_YEAR + 1))

for f in "${FILES_TO_REMAP[@]}"; do
    out="${f/${RES005}/0p25}"
    echo "  Remapping $f → $OUT025/$out"
    cdo -O remapcon,$REF025 "$f" "$OUT025/$out"

    # Check number of layers
    n_layers=$(cdo ntime "$f" 2>/dev/null | awk '{print $NF}')

    if [ "$n_layers" -eq "$EXPECTED_LAYERS" ]; then
        echo "    Restoring time axis on $out"
        cdo settaxis,${START_YEAR}-01-01,00:00:00,1years \
            "$OUT025/$out" \
            "$OUT025/${out%.nc}_time.nc"
        mv "$OUT025/${out%.nc}_time.nc" "$OUT025/$out"
    else
        echo "    Skipping time restoration: expected $EXPECTED_LAYERS layers, got $n_layers"
    fi
done

echo "Step 6b: Remapping slope (trend) rasters to 0.25°"

SLOPE_FILES=(
    "${VAR}_georef_yearmean_trend_slope_peryear_${RES005}.nc"
    "${VAR}_georef_yearmax_trend_slope_peryear_${RES005}.nc"
    "${VAR}_georef_yearamp_trend_slope_peryear_${RES005}.nc"
    "${VAR}_georef_yearmin_trend_slope_peryear_${RES005}.nc"
)

for f in "${SLOPE_FILES[@]}"; do
    if [[ -f "$f" ]]; then
        out="${f/${RES005}/0p25}"
        echo "  Remapping slope: $f → $OUT025/$out"
        cdo -O remapcon,$REF025 "$f" "$OUT025/$out"
    fi
done


echo "----------------------------------------------------------"
echo "All tasks completed."
echo "Products:"
echo "  - $OUT005  (0.05° full outputs)"
echo "  - $OUT025  (0.25° remapped analysis datasets)"
echo "=========================================================="
