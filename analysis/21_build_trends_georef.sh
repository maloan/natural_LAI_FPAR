#!/bin/bash
## =============================================================================
# 21_build_trends_georef.sh — Trends for georeferenced LAI/FPAR (0.05° → 0.25°)
#
# Usage:
#   ./21_build_trends_georef.sh VAR
#   VAR ∈ {LAI, FPAR}
#
# Inputs:
#   data/georef/georef_<var>_0p05/<VAR>_YYYYMM_0p05.tif
#
# Outputs:
#   analysis/georef/<VAR>/0p05/*         (native)
#   analysis/georef/<VAR>/0p25/*         (regridded)
#
# Requires:
#   config/config.yml  (for grid_025 reference)
# =============================================================================

set -euo pipefail

VAR=${1:?Provide VAR=LAI or FPAR}
ROOT=${SNU_LAI_FPAR_ROOT:-"$HOME/GitHub/natural_LAI_FPAR"}
CFG="$ROOT/config/config.yml"

# Reference 0.25° grid (NetCDF, created in 00_setup.R)
REF025="$(yq -r '.grids.grid_025.ref_raster' "$CFG")"

INDIR="$ROOT/data/georef/georef_${VAR,,}_0p05"
OUT05="$ROOT/analysis/georef/${VAR}/0p05"
OUT25="$ROOT/analysis/georef/${VAR}/0p25"

mkdir -p "$OUT05" "$OUT25"

cd "$INDIR"

echo "=========================================================="
echo "Processing GEORF 0.05°  VAR=$VAR"
echo "Regridding to 0.25° afterwards"
echo "=========================================================="

# ------------------------------------------------------------
# 1. Convert all tif → netcdf
# ------------------------------------------------------------
echo "Converting GeoTIFF → NetCDF..."
for f in ${VAR}_??????_0p05.tif; do
    base="${f%.tif}"
    gdal_translate -q -of NETCDF "$f" "${base}.nc"
done

# ------------------------------------------------------------
# 2. Merge into monthly stack
# ------------------------------------------------------------
echo "Merging monthly files..."
cdo -O mergetime ${VAR}_??????_0p05.nc \
    "$OUT05/${VAR}_monthly_0p05.nc"

# ------------------------------------------------------------
# 3. Add monthly time axis
# ------------------------------------------------------------
echo "Adding monthly time axis..."
cdo settaxis,1982-01-01,00:00,1mon \
    "$OUT05/${VAR}_monthly_0p05.nc" \
    "$OUT05/${VAR}_monthly_0p05_time.nc"

# ------------------------------------------------------------
# 4. Yearly mean + yearly max
# ------------------------------------------------------------
echo "Computing yearly mean..."
cdo yearmean "$OUT05/${VAR}_monthly_0p05_time.nc" \
    "$OUT05/${VAR}_yearmean_0p05.nc"

echo "Computing yearly max..."
cdo yearmax "$OUT05/${VAR}_monthly_0p05_time.nc" \
    "$OUT05/${VAR}_yearmax_0p05.nc"

# ------------------------------------------------------------
# 5. Trends @ 0.05°
# ------------------------------------------------------------
echo "Trend (yearmean)..."
cdo trend "$OUT05/${VAR}_yearmean_0p05.nc" \
    "$OUT05/${VAR}_trend_intercept_yearmean_0p05.nc" \
    "$OUT05/${VAR}_trend_slope_yearmean_0p05.nc"

echo "Trend (yearmax)..."
cdo trend "$OUT05/${VAR}_yearmax_0p05.nc" \
    "$OUT05/${VAR}_trend_intercept_yearmax_0p05.nc" \
    "$OUT05/${VAR}_trend_slope_yearmax_0p05.nc"

# per decade
cdo mulc,10 "$OUT05/${VAR}_trend_slope_yearmean_0p05.nc" \
    "$OUT05/${VAR}_trend_slope_yearmean_decade_0p05.nc"
cdo mulc,10 "$OUT05/${VAR}_trend_slope_yearmax_0p05.nc" \
    "$OUT05/${VAR}_trend_slope_yearmax_decade_0p05.nc"

# ------------------------------------------------------------
# 6. Regrid 0.05° → 0.25° for all summary products
# ------------------------------------------------------------
echo "Regridding to 0.25° resolution..."

for src in \
    monthly_0p05 \
    monthly_0p05_time \
    yearmean_0p05 \
    yearmax_0p05 \
    trend_slope_yearmean_0p05 \
    trend_slope_yearmax_0p05 \
    trend_slope_yearmean_decade_0p05 \
    trend_slope_yearmax_decade_0p05 \
    trend_intercept_yearmean_0p05 \
    trend_intercept_yearmax_0p05; do

    INFILE="$OUT05/${VAR}_${src}.nc"
    OUTFILE="$OUT25/${VAR}_${src/0p05/0p25}.nc"

    echo " → remapbil: $(basename "$INFILE")"
    cdo remapbil,"$REF025" "$INFILE" "$OUTFILE"
done

# ------------------------------------------------------------
# 7. Export final 0.25° GeoTIFF for main decadal slope map
# ------------------------------------------------------------
echo "Exporting GeoTIFF (0.25° decadal trend, yearmean)..."

gdal_translate -of GTiff \
    NETCDF:"$OUT25/${VAR}_trend_slope_yearmean_decade_0p25.nc":Band1 \
    "$OUT25/${VAR}_trend_slope_yearmean_decade_0p25.tif"

echo "DONE for VAR=$VAR (georef, 0.05° → 0.25°)"
