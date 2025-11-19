#!/bin/bash
## =============================================================================
# 20_build_trends.sh — Build monthly → yearly → trend products for all outputs
#
# Usage
#   ./20_build_trends.sh tau_0.1
#
# Inputs
#   output/<τ>/masked_0p25/<VAR>/masked_<VAR>_{CCI|GLC}/masked_<var>_0p25/*.tif
#
# Outputs (per mask + variable)
#   - monthly stack:    <var>_masked_monthly_0p25.nc
#   - yearly products:  <var>_yearmean_0p25.nc, <var>_yearmax_0p25.nc
#   - trends:           <var>_trend_slope_*.nc, <var>_trend_intercept_*.nc
#   - GeoTIFF export:   slope_decade.tif
#
# Notes
#   - Reads the reference 0.25° NC from config.yml (already created in 00_setup.R)
#   - Assumes filenames follow: VAR_masked_YYYYMM_0p25.tif
## =============================================================================

set -euo pipefail

TAU=${1:?Provide tau, example: tau_0.1}
ROOT=${SNU_LAI_FPAR_ROOT:-"$HOME/GitHub/natural_LAI_FPAR"}

CFG="$ROOT/config/config.yml"
REF025="$(yq -r '.grids.grid_025.ref_raster' "$CFG")"

# ----- variables & mask types --------------------------------------------------
VARS=("FPAR" "LAI")
MASKS=("CCI" "GLC")

# ----- loop --------------------------------------------------------------------
for VAR in "${VARS[@]}"; do
  for MASK in "${MASKS[@]}"; do

    echo "=============================================================="
    echo "Processing VAR=$VAR MASK=$MASK τ=$TAU"
    echo "=============================================================="

    INDIR="$ROOT/output/$TAU/masked_0p25/$VAR/masked_${VAR}_${MASK}/masked_${var,,}_0p25"
    OUTDIR="$ROOT/analysis/${TAU}_0p25/${VAR}_${MASK}"
    mkdir -p "$OUTDIR"

    cd "$INDIR"

    # --- Step 1: convert all tif → nc ----------------------------------------
    echo "Converting GeoTIFF → NetCDF..."
    for f in ${VAR}_masked_??????_0p25.tif; do
        base="${f%.tif}"
        gdal_translate -q -of NETCDF "$f" "${base}.nc"
    done

    # --- Step 2: merge time ---------------------------------------------------
    echo "Merging to monthly time-series..."
    cdo -O mergetime ${VAR}_masked_??????_0p25.nc \
        "$OUTDIR/${VAR}_masked_monthly_0p25.nc"

    # --- Step 3: add time axis ------------------------------------------------
    echo "Adding monthly time axis (starting 1982-01-01)..."
    cdo settaxis,1982-01-01,00:00,1mon \
        "$OUTDIR/${VAR}_masked_monthly_0p25.nc" \
        "$OUTDIR/${VAR}_masked_monthly_0p25_time.nc"

    # --- Step 4: yearly products ---------------------------------------------
    echo "Computing yearly mean..."
    cdo yearmean "$OUTDIR/${VAR}_masked_monthly_0p25_time.nc" \
        "$OUTDIR/${VAR}_yearmean_0p25.nc"

    echo "Computing yearly max..."
    cdo yearmax "$OUTDIR/${VAR}_masked_monthly_0p25_time.nc" \
        "$OUTDIR/${VAR}_yearmax_0p25.nc"

    # --- Step 5: trends -------------------------------------------------------
    echo "Computing trends for yearly mean..."
    cdo trend "$OUTDIR/${VAR}_yearmean_0p25.nc" \
        "$OUTDIR/${VAR}_trend_intercept_yearmean.nc" \
        "$OUTDIR/${VAR}_trend_slope_yearmean.nc"

    echo "Computing trends for yearly max..."
    cdo trend "$OUTDIR/${VAR}_yearmax_0p25.nc" \
        "$OUTDIR/${VAR}_trend_intercept_yearmax.nc" \
        "$OUTDIR/${VAR}_trend_slope_yearmax.nc"

    # per decade
    cdo mulc,10 "$OUTDIR/${VAR}_trend_slope_yearmean.nc" \
        "$OUTDIR/${VAR}_trend_slope_yearmean_decade.nc"
    cdo mulc,10 "$OUTDIR/${VAR}_trend_slope_yearmax.nc" \
        "$OUTDIR/${VAR}_trend_slope_yearmax_decade.nc"

    # --- Step 6: export GeoTIFF ----------------------------------------------
    echo "Exporting decadal mean-trend GeoTIFF..."
    gdal_translate -of GTiff \
        NETCDF:"$OUTDIR/${VAR}_trend_slope_yearmean_decade.nc":${VAR} \
        "$OUTDIR/${VAR}_trend_slope_yearmean_decade.tif"

    echo "Finished VAR=$VAR MASK=$MASK"

  done
done

echo "All trend products completed."
