#!/bin/bash
## =============================================================================
# compute_gsmean_trends.sh — Growing-Season Mean & Trend from monthly 0.25° data
#
# Purpose:
#   Given a monthly 0.25° NetCDF file (e.g. LAI_georef_monthly_0p25.nc),
#   compute the growing-season mean (GSmean) for each year using a threshold,
#   then compute the linear trend per year.
#
# Usage:
#   ./compute_gsmean_trends.sh VAR THRESH
#   VAR    = LAI | FPAR
#   THRESH = growing season threshold (e.g. 0.5)
#
# Requirements:
#   cdo >= 1.9.8
#
## =============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Read arguments
# ------------------------------------------------------------------------------
VAR=${1:?Provide variable: LAI or FPAR}
THR=${2:?Provide threshold for GS detection (e.g. 0.5)}

echo "=========================================================="
echo "Growing Season Mean & Trend for: $VAR"
echo "Threshold for growing season: $THR"
echo "=========================================================="

# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
IN_DIR="$ROOT/analysis/unmasked/0p25"
OUT_DIR="$ROOT/analysis/unmasked/0p25"
mkdir -p "$OUT_DIR"

IN_MONTHLY="${IN_DIR}/${VAR}_georef_monthly_0p25.nc"

if [[ ! -f "$IN_MONTHLY" ]]; then
    echo "ERROR: Monthly file not found:"
    echo "  $IN_MONTHLY"
    exit 1
fi

echo "Input monthly file: $IN_MONTHLY"
echo "Output directory:   $OUT_DIR"
echo "----------------------------------------------------------"

# temporary directory
TMPDIR=$(mktemp -d)
echo "Using temp dir: $TMPDIR"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 1. Create binary mask of growing-season months (X > THRESH)
# ------------------------------------------------------------------------------
echo "Step 1: Mask months above threshold ($THR)"

cdo -O gtc,$THR "$IN_MONTHLY" "$TMPDIR/gs_mask.nc"
# result: 1 = GS month, 0 = non-GS month

# active values: x * mask
cdo -O mul "$IN_MONTHLY" "$TMPDIR/gs_mask.nc" "$TMPDIR/active.nc"

# ------------------------------------------------------------------------------
# 2. Yearly sum and yearly GS-month count
# ------------------------------------------------------------------------------
echo "Step 2: Yearly aggregation"

cdo -O yearsum "$TMPDIR/active.nc"      "$TMPDIR/gs_sum.nc"
cdo -O yearsum "$TMPDIR/gs_mask.nc"     "$TMPDIR/gs_nmonths.nc"

# ------------------------------------------------------------------------------
# 3. Growing-season mean = yearly_sum / n_months
# ------------------------------------------------------------------------------
echo "Step 3: Compute GS-mean (annual GSmean)"

GS_MEAN="${OUT_DIR}/${VAR}_georef_gsmean_0p25.nc"

cdo -O div "$TMPDIR/gs_sum.nc" "$TMPDIR/gs_nmonths.nc" "$GS_MEAN"

echo "GS-mean written to: $GS_MEAN"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# 4. Trend on GS-mean
# ------------------------------------------------------------------------------
echo "Step 4: Computing GS-mean trend"

TREND_INT="${OUT_DIR}/${VAR}_georef_gsmean_trend_intercept_0p25.nc"
TREND_SLP="${OUT_DIR}/${VAR}_georef_gsmean_trend_slope_peryear_0p25.nc"

cdo -O trend "$GS_MEAN" "$TREND_INT" "$TREND_SLP"

echo "GS-mean trend slope written to:   $TREND_SLP"
echo "GS-mean trend intercept written to: $TREND_INT"
echo "----------------------------------------------------------"

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
rm -rf "$TMPDIR"

echo "=========================================================="
echo "Growing-season processing completed for $VAR"
echo "Threshold used: $THR"
echo "GSmean file:    $GS_MEAN"
echo "Trend file:     $TREND_SLP"
echo "=========================================================="

