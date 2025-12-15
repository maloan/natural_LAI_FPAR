#!/bin/bash
## =============================================================================
# remap_monthly_to_0p25.sh — Remap monthly 0.05° NetCDF to 0.25° and restore time
#
# Usage:
#   ./remap_monthly_to_0p25.sh VAR RES005
#
# Arguments:
#   VAR     = LAI | FPAR
#   RES005  = 0p05   (resolution suffix of input monthly file)
#
# Input expected at:
#   $ROOT/analysis/unmasked/$RES005/${VAR}_georef_monthly_${RES005}.nc
#
# Output:
#   $ROOT/analysis/unmasked/0p25/${VAR}_georef_monthly_0p25.nc
#
# Requirements:
#   cdo, bash, NetCDF tools
#
## =============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------
VAR=${1:?Provide variable: LAI or FPAR}
RES005=${2:?Provide input resolution, e.g. 0p05}

echo "=============================================================="
echo "Remapping monthly product for: $VAR at input resolution $RES005"
echo "=============================================================="

# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
IN_DIR="$ROOT/analysis/unmasked/${RES005}"
OUT025="$ROOT/analysis/unmasked/0p25"
REF025="$ROOT/src/ref_0p25.nc"

mkdir -p "$OUT025"

echo "Input directory:  $IN_DIR"
echo "Output directory: $OUT025"
echo "Reference grid:   $REF025"
echo "--------------------------------------------------------------"

# ------------------------------------------------------------------------------
# Input/output filenames
# ------------------------------------------------------------------------------
monthly_in="${IN_DIR}/${VAR}_georef_monthly_${RES005}.nc"
monthly_out="${OUT025}/${VAR}_georef_monthly_0p25.nc"

if [[ ! -f "$monthly_in" ]]; then
    echo "ERROR: monthly input file not found:"
    echo "  $monthly_in"
    exit 1
fi

echo "Found monthly input: $monthly_in"
echo "--------------------------------------------------------------"

# ------------------------------------------------------------------------------
# Step 1: Remap monthly to 0.25°
# ------------------------------------------------------------------------------
echo "Step 1: Remapping monthly file to 0.25°"

cdo -O remapcon,$REF025 "$monthly_in" "$monthly_out"

echo "Remapped file created:"
echo "  $monthly_out"
echo "--------------------------------------------------------------"

# ------------------------------------------------------------------------------
# Step 2: Restore monthly time axis (fixed monthly sequence)
# ------------------------------------------------------------------------------
echo "Step 2: Restoring monthly time axis"

tmp="${monthly_out%.nc}_time.nc"

cdo -settaxis,1982-01-01,00:00:00,1month \
    "$monthly_out" \
    "$tmp"

mv "$tmp" "$monthly_out"

echo "Time axis restored in: $monthly_out"
echo "--------------------------------------------------------------"

echo "DONE."
echo "=============================================================="

