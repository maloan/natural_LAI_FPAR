#!/bin/bash
# =============================================================================
# rescale_lai_all.sh — Fix double scaling of LAI GeoTIFFs (×100)
# =============================================================================
set -euo pipefail

ROOT="$HOME/GitHub/natural_LAI_FPAR/output"
OUTSUFFIX="_rescaled"

# Find all LAI GeoTIFFs (0.05° + 0.25°, all τ)
find "$ROOT" -type f -path "*/masked_LAI_*/*_0p*/LAI_*_0p*.tif" | while read -r f; do
  dir=$(dirname "$f")
  base=$(basename "$f")
  outdir="${dir}${OUTSUFFIX}"
  mkdir -p "$outdir"

  echo "→ Scaling: $base"
  gdal_calc.py -A "$f" \
    --outfile="$outdir/$base" \
    --calc="A*100" \
    --NoDataValue=-9999 \
    --type=Float32 \
    --creation-option="COMPRESS=DEFLATE" \
    --overwrite --quiet
done

echo "Done. All rescaled LAI GeoTIFFs written with suffix '${OUTSUFFIX}'."
