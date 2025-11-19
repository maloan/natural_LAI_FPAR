#!/bin/bash
# =============================================================================
# rescale_fpar_all.sh — Fix double scaling of FPAR GeoTIFFs (×100)
# =============================================================================
set -euo pipefail

ROOT="$HOME/GitHub/natural_LAI_FPAR/output"
OUTSUFFIX="_rescaled"

# Find all FPAR GeoTIFFs in masked outputs (both 0p05 and 0p25)
find "$ROOT" -type f -path "*/masked_FPAR_*/*_0p*/FPAR_*_0p*.tif" | while read -r f; do
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

echo "Done. All rescaled GeoTIFFs written with suffix '${OUTSUFFIX}'."
