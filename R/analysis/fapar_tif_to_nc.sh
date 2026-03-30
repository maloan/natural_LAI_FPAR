#!/usr/bin/env bash
# ==============================================================================
# fapar_tif_to_nc.sh
# Convert monthly masked fPAR GeoTIFFs (0.5°) -> single monthly NetCDF (1982–2024)
# ==============================================================================

set -euo pipefail

OUT="FPAR_masked_0p5_monthly_1982-2024.nc"
REF_DATE="1982-01-01"
TMPDIR="${TMPDIR:-./_tmp_fpar_nc}"

mkdir -p "$TMPDIR" "$TMPDIR/nc_monthly" "$TMPDIR/nc_monthly_notime"
rm -f "$OUT"

# --- optional safety copy
if [[ -f "$OUT" ]]; then
  cp -p "$OUT" "${OUT%.nc}_OLD.nc"
fi

# ------------------------------------------------------------------------------
# 1) GeoTIFF -> per-month NetCDF (one timestep each)
# ------------------------------------------------------------------------------
rm -f "$TMPDIR/nc_monthly"/FPAR_masked_*_0p5.nc
ls -1 FPAR_masked_*_0p5.tif >/dev/null 2>&1 || {
  echo "ERROR: No inputs matching FPAR_masked_*_0p5.tif in $(pwd)" >&2
  exit 1
}

while IFS= read -r f; do
  gdal_translate -q -of NetCDF -co FORMAT=NC4 -co COMPRESS=DEFLATE -co ZLEVEL=4 \
    "$f" "$TMPDIR/nc_monthly/${f%.tif}.nc"
done < <(ls -1 FPAR_masked_*_0p5.tif | sort)

# ------------------------------------------------------------------------------
# 2) Remove GDAL's bad coordinate variable 'time' (keep the time dimension)
# ------------------------------------------------------------------------------
rm -f "$TMPDIR/nc_monthly_notime"/FPAR_masked_*_0p5.nc
while IFS= read -r f; do
  ncks -O -C -x -v time \
    "$f" "$TMPDIR/nc_monthly_notime/$(basename "$f")"
done < <(ls -1 "$TMPDIR/nc_monthly"/FPAR_masked_*_0p5.nc | sort)

# ------------------------------------------------------------------------------
# 3) Concatenate and create a fresh monthly time axis
#    (CDO often stores time internally as hours; we'll convert to days)
# ------------------------------------------------------------------------------
cdo -O \
  -setcalendar,standard \
  -settaxis,"${REF_DATE}",00:00:00,1mon \
  -cat "$TMPDIR"/nc_monthly_notime/FPAR_masked_*_0p5.nc \
  "$TMPDIR/concat_hours.nc"

# ------------------------------------------------------------------------------
# 4) Convert time values: hours -> days, enforce CF attributes
# ------------------------------------------------------------------------------
ncap2 -O -s 'time=time/24.0' \
  "$TMPDIR/concat_hours.nc" "$TMPDIR/axis_fixed.nc"
rm -f "$TMPDIR/concat_hours.nc"

ncatted -O \
  -a units,time,o,c,"days since ${REF_DATE} 00:00:00" \
  -a calendar,time,o,c,"standard" \
  -a standard_name,time,o,c,"time" \
  -a axis,time,o,c,"T" \
  "$TMPDIR/axis_fixed.nc"

# ------------------------------------------------------------------------------
# 5) Rename variable and cast types to match MODIS-style choices
# ------------------------------------------------------------------------------
# Rename Band1 -> fpar (if needed)
ncrename -O -v Band1,fpar "$TMPDIR/axis_fixed.nc" 2>/dev/null || true

# Cast fpar/lon/lat to float IN PLACE (no renaming of dimscales -> avoids HDF5 error)
ncap2 -O -s 'fpar=float(fpar); lon=float(lon); lat=float(lat);' \
  "$TMPDIR/axis_fixed.nc" "$TMPDIR/typed.nc"
rm -f "$TMPDIR/axis_fixed.nc"

# ------------------------------------------------------------------------------
# 6) Metadata (minimal but clean)
# ------------------------------------------------------------------------------
# fpar metadata + fill
ncatted -O \
  -a long_name,fpar,o,c,"Fraction of Photosynthetically Active Radiation (FPAR) (masked)" \
  -a units,fpar,o,c,"1" \
  -a _FillValue,fpar,o,f,-9999.0 \
  -a missing_value,fpar,o,f,-9999.0 \
  "$TMPDIR/typed.nc"

# lon/lat metadata (MODIS-like text)
ncatted -O \
  -a standard_name,lon,o,c,"longitude" \
  -a long_name,lon,o,c,"longitude coordinate" \
  -a units,lon,o,c,"degrees_east" \
  -a axis,lon,o,c,"X" \
  -a standard_name,lat,o,c,"latitude" \
  -a long_name,lat,o,c,"latitude coordinate" \
  -a units,lat,o,c,"degrees_north" \
  -a axis,lat,o,c,"Y" \
  "$TMPDIR/typed.nc"

# ------------------------------------------------------------------------------
# 7) Write final + quick checks
# ------------------------------------------------------------------------------
mv -f "$TMPDIR/typed.nc" "$OUT"
rm -r _tmp_fpar_nc
echo "Wrote: $OUT"
ncdump -h "$OUT" | egrep 'float lon|float lat|float fpar|double time'
ncdump -v time "$OUT" | sed -n '/data:/,$p' | head -n 5
