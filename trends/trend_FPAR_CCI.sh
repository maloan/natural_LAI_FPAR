#!/usr/bin/env bash
set -euo pipefail

# --- config ---
VAR="FPAR"
RES="0p25"
PATTERN="${VAR}_masked_*.tif"         # e.g., FPAR_masked_198201_0p25.tif
OUT_MONTHLY="${VAR}_masked_monthly_${RES}.nc"

# thresholds for relative trend (paper definition)
EPS_FPAR="${EPS_FPAR:-0.02}"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need_cmd gdal_translate
need_cmd cdo

# --- 0) build monthly NetCDF with time axis ---
rm -rf tmp_nc
mkdir -p tmp_nc

n_tif=$(ls -1 ${PATTERN} 2>/dev/null | wc -l | awk '{print $1}')
if [[ "$n_tif" -eq 0 ]]; then
  echo "No files match: ${PATTERN}"
  exit 1
fi

echo "Converting ${n_tif} GeoTIFFs -> NetCDF (temporary)..."

for tif in ${PATTERN}; do
  base=$(basename "$tif")
  yyyymm=${base#${VAR}_masked_}
  yyyymm=${yyyymm%_${RES}.tif}

  year=${yyyymm:0:4}
  mon=${yyyymm:4:2}

  raw_nc="tmp_nc/${VAR}_${yyyymm}.nc"
  dated_nc="tmp_nc/${VAR}_${yyyymm}_dated.nc"

  echo "  ${tif} -> ${dated_nc}"
  gdal_translate -of NETCDF "$tif" "$raw_nc"
  cdo -O setdate,${year}-${mon}-01 "$raw_nc" "$dated_nc"
done

echo "Merging monthly time series -> ${OUT_MONTHLY}"
cdo -O mergetime tmp_nc/${VAR}_*_dated.nc "${OUT_MONTHLY}"
rm -rf tmp_nc

# --- 1) annual metrics ---
echo "Computing annual metrics..."
cdo -O yearmean "${OUT_MONTHLY}" "${VAR}_yearmean_${RES}.nc"
cdo -O yearmax  "${OUT_MONTHLY}" "${VAR}_yearmax_${RES}.nc"
cdo -O yearmin  "${OUT_MONTHLY}" "${VAR}_yearmin_${RES}.nc"
cdo -O sub "${VAR}_yearmax_${RES}.nc" "${VAR}_yearmin_${RES}.nc" "${VAR}_yearamp_${RES}.nc"

# --- 2) trends (OLS slope per year + intercept) ---
echo "Computing trends..."
for MET in yearmean yearmax yearmin yearamp; do
  cdo -O trend "${VAR}_${MET}_${RES}.nc" \
    "${VAR}_${MET}_trend_intercept_${RES}.nc" \
    "${VAR}_${MET}_trend_slope_peryear_${RES}.nc"
done

# --- 3) relative trends (recommended; matches your Methods) ---
echo "Computing relative trends (slope / local mean level)..."
for MET in yearmean yearmax yearmin yearamp; do
  ts="${VAR}_${MET}_${RES}.nc"
  sl="${VAR}_${MET}_trend_slope_peryear_${RES}.nc"
  mu="${VAR}_${MET}_meanlevel_${RES}.nc"
  mask="${VAR}_${MET}_meanlevel_geEPS_${RES}.nc"
  out="${VAR}_${MET}_trend_relative_peryear_${RES}.nc"

  # mean level over time (per pixel)
  cdo -O timmean "$ts" "$mu"
  # validity mask (avoid unstable ratios)
  cdo -O gec,"${EPS_FPAR}" "$mu" "$mask"
  # relative trend
  cdo -O setname,trend_relative -ifthen "$mask" -div "$sl" "$mu" "$out"

  rm -f "$mu" "$mask"
done

echo "Done."
