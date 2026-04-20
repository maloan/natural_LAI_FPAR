# ==============================================================================
# build_trends_masked_0p25.sh — Run masked 0.25° trend workflow from repo root
#
# Examples
#   ./build_trends_masked_0p25.sh tau_0.2 FPAR GLC
#   ./build_trends_masked_0p25.sh tau_0.1 LAI  CCI
#
# ==============================================================================
set -euo pipefail

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need_cmd gdal_translate
need_cmd cdo

# ------------------------------------------------------------------------------
# Args
# ------------------------------------------------------------------------------
TAU=${1:?Provide TAU folder, e.g. tau_0.2}
VAR=${2:?Provide VAR: LAI or FPAR}
MASKTAG=${3:?Provide MASKTAG: CCI or GLC}
RES="0p25"

# ------------------------------------------------------------------------------
# Repo root + derived paths
# ------------------------------------------------------------------------------
ROOT="${SNU_LAI_FPAR_ROOT:-$HOME/GitHub/natural_LAI_FPAR}"
IN_DIR="${ROOT}/output/${TAU}/masked_0p25/${VAR}/masked_${VAR}_${MASKTAG}"
OUT_EVAL="${ROOT}/output/${TAU}/eval/trend_${VAR}_${MASKTAG}"

# ------------------------------------------------------------------------------
# Thresholds for relative trends and significance
# ------------------------------------------------------------------------------
EPS_LAI="${EPS_LAI:-0.05}"
EPS_FPAR="${EPS_FPAR:-0.02}"
ALPHA="${ALPHA:-0.05}"

case "$VAR" in
  LAI)  EPS="$EPS_LAI" ;;
  FPAR) EPS="$EPS_FPAR" ;;
  *) echo "VAR must be LAI or FPAR, got: $VAR" >&2; exit 1 ;;
esac

# ------------------------------------------------------------------------------
# Check folders and move into input directory
# ------------------------------------------------------------------------------
if [[ ! -d "$IN_DIR" ]]; then
  echo "Input directory not found:" >&2
  echo "  $IN_DIR" >&2
  echo "Check TAU/VAR/MASKTAG or your folder naming." >&2
  exit 1
fi

mkdir -p "$OUT_EVAL"
cd "$IN_DIR"

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------
PATTERN="${VAR}_masked_*.tif"          # expects: VAR_masked_YYYYMM_0p25.tif
MONTHLY="${VAR}_masked_monthly_${RES}.nc"

# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------
ANNUAL_FILES=(
  "${VAR}_yearmean_${RES}.nc"
  "${VAR}_yearmax_${RES}.nc"
  "${VAR}_yearmin_${RES}.nc"
  "${VAR}_yearamp_${RES}.nc"
)

TREND_FILES=(
  "${VAR}_yearmean_trend_intercept_${RES}.nc" "${VAR}_yearmean_trend_slope_peryear_${RES}.nc"
  "${VAR}_yearmax_trend_intercept_${RES}.nc"  "${VAR}_yearmax_trend_slope_peryear_${RES}.nc"
  "${VAR}_yearmin_trend_intercept_${RES}.nc"  "${VAR}_yearmin_trend_slope_peryear_${RES}.nc"
  "${VAR}_yearamp_trend_intercept_${RES}.nc"  "${VAR}_yearamp_trend_slope_peryear_${RES}.nc"
)

REL_FILES=(
  "${VAR}_yearmean_trend_relative_peryear_${RES}.nc"
  "${VAR}_yearmax_trend_relative_peryear_${RES}.nc"
  "${VAR}_yearmin_trend_relative_peryear_${RES}.nc"
  "${VAR}_yearamp_trend_relative_peryear_${RES}.nc"
)

MK_FILES=(
  "${VAR}_yearmean_mk_sig_pLT${ALPHA}_${RES}.nc"
  "${VAR}_yearmax_mk_sig_pLT${ALPHA}_${RES}.nc"
  "${VAR}_yearmin_mk_sig_pLT${ALPHA}_${RES}.nc"
  "${VAR}_yearamp_mk_sig_pLT${ALPHA}_${RES}.nc"
)

REL_SIG_FILES=(
  "${VAR}_yearmean_trend_relative_peryear_sig_pLT${ALPHA}_${RES}.nc"
  "${VAR}_yearmax_trend_relative_peryear_sig_pLT${ALPHA}_${RES}.nc"
  "${VAR}_yearmin_trend_relative_peryear_sig_pLT${ALPHA}_${RES}.nc"
  "${VAR}_yearamp_trend_relative_peryear_sig_pLT${ALPHA}_${RES}.nc"
)

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
rm -f "$MONTHLY" "${ANNUAL_FILES[@]}" "${TREND_FILES[@]}" "${REL_FILES[@]}" "${MK_FILES[@]}" "${REL_SIG_FILES[@]}"
rm -rf tmp_nc
mkdir -p tmp_nc

# ------------------------------------------------------------------------------
# 1) GeoTIFF to dated NetCDF
# ------------------------------------------------------------------------------
n_tif=$(ls -1 ${PATTERN} 2>/dev/null | wc -l | awk '{print $1}')
if [[ "$n_tif" -eq 0 ]]; then
  echo "No files match: ${PATTERN}" >&2
  echo "In directory: $IN_DIR" >&2
  exit 1
fi

echo "============================================================"
echo "ROOT    = ${ROOT}"
echo "TAU     = ${TAU}"
echo "VAR     = ${VAR}"
echo "MASKTAG = ${MASKTAG}"
echo "IN_DIR  = ${IN_DIR}"
echo "OUT_EVAL= ${OUT_EVAL}"
echo "RES     = ${RES}"
echo "EPS     = ${EPS}"
echo "ALPHA   = ${ALPHA}"
echo "N_TIF   = ${n_tif}"
echo "============================================================"

echo "Step 1: Converting GeoTIFFs to dated NetCDF (tmp_nc/)"
for tif in ${PATTERN}; do
  base=$(basename "$tif")
  yyyymm=${base#${VAR}_masked_}
  yyyymm=${yyyymm%_${RES}.tif}

  year=${yyyymm:0:4}
  mon=${yyyymm:4:2}

  raw_nc="tmp_nc/${VAR}_${yyyymm}.nc"
  dated_nc="tmp_nc/${VAR}_${yyyymm}_dated.nc"

  echo "  ${tif} -> ${dated_nc}"
  gdal_translate -of NETCDF "$tif" "$raw_nc" >/dev/null
  cdo -O setdate,${year}-${mon}-01 "$raw_nc" "$dated_nc" >/dev/null
done

echo "Step 2: Merging monthly time series -> ${MONTHLY}"
cdo -O mergetime tmp_nc/${VAR}_*_dated.nc "${MONTHLY}" >/dev/null
rm -rf tmp_nc

# ------------------------------------------------------------------------------
# 2) Annual metrics
# ------------------------------------------------------------------------------
echo "Step 3: Computing annual metrics"
cdo -O yearmean "${MONTHLY}" "${VAR}_yearmean_${RES}.nc" >/dev/null
cdo -O yearmax  "${MONTHLY}" "${VAR}_yearmax_${RES}.nc"  >/dev/null
cdo -O yearmin  "${MONTHLY}" "${VAR}_yearmin_${RES}.nc"  >/dev/null
cdo -O sub "${VAR}_yearmax_${RES}.nc" "${VAR}_yearmin_${RES}.nc" "${VAR}_yearamp_${RES}.nc" >/dev/null

# ------------------------------------------------------------------------------
# 3) OLS trends on annual metrics
# ------------------------------------------------------------------------------
echo "Step 4: Computing OLS trends (cdo trend)"
for MET in yearmean yearmax yearmin yearamp; do
  cdo -O trend "${VAR}_${MET}_${RES}.nc" \
    "${VAR}_${MET}_trend_intercept_${RES}.nc" \
    "${VAR}_${MET}_trend_slope_peryear_${RES}.nc" >/dev/null
done

# ------------------------------------------------------------------------------
# 4) Relative trends
# ------------------------------------------------------------------------------
echo "Step 5: Computing relative trends"
for MET in yearmean yearmax yearmin yearamp; do
  ts="${VAR}_${MET}_${RES}.nc"
  sl="${VAR}_${MET}_trend_slope_peryear_${RES}.nc"

  mu="tmp_${VAR}_${MET}_meanlevel_${RES}.nc"
  mueps="tmp_${VAR}_${MET}_meanlevel_ge${EPS}_mask_${RES}.nc"
  out="${VAR}_${MET}_trend_relative_peryear_${RES}.nc"

  # mean level
  cdo -O timmean "$ts" "$mu" >/dev/null

  # EPS mask as a single-layer 0/1 field
  cdo -O gec,"${EPS}" "$mu" "$mueps" >/dev/null

  # relative trend (defined only where meanlevel >= EPS)
  cdo -O setname,trend_relative \
    -ifthen "$mueps" \
    -div "$sl" "$mu" \
    "$out" >/dev/null
done

# ------------------------------------------------------------------------------
# 5) Move outputs to output/<TAU>/eval/
# ------------------------------------------------------------------------------
mv "$MONTHLY" "${ANNUAL_FILES[@]}" "${TREND_FILES[@]}" "${REL_FILES[@]}" "${MK_FILES[@]}" "${REL_SIG_FILES[@]}" "$OUT_EVAL/"

