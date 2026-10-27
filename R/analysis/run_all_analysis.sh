#!/usr/bin/env bash
# Run the complete Chapter 1 analysis in dependency order.

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_dir"

export RUN_TAG=alpha_0.1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export GDAL_NUM_THREADS="${GDAL_NUM_THREADS:-1}"
workflow_nice_level="${WORKFLOW_NICE_LEVEL:-15}"

run_r() {
  echo
  echo "==> Rscript $*"
  nice -n "$workflow_nice_level" ionice -c 3 Rscript "$@"
}

# Refresh the water/ice-only LAI baseline when its inputs or code have changed.
unmasked_lai="analysis/unmasked/0p25/LAI_georef_yearmean_trend_relative_percent_peryear_0p25.nc"
georef_lai_dir="data/georef/georef_lai_0p05"
georef_lai_count=$(find "$georef_lai_dir" -maxdepth 1 -name 'LAI_*_0p05.tif' | wc -l)
[[ "$georef_lai_count" -eq 516 ]] || {
  echo "Expected 516 georeferenced LAI months, found $georef_lai_count" >&2
  exit 1
}
baseline_sources=(
  R/07_apply_nonveg_only_0p05.R
  R/08_agg_nonveg_0p25.R
  trends/01_build_unmasked_0p25.sh
  trends/compute_mk_pval.R
  config/config_alpha_0.1.yml
  output/alpha_0.1/masks/mask_nonvegetated/mask_nonvegetated_CCI_2007_alphaW0p05_alphaI0p05_0p05.tif
)
refresh_unmasked=0
if [[ ! -s "$unmasked_lai" ]]; then
  refresh_unmasked=1
else
  for source in "${baseline_sources[@]}"; do
    [[ -f "$source" ]] || { echo "Missing baseline dependency: $source" >&2; exit 1; }
    if [[ "$source" -nt "$unmasked_lai" ]]; then
      refresh_unmasked=1
    fi
  done
  if [[ -n "$(find "$georef_lai_dir" -maxdepth 1 -name 'LAI_*_0p05.tif' -newer "$unmasked_lai" -print -quit)" ]]; then
    refresh_unmasked=1
  fi
fi

if [[ "$refresh_unmasked" -eq 1 ]]; then
  VAR=LAI run_r R/07_apply_nonveg_only_0p05.R
  VAR=LAI run_r R/08_agg_nonveg_0p25.R
  VAR=LAI nice -n "$workflow_nice_level" ionice -c 3 \
    bash trends/01_build_unmasked_0p25.sh
fi

# Domain and mask diagnostics. A path-only migration can use the saved,
# previously validated CCI-derived nonvegetated mask when raw CCI is absent.
if [[ "${USE_SAVED_NONVEG:-0}" == "1" ]]; then
  nonveg_mask="output/alpha_0.1/masks/mask_nonvegetated/mask_nonvegetated_CCI_2007_alphaW0p05_alphaI0p05_0p05.tif"
  [[ -s "$nonveg_mask" ]] || {
    echo "Missing saved nonvegetated mask: $nonveg_mask" >&2
    exit 1
  }
  echo "Using saved nonvegetated mask: $nonveg_mask"
else
  run_r R/06_nonveg_static_from_cci_0p05.R
fi
run_r R/analysis/00_area_validdomain_after_nonvegetated.R
run_r R/analysis/01_masking_footprint_summary.R

# Global and zonal trend summaries.
run_r R/analysis/02_global_trends_summary.R
run_r R/analysis/03_matched_period_annual_mean_trends.R
run_r R/analysis/04_global_absolute_trends_timeseries.R
run_r R/analysis/05_zonal_seasonal_amplitude.R
run_r R/analysis/06_zonal_absolute_trends_all_masks.R
run_r R/analysis/07_nonveg_snapshot_sensitivity.R

# Maps and combined diagnostics.
run_r R/analysis/08_lai_yearmean_trend_maps.R
run_r R/analysis/09_lai_yearmax_trend_maps.R
run_r R/analysis/10_cropland_pasture_trend_diagnostics.R
run_r R/analysis/11_zonal_diagnostics_combined_figure.R

# Land-cover summaries use a dominant 0.25-degree class and the fixed
# post-nonvegetated support area.
echo
run_r R/tests/test_landcover_classification.R
if [[ "${USE_SAVED_LC025:-0}" == "1" ]]; then
  lc_fraction_dir="analysis/tmp/lc025_fraction_yearly_primary22_v1"
  lc_fraction_count=$(find "$lc_fraction_dir" -maxdepth 1 \
    -type f -name 'lc025_fraction_*.tif' -size +0c | wc -l)
  [[ "$lc_fraction_count" -eq 31 ]] || {
    echo "Expected 31 saved primary22_v1 land-cover fractions, found $lc_fraction_count" >&2
    exit 1
  }
  echo "Using 31 saved primary22_v1 annual land-cover fractions: $lc_fraction_dir"
else
  run_r R/12_make_lc025_fractions.R
fi
run_r R/analysis/12_0_landcover_dominant_class.R
run_r R/analysis/12_a_landcover_trend_summary.R use_relative=false
run_r R/analysis/12_a_landcover_trend_summary.R use_relative=true

# Produce absolute and relative climate-zone summaries from area-dominant
# aggregation of the nominal 100-arc-second Köppen-Geiger classification.
run_r R/analysis/13_kg_trend_summary.R use_relative=false
run_r R/analysis/13_kg_trend_summary.R use_relative=true
run_r R/analysis/14_climate_landcover_trend_figure.R

echo
echo "All Chapter 1 analysis scripts completed."
