#!/usr/bin/env bash
# Run the complete Chapter 1 analysis in dependency order.

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_dir"

export RUN_TAG=alpha_0.1
run_r() {
  echo
  echo "==> Rscript $*"
  Rscript "$@"
}

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
run_r R/12_make_lc025_fractions.R
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
