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

# Build the water/ice-only baseline once when it is absent.
unmasked_lai="analysis/unmasked/0p25/LAI_georef_yearmean_trend_relative_percent_peryear_0p25.nc"
unmasked_fpar="analysis/unmasked/0p25/FPAR_georef_yearmean_trend_relative_percent_peryear_0p25.nc"

if [[ ! -f "$unmasked_lai" || ! -f "$unmasked_fpar" ]]; then
  VAR=LAI run_r R/07_apply_nonveg_only_0p05.R
  VAR=FPAR run_r R/07_apply_nonveg_only_0p05.R
  run_r R/08_agg_nonveg_0p25.R
  bash trends/01_build_unmasked_0p25.sh
fi

# Domain and mask diagnostics.
run_r R/06_nonveg_static_from_cci_0p05.R
run_r R/analysis/00_area_validdomain_after_nonvegetated.R
run_r R/analysis/01_masking_footprint_summary.R

# Global and zonal trend summaries.
run_r R/analysis/02_global_relative_trends_summary.R
run_r R/analysis/03_global_absolute_trends_summary.R
run_r R/analysis/04_global_absolute_trends_timeseries.R
run_r R/analysis/05_zonal_seasonal_amplitude.R
run_r R/analysis/06_zonal_absolute_trends_all_masks.R
run_r R/analysis/07_nonveg_snapshot_sensitivity.R

# Maps and combined diagnostics.
run_r R/analysis/08_lai_yearmean_trend_maps.R
run_r R/analysis/09_lai_yearmax_trend_maps.R
run_r R/analysis/10_cropland_pasture_trend_diagnostics.R
run_r R/analysis/11_zonal_diagnostics_combined_figure.R

# Land-cover summaries require both absolute and relative variants.
echo
echo "==> python3 R/12_make_lc025_fractions.py"
python3 R/12_make_lc025_fractions.py
run_r R/analysis/12_a_landcover_trend_summary.R use_relative=false
run_r R/analysis/12_a_landcover_trend_summary.R use_relative=true
run_r R/analysis/12_b_landcover_abs_vs_rel_trend.R

# Produce both absolute and relative climate-zone summaries.
run_r R/analysis/13_kg_trend_summary.R use_relative=false
run_r R/analysis/13_kg_trend_summary.R use_relative=true

echo
echo "All Chapter 1 analysis scripts completed."
