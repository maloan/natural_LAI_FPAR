#!/usr/bin/env bash
# Build the data products and analysis results.

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$repo_dir"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export GDAL_NUM_THREADS=1

run() {
  echo
  echo "==> $*"
  nice -n 15 ionice -c 3 "$@"
}

scenarios=(alpha_0.05 alpha_0.1 alpha_0.2 glc)
variables=(LAI FPAR)

# Configuration and shared grids.
for run_tag in "${scenarios[@]}"; do
  RUN_TAG="$run_tag" run Rscript R/00_setup.R
done

# Shared monthly vegetation inputs.
for variable in "${variables[@]}"; do
  RUN_TAG=alpha_0.1 VAR="$variable" run Rscript R/01_georef_0p05.R
done

# Land-cover fractions and static masks.
RUN_TAG=alpha_0.1 run Rscript R/02_cci_frac_0p05.R
for run_tag in alpha_0.05 alpha_0.1 alpha_0.2; do
  RUN_TAG="$run_tag" run Rscript R/03_cci_mask_0p05.R
done
run Rscript R/04_glc_native_to_0p05.R
RUN_TAG=glc run Rscript R/04_glc_stack_0p05.R
RUN_TAG=glc run Rscript R/05_glc_mask_0p05.R
RUN_TAG=alpha_0.1 run Rscript R/06_nonveg_static_from_cci_0p05.R

for variable in "${variables[@]}"; do
  RUN_TAG=alpha_0.1 VAR="$variable" run Rscript R/07_apply_nonveg_only_0p05.R
done
RUN_TAG=alpha_0.1 run Rscript R/08_agg_nonveg_0p25.R
run bash trends/01_build_unmasked_0p25.sh

# LUH2 pasture masks, monthly masking, and 0.25-degree aggregation.
for run_tag in "${scenarios[@]}"; do
  RUN_TAG="$run_tag" run Rscript R/09_luh_pasture_overlap_0p25.R
  for variable in "${variables[@]}"; do
    RUN_TAG="$run_tag" VAR="$variable" run Rscript R/10_apply_mask_0p05.R
    RUN_TAG="$run_tag" VAR="$variable" run Rscript R/11_agg_0p25.R
  done
done

# Trend products, Chapter 2 inputs, and Chapter 1 results.
run bash trends/02_batch_build_trends_masked.sh
RUN_TAG=alpha_0.1 run Rscript R/11_agg_0p5.R
RUN_TAG=alpha_0.1 run Rscript R/12_chapter2_inputs_0p5.R
run bash R/analysis/run_all_analysis.sh

echo
echo "Data and analysis workflow finished successfully."
