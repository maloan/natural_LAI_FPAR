# Analysis

This folder contains the post-processing results and diagnostics built from the masked and unmasked LAI/FPAR products.
Most analyses are based on the 0.25 degree outputs and are designed to be reproducible.

## What this folder is for

The analysis layer covers:

- Global and regional trend estimates.
- Masked vs unmasked comparisons.
- LAI and FPAR consistency checks.
- Greening and browning summaries.
- Pixel-level significance testing.
- Global and zonal time-series summaries.

## Main subfolders

### Baseline data

- unmasked/: Unmasked LAI/FPAR fields used as reference baselines.
- tmp/: Temporary and intermediate analysis artifacts.

### Current results (results)

- global_mean_relative_trends/: Global mean trend and relative-change summaries.
- kg_trends/: Trend summaries by Koppen-Geiger climate zone.
- lc_trends/: Trend summaries by land-cover class.
- mask_sensitivity_tau/: Sensitivity results across masking thresholds.
- masks/: Processed mask combinations used in analysis.
- metrics_mask_effect_zonal/: Zonal metrics of masking effects.
- paper_figures/: Figure-ready outputs for manuscripts.
- relative_trend_distribution/: Spatial distribution summaries of relative trends.
- spatial_redistribution_after_masking/: How masking shifts spatial trend patterns.
- zonal_relative_trends/: Zonal profiles of relative trend change.
- dropped_region/: Diagnostics for regions removed by masking.

## Notes

- These results are derived from output/<run_tag>/masked_0p25 products.
- Naming usually encodes variable, metric, mask setup, and run tag.
- This folder stores results and diagnostics only. Processing scripts live in R.
