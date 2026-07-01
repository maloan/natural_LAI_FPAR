# Analysis

This folder contains the post-processing results and diagnostics built from the masked and unmasked LAI/FPAR products. Most analyses are based on the 0.25 degree outputs.

## What this folder is for

The analysis layer covers:

- Global and regional trend estimates.
- Masked vs unmasked comparisons.
- LAI and FPAR consistency checks.
- Pixel-level significance testing.
- Global and zonal time-series summaries.

## Main subfolders

### Baseline data

- unmasked/: Unmasked 0.25° LAI/FPAR products used as reference baselines.
- tmp/: Temporary and intermediate analysis artifacts.

### Current results

- figures/: Summary and diagnostic figures.
- tables/: CSV outputs used for diagnostics and manuscript assets.

## Notes

- These results are derived from output/<run_tag>/masked_0p25 products.
- This folder stores results and diagnostics only. Processing scripts live in R/analysis.
