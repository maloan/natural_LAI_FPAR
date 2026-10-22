# Analysis

This folder contains post-processing results and diagnostics derived from the masked and unmasked LAI products. Most analyses use the 0.25° products. Existing FPAR trend files under `unmasked/` are retained, but the current Chapter 1 analysis runner summarizes LAI.

## Contents

The analysis includes:

- global and regional trend summaries
- masked and unmasked comparisons
- LAI diagnostics
- grid-cell significance results
- global and zonal time-series summaries

## Structure

- `unmasked/` contains the 0.25° baseline LAI and FPAR products used as the unmasked reference.
- `tmp/` contains temporary analysis files.
- `results/figures/` contains summary and diagnostic figures.
- `results/tables/` contains tabular outputs used for diagnostics and manuscript results.

Analysis products are derived from the processed outputs under `output/`.

The corresponding analysis scripts are stored in `R/analysis/`.
