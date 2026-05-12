---
editor_options: 
  markdown: 
    wrap: 72
---
# R Processing and Analysis Pipeline for Natural LAI / FPAR

This folder contains the main R scripts for building and analyzing the
natural-vegetation LAI/FPAR products.

In short, this is where the pipeline is executed: setup, preprocessing,
masking, aggregation, and analysis.

The workflow is built to be reproducible and grid-consistent, with land-use
exclusion handled explicitly.

## Workflow summary

The pipeline follows six broad steps:

1. Setup and reference-grid creation.
2. Georeferencing of raw LAI/FPAR inputs.
3. Land-cover preprocessing (CCI and GLC).
4. Mask construction (used land plus non-vegetated filters).
5. Mask application and spatial aggregation.
6. Trend and diagnostic analysis.

All rasters are aligned to shared global grids (0.05 degree native,
0.25 degree analysis).

## Folder structure

```text
R/
├── 00_setup.R
├── 01_georef_0p05.R
├── 02_cci_frac_0p05.R
├── 03_cci_mask_0p05.R
├── 04_glc_stack_0p05.R
├── 05_glc_mask_0p05.R
├── 06_nonveg_static_from_cci_0p05.R
├── 07_apply_nonveg_only_0p05.R
├── 08_luh_use_masks.R
├── 09_luh_pasture_overlap_0p25.R
├── 10_apply_mask_0p05.R
├── 11_agg_0p25.R
├── 11_agg_0p5.R
├── 12_make_lc025_majority.R
├── Makefile
├── analysis/
└── helpers/
```

## What the main scripts do

### Setup

- 00_setup.R: Builds reference grids and area layers, sets paths, and writes
  config/config.yml for the selected run tag.

### Georeferencing

- 01_georef_0p05.R: Converts LAI/FPAR NetCDF inputs into aligned 0.05 degree
  global rasters.

### Land-cover preprocessing

- 02_cci_frac_0p05.R: Builds fractional cover layers from ESA-CCI/C3S.
- 04_glc_stack_0p05.R: Harmonizes and stacks GLC_FCS30D maps on the project
  grid.
- 12_make_lc025_majority.R: Generates 0.25° landcover majority maps from annual
  ESACCI data (1992-2020) for trend analysis by land-cover class. Auto-triggered
  by analysis scripts if needed.

### Mask construction

- 03_cci_mask_0p05.R: Creates CCI-based used-land masks.
- 05_glc_mask_0p05.R: Creates GLC-based persistence masks.
- 06_nonveg_static_from_cci_0p05.R: Builds static non-vegetated masks.
- 09_luh_pasture_overlap_0p25.R: Adds LUH2 pasture-overlap diagnostics.

### Masking and aggregation

- 07_apply_nonveg_only_0p05.R: Applies non-vegetated exclusions.
- 08_luh_use_masks.R: Processes LUH2 layers used by masking.
- 10_apply_mask_0p05.R: Applies selected masks to monthly LAI/FPAR.
- 11_agg_0p5.R: Aggregates to 0.5 degree.
- 11_agg_0p25.R: Aggregates to 0.25 degree for analysis.

Mask convention is consistent across scripts:

- 1 = drop
- 0 = keep
- NA = undefined

## Helpers

The helper scripts in helpers/ are shared utilities and are not intended to be
run directly.

They cover common tasks like file I/O, NetCDF handling, options parsing,
visualization helpers, and utility wrappers used throughout the pipeline.

## Makefile Usage

Run from this folder:

```bash
make analysis
```

This runs the full chain from setup through analysis.

Common targets:

```bash
make pipeline
make pipeline MASKS=CCI
make ql
```

Common variables:

- RUN_TAG: output namespace, for example tau_0.2.
- MASKS: mask source selection (CCI, GLC, or both).
