# R Processing and Analysis Pipeline for Managed-Land-Masked LAI / FPAR

This folder contains the main R scripts for building and analyzing LAI and FPAR products after excluding mapped managed areas.

## Workflow

1. Create reference grids and configuration.
2. Georeference LAI and FPAR inputs.
3. Preprocess CCI and GLC land cover.
4. Construct land-use, non-vegetated, and pasture masks.
5. Apply masks and aggregate to the analysis grid.
6. Compute trends, summaries, and diagnostics.

All rasters are aligned to common global grids at 0.05° for processing and 0.25° for analysis.

## Main scripts

### Setup and georeferencing

- `00_setup.R` creates reference grids, area layers, and run configuration.
- `01_georef_0p05.R` converts LAI and FPAR inputs to the common 0.05° grid.

### Land-cover preprocessing

- `02_cci_frac_0p05.R` creates fractional CCI land-cover layers.
- `04_glc_native_to_0p05.R` aggregates native GLC_FCS30D data to the 0.05° grid.
- `04_glc_stack_0p05.R` harmonizes and stacks GLC maps.
- `12_make_lc025_fractions.R` creates annual 0.25° CCI land-cover fractions.

### Mask construction

- `03_cci_mask_0p05.R` creates CCI-based managed-land masks.
- `05_glc_mask_0p05.R` creates GLC-based persistence masks.
- `06_nonveg_static_from_cci_0p05.R` creates the static non-vegetated mask and snapshot-year sensitivity masks.
- `09_luh_pasture_overlap_0p25.R` creates LUH2 pasture-overlap masks.

### Masking and aggregation

- `07_apply_nonveg_only_0p05.R` applies the common water and permanent snow or ice exclusions.
- `08_agg_nonveg_0p25.R` aggregates the unmasked baseline to 0.25°.
- `10_apply_mask_0p05.R` applies the combined masks to monthly LAI and FPAR.
- `11_agg_0p25.R` aggregates masked products to 0.25°.

### Chapter 2 inputs

- `12_chapter2_inputs_0p5.R` creates the 0.5° FPAR and mask inputs used for Chapter 2.

## Mask convention

```text id="5w09dd"
1  = exclude
0  = retain
NA = undefined
```

## Analysis

Analysis scripts are stored under `R/analysis/`. After the masked trend products have been generated, run:

```bash id="tjtw2e"
R/analysis/run_all_analysis.sh
```

The analysis workflow uses the completed masked and unmasked trend products and produces the statistical summaries, tables, diagnostics, and figures.

## Helpers

Shared functions for raster I/O, area weighting, bootstrap confidence intervals, climate classification, and plotting are stored under `R/helpers/`.