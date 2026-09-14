
# R Processing and Analysis Pipeline for Natural LAI / FPAR

This folder contains the main R scripts for building and analyzing the natural-vegetation LAI/FPAR products.

## Workflow summary

1. Setup and reference-grid creation.
2. Georeferencing of raw LAI/FPAR inputs.
3. Land-cover preprocessing (CCI and GLC).
4. Mask construction (used land plus non-vegetated filters).
5. Mask application and spatial aggregation.
6. Trend and diagnostic analysis.

All rasters are aligned to shared global grids (0.05 degree native, 0.25 degree analysis).

## Folder structure

```text
R/
├── 00_setup.R
├── 01_georef_0p05.R
├── 02_cci_frac_0p05.R
├── 03_cci_mask_0p05.R
├── 04_glc_native_to_0p05.R
├── 04_glc_stack_0p05.R
├── 05_glc_mask_0p05.R
├── 06_nonveg_static_from_cci_0p05.R
├── 07_apply_nonveg_only_0p05.R
├── 08_agg_nonveg_0p25.R
├── 09_luh_pasture_overlap_0p25.R
├── 10_apply_mask_0p05.R
├── 11_agg_0p25.R
├── 12_make_lc025_fractions.py
├── 12_chapter2_inputs_0p5.R
├── analysis/
└── helpers/
```

## Scripts

### Setup

- 00_setup.R: Builds reference grids and area layers, sets paths, and writes the exact `config/config_<RUN_TAG>.yml` file for the selected scenario.

### Georeferencing

- 01_georef_0p05.R: Converts LAI/FPAR NetCDF inputs into aligned 0.05 degree global rasters.

### Land-cover preprocessing

- 02_cci_frac_0p05.R: Builds fractional cover layers from ESA-CCI/C3S.
- 04_glc_native_to_0p05.R: Aggregates the native GLC_FCS30D v2 tiles to categorical mode and fractional grass cover on the exact 0.05° grid. 
- 04_glc_stack_0p05.R: Harmonizes and stacks GLC_FCS30D maps on the project grid.
- 12_make_lc025_fractions.py: Generates annual 0.25° land-cover fractions from ESACCI classes (1992–2022). 

### Mask construction

- 03_cci_mask_0p05.R: Creates CCI-based used-land masks.
- 05_glc_mask_0p05.R: Creates GLC-based persistence masks from all 26 available maps (1985, 1990, 1995, and annually from 2000 to 2022).
- 06_nonveg_static_from_cci_0p05.R: Builds the 2007 static non-vegetated mask and the 1995/2022 snapshot-sensitivity masks.
- 09_luh_pasture_overlap_0p25.R: Adds LUH2 pasture-overlap diagnostics.

### Masking and aggregation

- 07_apply_nonveg_only_0p05.R: Applies non-vegetated exclusions.
- 08_agg_nonveg_0p25.R: Area-weights the water/ice-only monthly baseline to 0.25 degree before annual diagnostics and trends are calculated.
- 10_apply_mask_0p05.R: Applies selected masks to monthly LAI/FPAR.
- 11_agg_0p25.R: Aggregates to 0.25 degree for analysis.

### Chapter 2 inputs

- 12_chapter2_inputs_0p5.R: Area-weights unmasked monthly fAPAR to 0.5 degree and builds a 0.5-degree mask combining the CCI alpha 0.1 and pasture masks. The binary mask uses 1=drop and excludes a coarse cell if any contributing 0.05-degree cell is excluded. A fractional exclusion layer is also retained.

Mask convention is consistent across scripts:

- 1 = drop
- 0 = keep
- NA = undefined

## Helpers

The helper scripts in helpers/ contain operations that are reused across scripts, such as raster I/O, area-weighted aggregation, bootstrap intervals, and plotting.

## Complete analysis

After the masked trend products are complete, run all analysis scripts in their required order with:

```bash
R/analysis/run_all_analysis.sh
```

After changing a GLC input or mask setting, rebuild the complete GLC branch and
all dependent analysis with:

```bash
R/run_glc_workflow.sh
```
