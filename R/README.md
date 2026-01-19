# `R/` — Processing and Analysis Pipeline for Natural LAI / fAPAR

This directory contains the **core processing, masking, aggregation, and
analysis scripts** used to construct and analyse global
**natural-vegetation LAI and fAPAR datasets** from satellite
observations. Execution is orchestrated via the accompanying
**Makefile**, which defines the recommended order of operations.

The workflow is designed to be **reproducible**, **grid-consistent**,
and **explicit about land-use exclusion**.

------------------------------------------------------------------------

## Conceptual overview

The pipeline proceeds in six logical stages:

1.  **Project setup & reference grids**
2.  **Georeferencing of raw LAI / fAPAR**
3.  **Land-cover preprocessing (CCI, GLC)**
4.  **Mask construction (cropland, urban, pasture, abiotic)**
5.  **Mask application and spatial aggregation**
6.  **Trend analysis and diagnostics**

All rasters are aligned to canonical global grids (0.05° native, 0.25°
analysis).

------------------------------------------------------------------------

## Directory structure (this folder)

```         
R/
├── 00_setup.R
├── 01_georef_0p05.R
├── 02_apply_abiotic_only_0p05.R
├── 03_cci_frac_0p05.R
├── 04_cci_mask_0p05.R
├── 05_glc_stack_0p05.R
├── 06_glc_mask_0p05.R
├── 07_abiotic_static_from_cci_0p05.R
├── 10_luh_pasture_overlap_025.R
├── 11_apply_mask_0p05.R
├── 12_agg_0p25.R
├── 13_analyse_georef_0p25.R
├── 14_compare_lai_fpar_trends.R
├── 15_analyse_masked_trends.R
├── 16_compare_masked_unmasked.R
├── 17_analyse_trend_mask.R
├── Makefile
│
├── helpers /
│   ├── utils.R
│   ├── io.R
│   ├── geom.R
│   ├── viz.R
│   └── options.R
```

------------------------------------------------------------------------

## Script roles (high level)

### 0. Setup

-   **`00_setup.R`** Initializes reference grids (0.05°, 0.25°), area
    rasters, AOIs, and writes `config.yml`. Must be run **once per run
    tag**.

------------------------------------------------------------------------

### 1. Georeferencing

-   **`01_georef_0p05.R`** Converts monthly LAI / fAPAR NetCDFs to
    aligned 0.05° global rasters.

------------------------------------------------------------------------

### 2. Land-cover preprocessing

-   **`03_cci_frac_0p05.R`** Builds fractional land-cover fields
    (cropland, urban, grass, bare) from ESA-CCI.
-   **`05_glc_stack_0p05.R`** Stacks yearly GLC-FCS30D categorical maps
    onto the 0.05° grid.

------------------------------------------------------------------------

### 3. Mask construction

-   **`04_cci_mask_0p05.R`** Natural-vegetation mask from ESA-CCI
    time-window logic.
-   **`06_glc_mask_0p05.R`** Equivalent mask derived from GLC-FCS30D.
-   **`07_abiotic_static_from_cci_0p05.R`** Static water / snow / ice
    mask.
-   **`10_luh_pasture_overlap_025.R`** LUH2 pasture overlap diagnostics
    at 0.25°.

Masks follow the convention: **1 = drop**, **0 = keep**, **NA =
undefined**.

------------------------------------------------------------------------

### 4. Mask application & aggregation

-   **`11_apply_mask_0p05.R`** Applies selected masks to monthly LAI /
    fAPAR at 0.05°.
-   **`12_agg_0p25.R`** Area-weighted aggregation to 0.25° for analysis.

------------------------------------------------------------------------

### 5. Analysis scripts

-   **`13_analyse_georef_0p25.R`** Baseline global and zonal trends
    (unmasked).
-   **`14_compare_lai_fpar_trends.R`** Structural comparison of LAI vs
    fAPAR trends.
-   **`15_analyse_masked_trends.R`** Trends for masked (natural-only)
    datasets.
-   **`16_compare_masked_unmasked.R`** Direct comparison of masked vs
    unmasked trends.
-   **`17_analyse_trend_mask.R`** Analysis of *removed* (masked-out)
    regions.

------------------------------------------------------------------------

## Helper modules

The `helpers/` scripts are **not meant to be run directly**.

-   **`utils.R`** – general helpers (timing, config, parallelism)
-   **`io.R`** – raster I/O, GDAL options, provenance
-   **`geom.R`** – grid alignment, NetCDF handling, aggregation
-   **`viz.R`** – consistent plotting utilities
-   **`options.R`** – environment-variable parsing

------------------------------------------------------------------------

## Makefile usage (recommended)

From this directory:

``` bash
make analysis
```

This runs the **full pipeline**:

-   setup → preprocessing → masks → masking → aggregation → analyses

Common variants:

``` bash
make pipeline            # data prep only
make pipeline MASKS=CCI  # ESA-CCI mask only
make ql                  # regenerate quicklooks
```

Variables:

-   `RUN_TAG` — output namespace (e.g. `tau_0.2`)
-   `MASKS` — which masks to run (`CCI`, `GLC`, or both)
