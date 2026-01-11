# **Natural LAI / fAPAR Processing Pipeline**

This repository implements a reproducible workflow for deriving **global
natural-vegetation LAI and fAPAR datasets** from satellite observations
by explicitly excluding anthropogenic land-use signals.

Independent land-use masks derived from **ESA-CCI** and **GLC-FCS30D**
are applied to monthly LAI/fAPAR fields at **0.05°** resolution. Masked
fields are subsequently aggregated to **0.25°** using explicit area
weighting, enabling consistent global and regional trend analyses.

------------------------------------------------------------------------

## **Scientific Motivation**

Satellite vegetation products integrate multiple signals, including:

1.  **Ecosystem responses to climate variability and rising CO₂**, and
2.  **Land-use change and management effects**, such as cropland
    expansion, irrigation, and urbanisation.

For attribution studies and model evaluation, these signals must be
disentangled. This workflow isolates the **natural-vegetation
component** of observed greening and browning by removing anthropogenic
land cover prior to aggregation and trend estimation.

------------------------------------------------------------------------

## **Key Outputs**

-   Monthly **natural-vegetation LAI and fAPAR** fields (0.05° and
    0.25°)
-   Binary land-use masks derived independently from CCI and GLC
-   Area-weighted global and zonal mean time series
-   Pixel-wise trend estimates (year-mean and year-maximum)
-   Masked vs. unmasked trend comparisons
-   Diagnostics of masked (anthropogenic) regions
-   Quality-control quicklooks at all major processing stages

All outputs are written to:

```         
output/<RUN_TAG>/
```

where each run tag corresponds to a specific masking configuration.

------------------------------------------------------------------------

## **Repository Structure**

```         
R/            Processing and analysis scripts
config/       Central configuration (paths, grids, AOIs, thresholds)
data-raw/     External datasets (not version controlled)
data/         Harmonised intermediate products
output/       Masked fields, aggregates, evaluation products
analysis/     Trend results and figures
src/          Reference grids, area rasters, AOI masks
vignettes/    Extended workflow documentation
```

------------------------------------------------------------------------

## **Pipeline Overview**

1.  **Setup** Definition of global reference grids (0.05°, 0.25°),
    cell-area rasters, and AOIs *(00_setup.R)*

2.  **Georeferencing** Alignment of LAI and fAPAR to the canonical 0.05°
    grid *(01_georef_0p05.R)*

3.  **Land-cover preprocessing**

    -   ESA-CCI: fractional cropland, urban, grass, and bare cover
    -   GLC-FCS30D: categorical land-cover stacks *(03_cci_frac_0p05.R,
        05_glc_stack_0p05.R)*

4.  **Mask construction**

    -   Multi-year persistence-based natural-vegetation masks
    -   Independent implementations for CCI and GLC
    -   Static abiotic exclusions (water, snow/ice)
        *(04_cci_mask_0p05.R, 06_glc_mask_0p05.R,
        07_abiotic_static_from_cci_0p05.R)*

5.  **Mask application** Monthly masking of LAI and fAPAR at 0.05°
    *(11_apply_mask_0p05.R)*

6.  **Aggregation** Explicit area-weighted aggregation to 0.25°
    *(12_agg_0p25.R)*

7.  **Analysis** Global, zonal, and pixel-level trend analyses;
    masked–unmasked comparisons; diagnostics of excluded regions
    *(Scripts 13–17)*

------------------------------------------------------------------------

## **Quick Start**

### 1. Set the project root

``` r
Sys.setenv(SNU_LAI_FPAR_ROOT = "~/path/to/natural_LAI_FPAR")
```

### 2. Run the full pipeline and analyses

``` bash
cd R
make full-analysis
```

This executes setup, georeferencing, mask construction, masking,
aggregation, and trend analysis.

### 3. Run selected stages

``` bash
make georef
make apply-all
make agg-all
make analysis
```

Mask source and thresholds can be overridden, for example:

``` bash
MASKS=CCI TAU_CCI=0.1 make full-pipeline
```

------------------------------------------------------------------------

## **Configuration**

All paths, grids, AOIs, class definitions, and thresholds are specified
in:

```         
config/config.yml
```

Core design principles include:

-   Fixed global grids at 0.05° and 0.25°
-   Nearest-neighbour resampling for categorical data
-   Bilinear resampling for continuous variables
-   Explicit area-weighted aggregation
-   Fully deterministic, script-driven execution

------------------------------------------------------------------------

## **Requirements**

-   **R ≥ 4.2**
-   Core packages: `terra`, `yaml`, `dplyr`, `ggplot2`, `stringr`,
    `glue`, `parallel`, `sf`
-   System libraries: **GDAL**, **PROJ**, **NetCDF**

Example (Ubuntu):

``` bash
sudo apt install gdal-bin libgdal-dev libproj-dev libnetcdf-dev
```

------------------------------------------------------------------------

## **External Data**

Raw datasets must be placed manually under:

```         
data-raw/
  ESACCI/
  GLC_FCS30D/
  LUH2_v2h/
  LAI/
  FPAR/
```

These files are not tracked by Git. All dataset paths are defined in
`config/config.yml`.

------------------------------------------------------------------------

## **Reproducibility and Quality Control**

-   Deterministic processing (no stochastic elements)
-   Run-tag–scoped outputs
-   Manifest files recording raster geometry and global summaries
-   Quicklook PNGs for visual inspection
-   Independent mask constructions to assess robustness

------------------------------------------------------------------------

## **Documentation**

A detailed workflow description is provided in:

```         
vignettes/vignette.Rmd
```

covering mask logic, AOIs, aggregation methodology, trend
interpretation, and sensitivity analyses.
