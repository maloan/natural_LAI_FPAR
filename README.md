# **Natural LAI/FPAR Processing Pipeline**

A reproducible workflow for generating **natural-vegetation LAI and fAPAR datasets** from satellite observations. The pipeline builds land-use masks (CCI, GLC-FCS30D, LUH2), applies them to monthly LAI/FPAR at 0.05°, and produces area-weighted aggregates at 0.25° for global analyses.

------------------------------------------------------------------------

## **Key Outputs**

-   Natural-only LAI and fAPAR time series (monthly)
-   Masked 0.05° fields (LAI/FPAR)
-   Aggregated 0.25° fields using explicit area weights
-   Global trend maps (slope per decade)
-   Difference maps, scatter diagnostics, latitudinal profiles
-   Mask sensitivity and agreement statistics
-   Quicklooks for QA at each stage

All products are stored under:

```         
output/<run_tag>/
```

------------------------------------------------------------------------

## **Repository Structure**

```         
R/          # Core processing scripts
config/     # Main config.yml with paths, grids, AOIs, variable limits
data-raw/   # Raw datasets (CCI, LUH2, GLC-FCS30D, MODIS LAI/FPAR)
data/       # Harmonized inputs (georeferenced LAI/FPAR, fractional cover)
output/     # Final masked fields and aggregated results
analysis/   # Trend computation, plotting scripts, reproducibility Makefiles
src/        # Reference grids, area rasters, coastline layers, AOIs
vignettes/  # Walkthrough documentation
```

------------------------------------------------------------------------

## **Quick Start**

### 1. Set the project root

``` r
Sys.setenv(SNU_LAI_FPAR_ROOT = "~/path/to/natural_LAI_FPAR")
```

### 2. Run the full pipeline

``` bash
cd R
make all-full
```

This executes mask construction, field masking, aggregation, and trend analysis.

### 3. Run individual steps

``` bash
Rscript R/04_glc_mask_0p05.R
Rscript R/07_apply_mask_0p05.R
Rscript R/08_agg_0p25.R
```

Environment variables override defaults:

``` bash
TAU_CCI=0.1 USED_N_YEARS=3 Rscript R/07_apply_mask_0p05.R
```

------------------------------------------------------------------------

## **Requirements**

-   R ≥ 4.2
-   Core packages: `terra`, `yaml`, `dplyr`, `ggplot2`, `stringr`, `parallel`, `glue`
-   Optional: Python tools for ESA-CCI downloads (`data-raw/ESACCI/`)
-   System libraries: GDAL, PROJ, NetCDF

------------------------------------------------------------------------

## **External Data**

Place raw datasets under:

```         
data-raw/
    ESACCI/
    GLC_FCS30D/
    LUH2_v2h/
    LAI/
    FPAR/
```

These files are **not tracked by Git**.

Paths to all external data are configured in:

```         
config/config.yml
```

------------------------------------------------------------------------

## **Reproducibility Principles**

-   Fixed reference grids at 0.05° and 0.25°
-   Categorical data → nearest-neighbour resampling
-   Continuous data → bilinear resampling
-   Explicit area-weighted aggregation (cell areas from `src/`)
-   Provenance metadata and manifest files per stage
-   Quicklook PNGs for masks, fractional cover, and aggregated products

------------------------------------------------------------------------

## **Troubleshooting**

### Missing GDAL / NetCDF

Install development headers. Example (Ubuntu):

``` bash
sudo apt install gdal-bin libgdal-dev libproj-dev libnetcdf-dev
```

### Path errors

Ensure `SNU_LAI_FPAR_ROOT` is set and that your `config.yml` matches your local directory layout.

### Missing raw files

Check that `data-raw/` contains the expected inputs and subfolders.

------------------------------------------------------------------------

## **Documentation**

A detailed workflow walkthrough is available in:

```         
vignettes/vignette.Rmd
```

It covers configuration, execution, AOIs, mask logic, aggregation, and exploratory analysis.
