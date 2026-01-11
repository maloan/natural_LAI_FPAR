# **Configuration (`config/`)**

This directory contains the **central configuration** for the
natural-vegetation LAI/FPAR workflow. The file `config.yml` is read by
all R scripts via `cfg_read()` and defines the project domain, grids,
paths, and dataset metadata.

------------------------------------------------------------------------

## **Contents**

-   **`config.yml`**\
    Single configuration file specifying:
    -   project metadata (run tag, CRS, time coverage)
    -   input and output paths
    -   canonical 0.05° and 0.25° grids (reference, area, AOI rasters)
    -   AOI definitions for quicklooks and regional analyses
    -   land-cover class mappings (ESA-CCI, GLC-FCS30D, LUH2)
    -   standard file-naming templates
-   **`README.md`**\
    Brief description of how configuration is used by the pipeline.

------------------------------------------------------------------------

## **Typical adjustments**

Most changes are limited to:

-   **Paths**: point to local or HPC data locations.

-   **AOIs**: add or modify regional bounding boxes.

-   **Time ranges**: adjust LAI/FPAR, CCI, or GLC year windows.

-   **Class mappings**: update land-cover codes if sources change.

------------------------------------------------------------------------

## **Usage**

All scripts load configuration as:

``` r
CFG <- cfg_read()
```

The file `config.yml` is written and finalised by `R/00_setup.R` and
should be treated as the authoritative description of a given run tag.
