# **Configuration**

This directory contains the central configuration used by all components of the LAI/FPAR workflow. The file `config.yml` defines paths, grid metadata, temporal coverage, environmental thresholds, and AOI specifications. All R scripts load it at run time via `cfg_read()`.

## **Contents**

-   **`config.yml`** Primary project configuration:

    -   input data locations (CCI, GLC, LUH2, etc.)
    -   output directories (masks, quicklooks, aggregates, analysis)
    -   canonical grid definitions for 0.05° and 0.25° (reference rasters and area rasters)
    -   variable limits for LAI/FPAR clamping
    -   AOI definitions for global and regional quicklooks
    -   project time span (e.g., CCI 1992–2020; GLC yearstack range)

-   **`README.md`** Overview of how configuration variables interact with the pipeline.

## **Customization**

Typical adjustments include:

-   **Paths** Point to local or HPC data directories for raw CCI, GLC, LUH, and intermediate outputs.

-   **AOIs** Add or modify named AOIs to control which regions appear in mask quicklooks and evaluation plots.

-   **Thresholds and parameters**

    -   CCI mask thresholds (`tau_cci`, `k_cci`, band selection)
    -   GLC used-years threshold (`used_n_years`)
    -   Abiotic limits (`tau_water`, `tau_ice`, `tau_bare`)
    -   LUH grass–pasture rules (`g_min`, `p_min`, `alpha`)
    -   Temporal averaging window for LUH (`luh_avg_start`, `luh_avg_end`)

-   **Grid metadata** Update reference and area rasters if using alternative resolutions or custom grid definitions.

## **Usage**

All R scripts expect the environment variable:

``` bash
export SNU_LAI_FPAR_ROOT="path/to/natural_LAI_FPAR"
```

and load configuration via:

``` r
CFG <- cfg_read()
```

This ensures consistent behaviour across preprocessing, masking, aggregation, and analysis steps.

------------------------------------------------------------------------
