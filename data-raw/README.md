# **Raw Datasets and Source Inputs**

This directory contains unprocessed external datasets used as primary inputs to the LAI/FPAR workflow. Files here retain their original structure (or minimal preprocessing only) and are not modified by the main pipeline.

## **Subdirectories**

-   **`ESACCI/`** ESA CCI / C3S Land Cover products and helper scripts for downloading and yearstack preparation.

-   **`FPAR/`** Original FPAR time series from SNU before reprojection or clipping.

-   **`LAI/`** Original LAI fields from SNU. Used as the starting point for georeferencing and temporal stacking.

-   **`GLC_FCS30D/`** GLC_FCS30D land cover at 0.05.

-   **`LUH2_v2h/`** Raw LUH2 v2h NetCDF files, including state variables, transitions, and management layers.

## **Usage**

These datasets act as immutable inputs for:

-   Fractional cover derivation
-   Land-cover masks (CCI, GLC)
-   LUH2 pasture/rangeland fraction and overlap analysis
-   LAI/FPAR georeferencing and time-series construction

All subsequent steps read from the processed `data/` directory.

## **Versioning**

These files are **not tracked by Git**. Reproducibility depends on external dataset versioning and documented download procedures.

------------------------------------------------------------------------
