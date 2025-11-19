# **Reference Grids and Auxiliary Source Files**

This directory contains core spatial reference products, AOI definitions, grid descriptors, and auxiliary datasets required by the LAI/FPAR processing workflow. All files here are static and shared across scripts in the `R/` and `workflow/` directories.

## **Reference grids**

These define the canonical spatial geometry for all rasters:

-   **`ref_0p05.tif` / `ref_0p05.nc`** Global 0.05° lon–lat grid (EPSG:4326), used as the native resolution for masks and monthly fields.

-   **`ref_0p25.tif` / `ref_0p25.nc`** Global 0.25° grid, used for area-weighted aggregation and trend analysis.

-   **`ref_0p05_griddes.txt`, `ref_0p25_griddes.txt`** CDO/ESMF-compatible grid descriptor files for remapping operations.

## **Cell-area rasters**

Required for area-weighted computations:

-   **`area_0p05_km2.{tif,nc}`**
-   **`area_0p25_km2.{tif,nc}`**

Each contains per-cell surface area in km².

## **Areas of interest (AOIs)**

-   **`aoi_0p05.tif`**, **`aoi_0p25.tif`** Binary masks defining global and regional AOIs. Used by quicklook scripts and the `analysis/` pipeline.

## **Additional auxiliary data**

-   **`ne_110m_coastline.gpkg`** Optional coastline overlay for quicklook maps.

-   **`cor_twi_vegh_5km_mosaic.nc`** External ancillary dataset (e.g., topographic or hydrologic covariate), used in selected diagnostics.

-   **`valid_tiles_info_0p05_full_10deg.rds`** Internal metadata on tile validity and 10° tiling structure, used for chunked processing.

-   **`manifest_00.csv`** Initial-stage manifest summarizing grid properties and provenance.

## **Notes**

-   All reference files are assumed to be stable across runs.
-   No dynamic products are written to this directory; outputs instead go to `data/`, `output/`, or `analysis/`.

------------------------------------------------------------------------
