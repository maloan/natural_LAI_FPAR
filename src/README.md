# **Reference Grids and Auxiliary Source Files (`src/`)**

This directory contains **static spatial reference products and
auxiliary datasets** used throughout the natural-vegetation LAI/FPAR
processing pipeline. All files here define *geometry, area weighting,
AOIs, and provenance* and are shared across scripts in `R/`,
`analysis/`, and the Makefile-driven workflow.

No dynamic or run-specific outputs are written to this directory.

------------------------------------------------------------------------

## **Reference grids**

Canonical global grids used by all raster products:

- **`ref_0p05.tif` / `ref_0p05.nc`**: Global 0.05° longitude–latitude grid (EPSG:4326). Native resolution for land-use masks, georeferenced monthly LAI/FPAR, and all masking operations.
- **`ref_0p25.tif` / `ref_0p25.nc`**: Global 0.25° longitude–latitude grid (EPSG:4326). Target resolution for area-weighted aggregation, trend analysis, and zonal/global diagnostics.
- **`ref_0p05_griddes.txt`**, **`ref_0p25_griddes.txt`**: Grid descriptor files compatible with CDO / ESMF remapping tools. Used for external or command-line regridding.

------------------------------------------------------------------------

## **Cell-area rasters**

Per-cell surface area rasters used for explicit area weighting:

- **`area_0p05_km2.tif` / `area_0p05_km2.nc`**
- **`area_0p25_km2.tif` / `area_0p25_km2.nc`**

Values represent grid-cell area in km² and are derived consistently from the reference grids. These rasters are required for aggregation from 0.05° → 0.25°, global and regional means, and zonal statistics.

------------------------------------------------------------------------

## **Areas of Interest (AOIs)**

Binary AOI masks aligned to the canonical grids:

- **`aoi_0p05.tif`**
- **`aoi_0p25.tif`**

Used for quicklook generation, regional diagnostics, and optional masking during analysis. AOI definitions are configured in `config/config.yml` and rasterized during setup.

------------------------------------------------------------------------

## **Auxiliary datasets**

Additional static inputs used by selected workflow stages:

- **`ne_110m_coastline.gpkg`**: Natural Earth coastline dataset for map overlays and quicklooks.
- **`cor_twi_vegh_5km_mosaic.nc`**: External ancillary dataset (e.g., topographic or hydrologic covariate) used in optional diagnostics or exploratory analyses.
- **`valid_tiles_info_0p05_full_10deg.rds`**: Internal metadata describing valid 0.05° tiles and the 10° tiling scheme, used to support chunked or tiled processing.

------------------------------------------------------------------------

## **Manifests and provenance**

- **`manifest_00.csv`**: Automatically generated during initial setup. Summarizes grid dimensions/extents, resolution consistency, global and AOI surface areas, and creation timestamp. This file provides a provenance checkpoint for the spatial reference framework.

------------------------------------------------------------------------

## **Notes**

- All files in this directory are **static and version-stable**.
- They are shared across all runs and `RUN_TAG`s.
- No outputs or intermediate products should be written here.
- Modifications should be rare and treated as **pipeline-level changes**.
