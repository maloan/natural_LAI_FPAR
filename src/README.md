# Reference Grids and Auxiliary Files

This folder contains spatial reference files used throughout the processing and analysis pipeline.

## Reference grids

Canonical global grids:

- `ref_0p05.tif` / `ref_0p05.nc`
- `ref_0p25.tif` / `ref_0p25.nc`
- `ref_0p05_griddes.txt` / `ref_0p25_griddes.txt` for CDO grid definitions

## Area rasters

Grid-cell area rasters used for area-weighted aggregation and summaries:

- `area_0p05_km2.tif` / `area_0p05_km2.nc`
- `area_0p25_km2.tif` / `area_0p25_km2.nc`

Valid-domain area rasters:

- `area_0p05_validdomain_km2.nc`
- `area_0p25_validdomain_km2.nc`

## Manifest

`manifest_00.csv` is generated during setup and records grid geometry checks and summary statistics.