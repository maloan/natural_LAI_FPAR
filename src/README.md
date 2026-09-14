# Reference Grids and Auxiliary Files (src)

This folder stores spatial references used across the whole pipeline. 

## Reference grids

Canonical global grids used by processing and analysis:

- ref_0p05.tif / ref_0p05.nc
- ref_0p25.tif / ref_0p25.nc
- ref_0p05_griddes.txt and ref_0p25_griddes.txt for CDO-style remapping.

## Area rasters

Grid-cell area files used for area-weighted aggregation and summaries:

- area_0p05_km2.tif / area_0p05_km2.nc
- area_0p25_km2.tif / area_0p25_km2.nc

Additional valid-domain area rasters are also present:

- area_0p05_validdomain_km2.nc
- area_0p25_validdomain_km2.nc

## Manifest and provenance

- manifest_00.csv is generated during setup and records geometry checks and summary statistics.

