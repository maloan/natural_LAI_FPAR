# Reference Grids and Auxiliary Files (src)

This folder stores static spatial references and shared auxiliary inputs used
across the whole pipeline.

Nothing here should be run-tag specific. Think of src as common infrastructure
for geometry, area weighting, and provenance.

## Reference grids

Canonical global grids used by processing and analysis:

- ref_0p05.tif / ref_0p05.nc
- ref_0p25.tif / ref_0p25.nc
- ref_0p5.tif / ref_0p5.nc
- ref_0p05_griddes.txt and ref_0p25_griddes.txt for CDO-style remapping.

## Area rasters

Grid-cell area files used for area-weighted aggregation and summaries:

- area_0p05_km2.tif / area_0p05_km2.nc
- area_0p25_km2.tif / area_0p25_km2.nc
- area_0p5_km2.tif / area_0p5_km2.nc

Additional valid-domain area rasters are also present:

- area_0p05_validdomain_km2.nc
- area_0p25_validdomain_km2.nc

## AOI masks

- aoi_0p05.tif
- aoi_0p25.tif

These masks support quicklooks and regional diagnostics.

## Other auxiliary inputs

- ne_110m_coastline.gpkg: coastline geometry for maps.
- cor_twi_vegh_5km_mosaic.nc: ancillary raster used in optional diagnostics.
- valid_tiles_info_0p05_full_10deg.rds: tile metadata for chunked processing.

## Manifest and provenance

- manifest_00.csv is generated during setup and records geometry checks and
	summary statistics.

## Notes

- Files in src are intended to be stable across runs.
- Avoid writing temporary or run-specific outputs here.
- Any changes in src should be treated as pipeline-level changes.
