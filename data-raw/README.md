# Raw Datasets (data-raw)

This folder contains external source datasets used by the workflow.

Files here are treated as read-only inputs. The pipeline reads from data-raw
and writes derived products to data and output.

## Folder layout

```text
data-raw/
├── ESACCI/
├── GLC_FCS30D/
├── LAI/
├── FPAR/
└── LUH2_v2h/
```

## Datasets

### ESA-CCI / C3S land cover (ESACCI)

Used to build fractional land-cover layers and CCI-based masks.

Included content:

- Download scripts: download_landcover_1992_2015.py and
  download_landcover_2016_2022.py.
- Annual maps in ESACCI/ESACCI_1992-2020/.

Source:

- Copernicus Climate Data Store,
  https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=download

### GLC_FCS30D (GLC_FCS30D)

Used for the independent GLC-based masking branch.

Included content:

- Annual GeoTIFF maps.
- Export recipe in R/_ref_google_engine_glc_code.txt.

Sources:

- Zenodo DOI: https://doi.org/10.5281/zenodo.8239305
- GEE catalog page: https://gee-community-catalog.org/projects/glc_fcs/

### LAI and FPAR monthly inputs (LAI, FPAR)

Monthly NetCDF time series used before georeferencing, masking, and
aggregation.

Included content:

- LAI/lai_1982-2024/
- FPAR/fpar_1982-2024/

Source publication:

- Jeong et al. (2024), Remote Sensing of Environment,
  https://doi.org/10.1016/j.rse.2024.114282

### LUH2 v2h (LUH2_v2h)

Land-use state variables used for pasture/rangeland constraints and overlap
logic.

Included content:

- states.nc
- supporting LUH2 documentation files.

Source:

- https://luh.umd.edu/data.shtml

## How these inputs are used

- Build CCI and GLC land-cover fractions and masks.
- Add LUH2-based pasture/rangeland overlap constraints.
- Georeference and prepare monthly LAI/FPAR products.

## Reproducibility notes

- data-raw contents are not tracked by git.
- Reproducibility depends on exact upstream dataset versions.
- Paths, year windows, and class mappings are defined in config/config.yml.
