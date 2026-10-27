# Managed-Land-Masked LAI / FPAR Processing Pipeline

This repository builds global LAI and FPAR products with mapped managed areas excluded from satellite observations. Static masks based on ESA-CCI/C3S and GLC_FCS30D land cover combined with LUH2 pasture data are applied at 0.05° resolution. The retained domain represents areas with lower mapped management influence. Products are aggregated to coarser grids using area weighting for subsequent analysis.

Observed vegetation trends combine responses to environmental change with the effects of land use and management. This pipeline quantifies the sensitivity of trend estimates to excluding mapped managed areas and selecting the retained spatial domain.

## Main outputs

- Monthly masked LAI and FPAR products.
- Binary CCI-based and GLC-based mask layers.
- Area-weighted aggregates and time-series summaries.
- Grid-cell trend and significance products.
- Masked and unmasked comparisons and diagnostics.
- Chapter 2 unmasked FPAR and CCI--pasture mask inputs at 0.5° resolution.

Outputs are organized under `output/<RUN_TAG>/`.

## Repository layout

```text
R/            Processing and analysis scripts
config/       Run configuration
data-raw/     Local external source data
data/         Intermediate harmonized products
output/       Generated products by run tag
analysis/     Analysis outputs and figures
src/          Static reference grids and auxiliary files
```

## Run

Place the raw datasets under `data-raw/`, then run the complete workflow from the repository root:

```bash
make
```

The workflow builds both variables for all four masking scenarios and runs the analyses.

The analysis uses three CCI thresholds under `output/alpha_<threshold>/` and
one GLC-based branch under `output/glc/`. The GLC land-use and pasture masks
are independent of the CCI managed-land alpha threshold. Both scenarios use
the same CCI-derived nonvegetated exclusion.

Build the release archives with:

```bash
bash release/build_archives.sh
```

## Configuration

The four fixed configuration files under `config/` define paths, years, classes, thresholds, and grids.

## Requirements

- R 4.1 or later.
- R packages used by the processing and analysis scripts, including `terra`, `sf`, `ncdf4`, and `tidyverse`.
- R packages `Rcpp`, `sf`, and `lwgeom` for annual 0.25° land-cover-fraction aggregation.
- System libraries including GDAL, PROJ, and NetCDF.

Ubuntu example:

```bash
sudo apt install gdal-bin libgdal-dev libproj-dev libnetcdf-dev
```

## Data sources

Place the ESA-CCI/C3S land cover, GLC_FCS30D, LUH2 v2h, LAI, and FPAR source datasets under `data-raw/`.
