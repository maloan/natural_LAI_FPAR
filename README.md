# Natural LAI / FPAR Processing Pipeline

This repository builds global natural-vegetation LAI and FPAR products from satellite observations.
We remove areas dominated by human land use first, then analyze trends on the remaining natural-vegetation signal. Masks from ESA-CCI/C3S and GLC_FCS30D are applied at 0.05 degree, and products are aggregated with area weighting to coarser grids for analysis.



Observed vegetation trends mix climate-driven ecosystem responses with land-use effects (cropland expansion, urbanization, management). For attribution and evaluation work, those signals need to be separated. This pipeline focuses on that separation by masking anthropogenic land cover before trend estimation.

## Main outputs

- Monthly masked LAI and FPAR products.
- Binary CCI and GLC mask layers.
- Area-weighted aggregates and time-series summaries.
- Pixel-wise trend and significance products.
- Masked vs unmasked comparisons and diagnostics.
- Chapter 2 unmasked fAPAR and CCI--pasture mask inputs at 0.5 degree.

Outputs are organized under output/<RUN_TAG>/.

## Repository layout

```text
R/            Processing and analysis scripts
config/       Run configuration
data-raw/     External source data (not tracked)
data/         Intermediate harmonized products
output/       Generated products by run tag
analysis/     Analysis outputs and figures
src/          Static reference grids and auxiliary files
vignettes/    Extended documentation
```

## Analysis

Processing scripts are run directly in their numbered order. After the CCI and GLC trend products have been created, run the complete analysis with:

```bash
R/analysis/run_all_analysis.sh
```

The runner executes three CCI thresholds and a single GLC branch stored under `alpha_0.1`. The GLC method itself does not depend on the CCI alpha threshold.

Chapter 2 inputs are built separately with:

```bash
Rscript R/12_chapter2_inputs_0p5.R
```

## Configuration

Core settings are generated in `config/config_<RUN_TAG>.yml` (paths, years, classes, thresholds, and grid references).

## Requirements

- R 4.1 or newer.
- R packages used in scripts (for example terra, sf, ncdf4, tidyverse).
- Python 3 with numpy, netCDF4, rasterio, and pyproj for the single-pass
  land-cover fraction aggregation.
- System libraries: GDAL, PROJ, NetCDF.

Ubuntu example:

```bash
sudo apt install gdal-bin libgdal-dev libproj-dev libnetcdf-dev
```

## Data sources

Place raw datasets under data-raw/ (ESACCI, GLC_FCS30D, LUH2_v2h, LAI, FPAR). These inputs are not tracked by git.

GLC_FCS30D v2 is downloaded from Zenodo and processed locally. The workflow creates:

- categorical modal-class maps for the GLC land-use mask.
- fractional grass-cover maps for the LUH2 pasture consistency check.

The download and validation scripts are stored in `data-raw/GLC_FCS30D/`.

## More documentation

See vignettes/vignette.md and the README files in each subfolder for details.
