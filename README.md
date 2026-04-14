# Natural LAI / FPAR Processing Pipeline

This repository builds global natural-vegetation LAI and FPAR products from
satellite observations.

The main idea is simple: remove areas dominated by human land use first, then
analyze trends on the remaining natural-vegetation signal.

Masks from ESA-CCI/C3S and GLC_FCS30D are applied at 0.05 degree, and products
are aggregated with area weighting to coarser grids for analysis.

## Why this workflow exists

Observed vegetation trends mix climate-driven ecosystem responses with land-use
effects (cropland expansion, urbanization, management). For attribution and
evaluation work, those signals need to be separated.

This pipeline focuses on that separation by masking anthropogenic land cover
before trend estimation.

## Main outputs

- Monthly masked LAI and FPAR products.
- Binary CCI and GLC mask layers.
- Area-weighted aggregates and time-series summaries.
- Pixel-wise trend and significance products.
- Masked vs unmasked comparisons and diagnostics.

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

## Pipeline at a glance

1. Setup reference grids and metadata (00_setup.R).
2. Georeference monthly LAI/FPAR (01_georef_0p05.R).
3. Prepare land-cover inputs (02_cci_frac_0p05.R, 04_glc_stack_0p05.R).
4. Build masks (03_cci_mask_0p05.R, 05_glc_mask_0p05.R,
   06_abiotic_static_from_cci_0p05.R).
5. Apply masks and aggregate (10_apply_mask_0p05.R,
   11_agg_0p25.R and 11_agg_0p5.R).
6. Run trend and diagnostic analysis scripts in R/analysis/.

## Quick start

Set project root (optional):

```r
Sys.setenv(SNU_LAI_FPAR_ROOT = "~/path/to/natural_LAI_FPAR")
```

Run from R/:

```bash
make pipeline
make analysis
```

Useful variants:

```bash
make pipeline MASKS=CCI
make pipeline RUN_TAG=tau_0.2
make ql
```

## Configuration

Core settings are generated in config/config.yml (paths, years, classes,
thresholds, grid references).

## Requirements

- R 4.1 or newer.
- R packages used in scripts (for example terra, sf, ncdf4, tidyverse).
- System libraries: GDAL, PROJ, NetCDF.

Ubuntu example:

```bash
sudo apt install gdal-bin libgdal-dev libproj-dev libnetcdf-dev
```

## Data sources

Place raw datasets under data-raw/ (ESACCI, GLC_FCS30D, LUH2_v2h, LAI, FPAR).
These inputs are not tracked by git.

## More documentation

See vignettes/vignette.md and the README files in each subfolder for details.
