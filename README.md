# Managed-Land-Masked LAI / FPAR Processing Pipeline

This repository builds global LAI and FPAR products with mapped managed areas excluded from satellite observations. Static masks based on ESA-CCI/C3S and GLC_FCS30D land cover combined with LUH2 pasture data are applied at 0.05° resolution. The retained domain is less managed according to these datasets but is not necessarily free of human influence. Products are aggregated to coarser grids using area weighting for subsequent analysis.

Observed vegetation trends combine responses to environmental change with the effects of land use and management. This pipeline quantifies how trend estimates change when mapped managed areas are excluded. The comparison describes sensitivity to the retained spatial domain and does not attribute trend differences to individual drivers.

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
data-raw/     External source data (not tracked)
data/         Intermediate harmonized products
output/       Generated products by run tag
analysis/     Analysis outputs and figures
src/          Static reference grids and auxiliary files
```

## Analysis

Run the numbered processing scripts to build the CCI-based and GLC-based masked trend products first. Then run the downstream analysis with:

```bash
R/analysis/run_all_analysis.sh
```

The runner checks and rebuilds the unmasked baseline when needed and then produces summaries and figures from the completed masked trend products. It does not rebuild the masked processing branches.

The analysis uses three CCI thresholds under `output/alpha_<threshold>/` and
one GLC-based branch under `output/glc/`. The GLC land-use and pasture masks
do not depend on the CCI managed-land alpha threshold. Both scenarios use the
same CCI-derived nonvegetated exclusion.

For GLC processing, select `RUN_TAG=glc` and `MASK=GLC`. For example,
`RUN_TAG=glc MASK=GLC VAR=LAI Rscript R/10_apply_mask_0p05.R` writes the
masked monthly LAI series under `output/glc/`. Use `RUN_TAG=alpha_0.1` with
`MASK=CCI` for the CCI alpha 0.1 counterpart.

Chapter 2 inputs are built separately with:

```bash
Rscript R/12_chapter2_inputs_0p5.R
```

The CCI alpha 0.1 plus pasture masked 0.5° fAPAR series is a separate product:

```bash
RUN_TAG=alpha_0.1 Rscript R/11_agg_0p5.R
```

It writes monthly GeoTIFFs under
`output/alpha_0.1/masked_0p5/FPAR/masked_FPAR_CCI/` and the Chapter 2 input
`output/chapter2/fpar_masked/fpar_masked_0p5_monthly_1982-2024.nc`. The
existing unmasked Chapter 2 fAPAR product remains separate.

## Configuration

Core settings are stored in `config/config_<RUN_TAG>.yml`, including paths, years, land-cover classes, thresholds, and grid references.

## Requirements

- R 4.1 or newer.
- R packages used by the processing and analysis scripts, including `terra`, `sf`, `ncdf4`, and `tidyverse`.
- R packages `Rcpp`, `sf`, and `lwgeom` for annual 0.25° land-cover-fraction aggregation.
- System libraries including GDAL, PROJ, and NetCDF.

Ubuntu example:

```bash
sudo apt install gdal-bin libgdal-dev libproj-dev libnetcdf-dev
```

## Data sources

Place the raw datasets under `data-raw/`, including ESA-CCI/C3S land cover, GLC_FCS30D, LUH2 v2h, LAI, and FPAR. These source datasets are not tracked by Git.

`RAW_INPUT_MANIFEST_2026-10-21.tsv` records the file sizes and SHA-256
checksums of the 516 LAI and 516 fAPAR NetCDFs used for the current source
audit. This distinguishes the new LAI delivery from older files with the same
version label.

GLC_FCS30D v2 is downloaded from Zenodo and processed locally. The workflow derives:

- Categorical modal-class maps used for the GLC-based land-use mask.
- Fractional grass-cover maps used for the LUH2 pasture-overlap calculation.

Download and local-processing instructions are provided in `data-raw/GLC_FCS30D/README.md`.
