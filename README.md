# Natural LAI/FPAR Processing Pipeline

This project implements an end-to-end, reproducible workflow for processing, masking, and aggregating remote sensing products of LAI and FPAR at 0.05° and 0.25° resolution. It leverages ESA-CCI land cover, GLC-FCS30D, and LUH2 data sets to isolate "natural" vegetation signals over time.

## Structure Overview

-   `R/` – Core scripts, functions, and processing logic
-   `data-raw/` – Raw data archives and download utilities
-   `data/` – Intermediate processed inputs (e.g., frac, georef)
-   `output/` – Final outputs (masked rasters, aggregated stacks)
-   `config/` – Centralized YAML config
-   `src/` – Grids, masks, auxiliary geodata
-   `analysis/` – Figures, Makefile-driven reproducibility
-   `vignettes/` – Long-form documentation and example walkthroughs

## Setup

-   R ≥ 4.2, packages: `terra`, `yaml`, `stringr`, `ggplot2`, `dplyr`, `glue`
-   Optional Python utils in `data-raw/ESACCI/`

## Reproducibility

Use the Makefile in `analysis/` to run the full pipeline:

``` bash
make all
```
