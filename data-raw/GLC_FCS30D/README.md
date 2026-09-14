# GLC_FCS30D v2

The workflow uses the original 30 m GLC_FCS30D version 2 archives from Zenodo record [15063683](https://doi.org/10.5281/zenodo.15063683).

## Download

The record contains 36 ZIP archives (approximately 135.8 GB). To be placed into `data-raw/GLC_FCS30D/archives/`. This directory is ignored by Git; the archives remain local raw data and are not committed.

## Local aggregation

Each nominal 5-degree tile is stored as a three-band file for 1985, 1990, and 1995 and a 23-band annual file for 2000--2022. To run the global processing use

```bash
Rscript R/04_glc_native_to_0p05.R
Rscript R/04_glc_stack_0p05.R
```
