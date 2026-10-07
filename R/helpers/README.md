# R Helpers

Shared helper modules used by scripts in `R/`.

These modules provide reusable functions for raster and NetCDF handling, area-weighted summaries, bootstrap confidence intervals, climate classification, command-line arguments, and plotting.

## Modules

- `netcdf.R` for NetCDF handling and raster alignment.
- `io.R` for raster I/O, path construction, numeric utilities, and file writing.
- `plotting.R` for shared plotting functions and styles.
- `cli_args.R` for command-line argument parsing and scenario helpers.
- `weighted_means.R` for area-weighted aggregation.
- `bootstrap_ci.R` for bootstrap confidence intervals.
- `kg_classification.R` for Köppen--Geiger classification, grid lookup, and classification diagnostics.

## Usage

These files are intended to be sourced by other pipeline scripts and are not run directly.