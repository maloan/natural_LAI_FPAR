# R Helpers

This folder contains shared helper modules used by scripts in R/.
They provide reusable functions for paths, file discovery, raster I/O,
NetCDF handling, options parsing, and plotting support.

## Modules

- netcdf.R: NetCDF variable handling and raster alignment helpers.
- io.R: shared raster I/O, path builders, safe numeric utilities, and write helpers.
- options.R: environment-variable parsing for runtime flags.
- plotting.R: shared plotting styles and quicklook helpers.
- cli_args.R: command-line argument parsing plus standardized scenario table helpers.
- weighted_means.R: area-weighted aggregation utilities.
- bootstrap_ci.R: bootstrap confidence interval utilities.


## Usage Notes

- These scripts are intended to be sourced by pipeline scripts.
- They are not meant to be executed directly as standalone jobs.
