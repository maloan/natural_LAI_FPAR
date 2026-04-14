# R Helpers

This folder contains shared helper modules used by scripts in R/.
They provide reusable functions for paths, file discovery, raster I/O,
NetCDF handling, options parsing, and plotting support.

## Modules

- paths.R: path expansion and config loading.
- files.R: filename parsing and file discovery.
- netcdf.R: NetCDF variable handling and raster alignment helpers.
- io.R: GDAL creation options and write helpers.
- netcdf_raster.R: NetCDF-to-raster ingestion and time-axis handling.
- options.R: environment-variable parsing for runtime flags.
- plotting.R: shared plotting styles and quicklook helpers.
- utils.R: compatibility wrapper that sources core helper modules.
- viz.R: compatibility wrapper that sources plotting.R.

## Usage Notes

- These scripts are intended to be sourced by pipeline scripts.
- They are not meant to be executed directly as standalone jobs.
- Keep helper functions generic and side-effect free when possible.

mk_sig_mask.R is intentionally outside helpers because it is a standalone
command-line script, not a sourced helper module.
