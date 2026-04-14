# Trend Workflows (trends)

This folder contains standalone shell and R scripts for trend products,
significance testing, and relative-trend calculations.

You can run these scripts independently from the main Makefile pipeline when
you want focused trend processing.

## Scripts in this folder

### build_georef_products.sh

Builds trend-ready products from unmasked georeferenced LAI/FPAR inputs.

Example:

```bash
./build_georef_products.sh LAI 0p05
./build_georef_products.sh FPAR 0p05
```

### build_trends_masked_0p25.sh

Computes masked trend products at 0.25 degree for a specific run tag and mask
source.

Example:

```bash
./build_trends_masked_0p25.sh tau_0.2 FPAR GLC
./build_trends_masked_0p25.sh tau_0.1 LAI CCI
```

### generated_masked_trends.sh

Batch helper for generating masked trend products across multiple run/mask
combinations.

### make_relative_trends_georef.sh

Creates relative-trend products from unmasked georeferenced 0.25 degree data.

### compute_mk_pval.R

Computes pixel-wise Mann-Kendall p-values from annual stacks.
Supports masked and unmasked modes.

### mk_sig_mask.R

Command-line utility that converts a time-series raster to a significance mask.

Example:

```bash
Rscript mk_sig_mask.R in_ts.nc out_mask.nc alpha [mask.nc]
```

## Common dependencies

- cdo
- gdal_translate
- Rscript with terra and trend packages

## Typical use in this project

1. Build baseline unmasked trend products.
2. Build masked 0.25 degree trends by run tag and mask source.
3. Generate relative trends and MK significance products for diagnostics.

## Notes

- SNU_LAI_FPAR_ROOT can override the default repository root in scripts that
	support it.
- Relative-trend thresholds (for example EPS_LAI and EPS_FPAR) are configurable
	via environment variables in shell workflows.
