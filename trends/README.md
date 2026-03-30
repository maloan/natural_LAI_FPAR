# **Trend Analysis Workflows (`trends/`)**

This directory contains **standalone shell and R scripts** for computing temporal trends and derived products from LAI/FPAR datasets. These scripts support both **georeferenced (unmasked) baseline products** and **masked trend analysis**, providing flexible workflows for trend estimation, significance testing, and relative-trend normalization.

All scripts are designed to be **independent and reusable**, accepting command-line arguments or environment variables for configuration.

------------------------------------------------------------------------

## **Scripts**

### **`build_georef_products.sh`**

Builds georeferenced LAI/FPAR products at multiple resolutions from input biotic fields.

**Usage**
```bash
./build_georef_products.sh LAI 0p05
./build_georef_products.sh FPAR 0p05
```

**Arguments**
- `VAR`: Variable name (`LAI` or `FPAR`)
- `RES`: Input resolution (e.g., `0p05`)

**Processing**
- Reads georeferenced biotic fields from `output/georef_biotic/`
- Produces trend products at native 0.05° and aggregated 0.25° resolutions
- Computes Mann–Kendall significance masks
- Outputs to `analysis/unmasked/{0p05,0p25}/`

**Dependencies**: `gdal_translate`, `cdo`, `Rscript`

------------------------------------------------------------------------

### **`build_trends_masked_0p25.sh`**

Computes masked (land-use restricted) trend products at 0.25° resolution for a specific mask configuration and run tag.

**Usage**
```bash
./build_trends_masked_0p25.sh tau_0.2 FPAR GLC
./build_trends_masked_0p25.sh tau_0.1 LAI  CCI
```

**Arguments**
- `TAU`: Run tag / threshold folder (e.g., `tau_0.2`)
- `VAR`: Variable name (`LAI` or `FPAR`)
- `MASKTAG`: Mask source (`CCI` or `GLC`)

**Processing**
- Reads masked 0.25° fields from `output/{TAU}/masked_0p25/`
- Computes pixel-level trends and significance
- Derives absolute and relative trend products
- Applies thresholds for relative trends (default: LAI = 0.05, FPAR = 0.02)
- Outputs to `output/{TAU}/eval/trend_{VAR}_{MASKTAG}/`

**Environment variables**
- `EPS_LAI`: Threshold for LAI relative trends (default: 0.05)
- `EPS_FPAR`: Threshold for FPAR relative trends (default: 0.02)

**Dependencies**: `gdal_translate`, `cdo`

------------------------------------------------------------------------

### **`make_relative_trends_georef.sh`**

Creates relative (baseline-normalized) trend rasters from unmasked georeferenced 0.25° fields.

**Usage**
```bash
VARS="LAI FPAR" METRICS="yearmean yearmax yearamp" bash make_relative_trends_georef.sh
```

**Environment variables**
- `VARS`: Space-separated list of variables (default: `"LAI FPAR"`)
- `METRICS`: Space-separated list of metrics (default: `"yearmean yearmax yearmin yearamp"`)
- `EPS_LAI`: Threshold for LAI (default: 0.05)
- `EPS_FPAR`: Threshold for FPAR (default: 0.02)

**Processing**
- Normalizes trends by baseline values to compute relative changes
- Generates trend maps at multiple aggregation levels
- Outputs NetCDF files to `analysis/unmasked/0p25/`

**Dependencies**: `cdo`

------------------------------------------------------------------------

### **`compute_mk_pval.R`**

Computes pixel-level Mann–Kendall p-values for time series stacks, enabling statistical significance assessment of temporal trends.

**Configuration** (edit script header)
```R
mode <- "masked"        # or "unmasked"
var <- "LAI"            # or "FPAR"
metric <- "yearmean"    # time series metric (yearmean, yearmax, etc.)
tau <- "tau_0.1"        # run tag (for masked mode)
mask <- "CCI"           # mask source (CCI or GLC)
```

**Processing**
- Reads annual-aggregated time series from NetCDF
- Computes pixel-wise Mann–Kendall trend test
- Outputs p-value raster to NetCDF
- Supports both masked and unmasked modes

**Input paths**
- **Unmasked**: `analysis/unmasked/0p25/{VAR}_georef_{METRIC}_0p25.nc`
- **Masked**: `output/{TAU}/eval/trend_{VAR}_{MASK}/{VAR}_{METRIC}_0p25.nc`

**Output paths**
- **Unmasked**: `analysis/unmasked/0p25/{VAR}_georef_{METRIC}_trend_mk_pval_0p25.nc`
- **Masked**: `output/{TAU}/eval/trend_{VAR}_{MASK}/{VAR}_{METRIC}_trend_mk_pval_0p25.nc`

**Dependencies**: `terra`, `trend` (R packages); `Rscript`

------------------------------------------------------------------------

## **Workflow integration**

These scripts are typically invoked:

1. **During pipeline setup**: `build_georef_products.sh` generates unmasked baseline products for all variables and resolutions.
2. **After masked output**: `build_trends_masked_0p25.sh` computes trends for each run tag and mask configuration.
3. **For diagnostics**: `make_relative_trends_georef.sh` and `compute_mk_pval.R` generate derived products for analysis and visualization.

All outputs integrate with the analysis framework in `analysis/` and `output/`.

------------------------------------------------------------------------

## **Notes**

- Environment variable `SNU_LAI_FPAR_ROOT` can override the default repo root path (defaults to `$HOME/GitHub/natural_LAI_FPAR`).
- Trend thresholds (EPS_LAI, EPS_FPAR) should be tuned according to data characteristics and analysis goals.
