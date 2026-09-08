# Trend Workflows (trends)

Standalone scripts for computing OLS trends, Mann-Kendall significance, and relative trends on LAI/FPAR data at 0.25° resolution.

## Quick Start

**Unmasked trends** (one-time setup):
```bash
./01_build_georef_products.sh
```

**Masked trends** (for specific mask combination):
```bash
./build_trends_masked_0p25.sh alpha_0.1 LAI CCI
```

**Batch masked trends** (all combinations):
```bash
./02_batch_build_trends_masked.sh
```


## What These Scripts Do

### 1. Convert GeoTIFFs to NetCDF
Each script reads monthly GeoTIFFs from the masked input directories and converts them to a single monthly time series NetCDF file (1982–2024, 516 timesteps).

### 2. Compute Annual Metrics
Four annual metrics are calculated from the monthly data:
- **yearmean**: annual average
- **yearmax**: annual maximum
- **yearmin**: annual minimum
- **yearamp**: range (max - min)

### 3. OLS Trends
Linear regression (using CDO's `trend` function) calculates:
- Slope (change per year)
- Intercept

### 4. Relative Trends (%/year)
Relative trend = (slope / temporal_mean) × 100, but only computed where the temporal mean exceeds a threshold:
- LAI: EPS = 0.05
- FPAR: EPS = 0.02
- yearamp: EPS = 0.01

This avoids division by near-zero values in non-vegetated pixels.

### 5. Mann-Kendall Significance
Pixel-wise Mann-Kendall test (using R package `trend`):
- Tests whether each pixel shows a monotonic trend
- Outputs p-value (no multiple comparison correction)
- Significance at α = 0.05

## Scripts in this folder

### 01_build_georef_products.sh

Builds unmasked georeferenced trend products: from raw GeoTIFFs through annual metrics, OLS trends, Mann-Kendall p-values, and relative trends.

**Usage:**
```bash
./01_build_georef_products.sh
```

**Workflow:**
1. Convert monthly GeoTIFFs (0.05°) → time-stamped NetCDF
2. Compute annual metrics: yearmean, yearmax, yearmin, yearamp
3. OLS trends: slope (per-year) and intercept for each metric
4. Remap to 0.25°
5. Mann-Kendall p-values (pixel-wise, unmasked)
6. Relative trends: (slope / temporal_mean) × 100, only where mean ≥ EPS

**Outputs** (in `analysis/unmasked/0p25/`):
- Annual metrics: `*_georef_yearmean_0p25.nc`, etc.
- Trends: `*_georef_yearmean_trend_slope_peryear_0p25.nc`, etc.
- Relative trends: `*_georef_yearmean_trend_relative_percent_peryear_0p25.nc`, etc.
- Significance: `*_georef_yearmean_trend_mk_pval_0p25.nc`, etc.

### build_trends_masked_0p25.sh

Computes masked trend products at 0.25° for a specific run tag and mask source. Primary workflow for analysis.

**Usage:**
```bash
./build_trends_masked_0p25.sh alpha_0.2 FPAR GLC
./build_trends_masked_0p25.sh alpha_0.1 LAI CCI
```

**Arguments:**
- ALPHA: mask folder name (e.g., alpha_0.1, alpha_0.2)
- VAR: LAI or FPAR
- MASKTAG: CCI or GLC (mask source)

**Workflow:**
1. GeoTIFF → monthly NetCDF: convert monthly GeoTIFFs to time-stamped NetCDF
2. Annual metrics: compute yearmean, yearmax, yearmin, yearamp using CDO
3. OLS trends: linear regression slope (per-year) and intercept for each metric
4. Relative trends: (slope / temporal_mean) × 100 (only where mean ≥ EPS)
5. MK significance: pixel-wise Mann-Kendall p-values (parallel)

**Outputs** (in `output/<ALPHA>/eval/trend_<VAR>_<MASKTAG>/`):
- Monthly: `<VAR>_masked_monthly_0p25.nc`
- Annual metrics: `<VAR>_yearmean_0p25.nc` etc.
- Absolute trends: `<VAR>_yearmean_trend_slope_peryear_0p25.nc` etc.
- Relative trends: `<VAR>_yearmean_trend_relative_percent_peryear_0p25.nc` etc.
- MK p-values: `<VAR>_yearmean_mk_pval_0p25.nc` etc.

### 02_batch_build_trends_masked.sh

Batch wrapper for generating masked trend products across multiple (ALPHA, VAR, MASKTAG) combinations.

**Usage:**
```bash
./02_batch_build_trends_masked.sh
```

The script directly runs CCI thresholds 0.05, 0.1, and 0.2 for LAI and FPAR
with both CCI and GLC masks.

### compute_mk_pval.R

Computes pixel-wise Mann-Kendall p-values from annual stacks (called automatically by the main scripts).

This internal script receives explicit command-line arguments from the two
workflow scripts above.

## Common dependencies

- `cdo` (Climate Data Operators)
- `gdal_translate` (GDAL)
- `Rscript` with packages: terra, here, trend
