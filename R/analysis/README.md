# Analysis Workflows (R/analysis/)

This folder contains R scripts for statistical summaries, visualizations, and diagnostics of LAI/FPAR trends across different masking scenarios.

## Data flow

```
trends/ (shell scripts)
  ↓
  Computes annual metrics, OLS trends, relative trends, MK p-values
  ↓
output/<TAU>/eval/trend_<VAR>_<MASK>/*.nc (trend rasters)
  ↓
R/analysis/*.R (this folder)
  ↓
analysis/results/ (tables, figures)
```

## Statistical methods

**Global/zonal/class-level aggregations:**
1. Load trend rasters (slope in native units/year, or relative %/year)
2. Restrict to valid domain (nonmissing mask)
3. Weight by pixel area (0.25° = varying km² per latitude)
4. Compute area-weighted mean: Σ(trend × area) / Σ(area)
5. Bootstrap confidence interval: resample pixels with replacement (preserving weights)
6. Significance: mark with * if 95% CI does not cross zero
