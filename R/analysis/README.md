# Analysis Workflows (R/analysis/)

This folder contains R scripts for statistical summaries, visualizations, and diagnostics of LAI/FPAR trends across different masking scenarios.

## Run the complete analysis

After all masked trend products are complete and validated, run:

```bash
R/analysis/run_all_analysis.sh
```

The runner uses the paper's fixed analysis design: three CCI thresholds and one GLC branch stored under `alpha_0.1`. If the water/ice-only baseline is absent, it first applies the static non-vegetated mask, aggregates monthly fields using the same area-weighted calculation as the masked branches, and builds the unmasked trend products. It then generates the 1995, 2007, and 2022 sensitivity masks and runs scripts 00--13 in dependency order Land-cover fractions are aggregated by `R/12_make_lc025_fractions.R`. Land-cover and Köppen-Geiger summaries are run for both absolute and relative trends.

## Data flow

```
trends/ (shell scripts)
  ↓
  Computes annual metrics, OLS trends, relative trends, MK p-values
  ↓
output/<ALPHA>/eval/trend_<VAR>_<MASK>/*.nc (trend rasters)
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
5. Bootstrap confidence interval: resample 5° × 5° spatial blocks with replacement
6. Significance: mark with * if 95% CI does not cross zero
