# Analysis Workflows (R/analysis/)

This folder contains R scripts for statistical summaries, visualizations, and diagnostics of LAI/FPAR trends across different masking scenarios.

## Run the complete analysis

After all masked trend products are complete and validated, run:

```bash
R/analysis/run_all_analysis.sh
```

The runner uses the paper's fixed analysis design: three CCI thresholds and one GLC branch stored under `alpha_0.1`. If the water/ice-only baseline is absent, it first applies the static non-vegetated mask, aggregates monthly fields using the same area-weighted calculation as the masked branches, and builds the unmasked trend products. It then runs the required analysis scripts in dependency order. Exact unmasked and retained land-cover class areas are derived directly on the 0.05° mask grid by `12_0_landcover_class_area_weights.R`.Script 14 produces separate two-panel climate-zone and land-cover figures ordered by masked trend.

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
3. Weight the unmasked reference by post-nonvegetated area and each masked scenario by its 0.05° retained area aggregated to 0.25°
4. Compute area-weighted mean: Σ(trend × area) / Σ(area)
5. Bootstrap confidence interval: resample 5° × 5° spatial blocks with replacement
6. Significance: mark with * if 95% CI does not cross zero

`00_area_validdomain_after_nonvegetated.R` creates the base valid-domain area, the scenario-specific retained-area rasters, and their QA tables. 

`12_0_landcover_class_area_weights.R` first calculates the 1992--2022 mean class area on the 0.05° mask grid and then aggregates the retained class area to 0.25°. This preserves the association between each class and the mask inside partially retained cells. Retained-versus-excluded trend contrasts refer only to fully excluded 0.25° cells because trends of excluded fine-cell fractions inside partially retained cells are not available.
