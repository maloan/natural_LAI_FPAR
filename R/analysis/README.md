# Analysis Workflows (R/analysis/)

This folder contains R scripts for statistical summaries, visualizations, and diagnostics of LAI/FPAR trends across different masking scenarios.

## Run the complete analysis

After all masked trend products are complete and validated, run:

```bash
R/analysis/run_all_analysis.sh
```

The runner uses the paper's fixed analysis design: three CCI thresholds and one GLC branch stored under `alpha_0.1`. If the water/ice-only baseline is absent, it first applies the static non-vegetated mask, aggregates monthly fields using the same area-weighted calculation as the masked branches, and builds the unmasked trend products. It then runs the required analysis scripts in dependency order. Land-cover summaries assign each 0.25° cell to its dominant 1992--2022 CCI class and weight every valid cell by the fixed area remaining after the common non-vegetated mask. Script 14 produces separate two-panel climate-zone and land-cover figures ordered by masked trend.

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
3. Weight every valid 0.25° estimate by its fixed post-nonvegetated support area, independently of the land-use scenario and the number of subcells retained by the land-use mask
4. Compute area-weighted mean: Σ(trend × area) / Σ(area)
5. Bootstrap confidence interval: resample 5° × 5° spatial blocks with replacement
6. Significance: mark with * if 95% CI does not cross zero

`00_area_validdomain_after_nonvegetated.R` creates the 0.05° and 0.25° post-nonvegetated domain rasters used to define the analysis domain. The 0.25° raster is the fixed statistical weight for every spatial summary: land-use masks determine whether a coarse-cell estimate is present, but do not reduce its baseline support area.

`03_matched_period_annual_mean_trends.R` recalculates unmasked annual-mean LAI trends for the literature-comparison windows. It reads `analysis/unmasked/0p25/LAI_georef_yearmean_0p25.nc` and `src/area_0p25_validdomain_km2.nc`, estimates each grid-cell OLS slope with CDO, and reports the post-nonvegetated-area-weighted global mean and 95% spatial block-bootstrap CI. Its output is `analysis/results/tables/trends/matched_period_unmasked_annual_mean_LAI_trends.csv`. The identical CSV used to typeset Appendix Table 16 is stored in the Paper project under `tables/trends/`.

`12_0_landcover_dominant_class.R` assigns each 0.25° cell to the CCI class with the largest mean fractional cover over 1992--2022. Class summaries and retained-versus-excluded contrasts use the fixed post-nonvegetated area of valid 0.25° cells. Land-cover fractions are used only for dominant-class assignment.

## Climate classification

The analysis assigns one Köppen--Geiger class to each 0.25° LAI cell centre
using the nominal 100-arc-second classification bundled with `kgc`. The lookup
is categorical and does not interpolate climate classes. The LAI values,
trends, fixed post-nonvegetated weights, and bootstrap blocks remain on their established
analysis grids.
