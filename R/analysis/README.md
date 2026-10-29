# Analysis Workflows

This folder contains R scripts for statistical summaries, figures, and diagnostics of LAI trends across masking scenarios. Chapter 1 summaries use LAI. The workflow analyzes three CCI thresholds and one GLC-based masking scenario.

## Data flow

```text
trends/
  ↓
Annual diagnostics, OLS trends, relative trends, MK p-values
  ↓
output/<RUN_TAG>/eval/trend_<VAR>_<MASK>/
  ↓
R/analysis/
  ↓
analysis/results/
```

## Spatial weighting and uncertainty
Spatial summaries weight each valid 0.25° trend estimate by the cell's fixed baseline area after water and permanent snow or ice exclusions. Managed-land masks exclude 0.05° subcells before aggregation and can therefore change the resulting 0.25° values and their availability. The baseline area weights remain fixed. Uncertainty is estimated using a spatial block bootstrap with 5° × 5° blocks. Confidence intervals excluding zero are treated as statistically resolved.

## Key scripts

- `00_area_validdomain_after_nonvegetated.R` creates the fixed post-non-vegetated support rasters.
- `03_matched_period_annual_mean_trends.R` recalculates unmasked annual-mean LAI trends for alternative analysis periods.
- `12_0_landcover_dominant_class.R` assigns each 0.25° grid cell to its dominant grouped 1992--2022 CCI land-cover class.
- Script 14 produces the separate climate-zone and land-cover summary figures used in the manuscript.

## Climate classification

Köppen--Geiger classes are assigned to the 0.25° grid by the dominant area-weighted class among the 9 by 9 nominal 100-arc-second source cells provided by `kgc`. Ocean cells are excluded from the class totals. Detailed three-character classes and two-character climate groups are aggregated and assigned independently. Class assignment uses categorical area dominance.

## Land-cover summaries


`R/12_make_lc025_fractions.R` builds the 31 annual 22-layer fraction rasters from raw 300 m CCI maps.
The workflow combines the 22 primary-class fractions into 19 analysis groups, averages the grouped fractions over 1992--2022, and assigns each cell to the group with the largest mean fraction.
