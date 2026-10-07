# Analysis Workflows

This folder contains R scripts for statistical summaries, figures, and diagnostics of LAI and FPAR trends across masking scenarios.

## Run the analysis

After the masked trend products have been generated, run:

```bash
R/analysis/run_all_analysis.sh
```

The runner uses three CCI thresholds and one GLC-based branch stored under `alpha_0.1`. It also checks and rebuilds the unmasked baseline when required.

## Data flow

```text
trends/
  ↓
Annual diagnostics, OLS trends, relative trends, MK p-values
  ↓
output/<ALPHA>/eval/trend_<VAR>_<MASK>/
  ↓
R/analysis/
  ↓
analysis/results/
```

## Spatial weighting and uncertainty

Spatial summaries use the fixed post-non-vegetated support area of each valid 0.25° grid cell.

Land-use masks determine whether a grid-cell estimate is retained, but do not reduce its baseline support area.

Area-weighted means are calculated as:

```text
Σ(trend × area) / Σ(area)
```

Uncertainty is estimated using a spatial block bootstrap with 5° × 5° blocks. Confidence intervals excluding zero are treated as statistically resolved.

## Key scripts

- `00_area_validdomain_after_nonvegetated.R` creates the fixed post-non-vegetated support rasters.
- `03_matched_period_annual_mean_trends.R` recalculates unmasked annual-mean LAI trends for alternative analysis periods.
- `12_0_landcover_dominant_class.R` assigns each 0.25° grid cell to its dominant 1992--2022 CCI land-cover class.
- Script 14 produces climate-zone and land-cover summary figures.

## Climate classification

Köppen--Geiger classes are assigned categorically to 0.25° grid-cell centres using the nominal 100-arc-second classification provided by `kgc`. Climate classes are not spatially interpolated.

## Land-cover summaries

Each 0.25° grid cell is assigned to the CCI class with the largest mean fractional cover over 1992--2022. Land-cover fractions are used only for class assignment. Spatial summaries use the same fixed post-non-vegetated area weights as the global analysis.