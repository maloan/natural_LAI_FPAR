# Analysis Workflows

This folder contains R scripts for statistical summaries, figures, and diagnostics of LAI trends across masking scenarios. The current Chapter 1 analysis runner does not summarize FPAR trends.

## Run the analysis

After the masked trend products have been generated, run:

```bash
R/analysis/run_all_analysis.sh
```

The runner uses three CCI thresholds under `output/alpha_<threshold>/` and one
GLC-based branch under `output/glc/`. It checks the LAI water/ice-only baseline
against its inputs and code and rebuilds that baseline when required.

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
- `12_0_landcover_dominant_class.R` assigns each 0.25° grid cell to its dominant grouped 1992--2022 CCI land-cover class.
- Script 14 produces the separate climate-zone and land-cover summary figures used in the manuscript.

## Climate classification

Köppen--Geiger classes are assigned to the 0.25° grid by the dominant area-weighted class among the 9 by 9 nominal 100-arc-second source cells provided by `kgc`. Ocean cells are excluded from the class totals. Detailed three-character classes and two-character climate groups are aggregated and assigned independently. Climate classes are not spatially interpolated.

## Land-cover summaries


`R/12_make_lc025_fractions.R` builds the 31 annual 22-layer fraction rasters from raw 300 m CCI maps.
