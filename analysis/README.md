# **Analysis**

This directory contains **post-processing and diagnostic outputs**
derived from the masked and unmasked LAI/FPAR products after aggregation
to **0.25°**. The emphasis is on reproducible summaries, comparisons,
and visual diagnostics used in evaluation and manuscript figures.

------------------------------------------------------------------------

## **Scope**

The analysis layer covers:

- Global and regional trend estimation
- Masked vs. unmasked comparisons
- LAI/FPAR consistency diagnostics
- Greening/browning attribution
- Pixel-level trend significance
- Global and zonal time series summaries

All analyses are derived from precomputed 0.25° fields and are fully
scriptable and HPC-compatible.

------------------------------------------------------------------------

## **Key outputs**

### **Baseline data**

- **`unmasked/`**: Unmasked LAI/FPAR fields at 0.05° and 0.25° used as analysis baselines.
- **`tmp/`**: Intermediate processing outputs and temporary data.

### **Current results** (`results/`)

- **`global_mean_relative_trends/`**: Global mean trends and relative changes across mask configurations.
- **`kg_trends/`**: Köppen-Geiger climate zone trend analysis.
- **`lc_trends/`**: Land-cover class trend analysis.
- **`mask_sensitivity_tau/`**: Sensitivity analysis for different masking thresholds.
- **`masks/`**: Processed and combined mask products.
- **`metrics_mask_effect_zonal/`**: Zonal metrics showing mask effects on trends.
- **`paper_figures/`**: Publication-ready figures for manuscripts.
- **`relative_trend_distribution/`**: Distribution of relative trends across space.
- **`spatial_redistribution_after_masking/`**: Analysis of how masking redistributes trends spatially.
- **`zonal_relative_trends/`**: Zonal profiles of relative trends.
- **`dropped_region/`**: Analysis of regions removed by masking.

### **Archived analyses** (`archive/`)

- **`trends_masked/`**: Masked trend products organized by mask configuration and run tag.
- **`trends_dropped/`**: Trends computed only over regions removed by masking.
- **`trend_significance/`**: Pixel-level trend significance testing (masked and unmasked).
- **`pixel_significance/`**: Per-pixel slope, p-values, and FDR-corrected significance layers.
- **`comparison_masked_unmasked/`**: Direct comparisons between masked and unmasked trends.
- **`delta_trends/`**: Δ-trend products (e.g., yearmax − yearmean), with maps and NetCDFs.
- **`global_timeseries_masks/`**: Global mean LAI/FPAR time series across mask configurations.
- **`greening_browning/`**: Greening vs. browning statistics (global and zonal) with summary plots.
- **`lai_vs_fpar_rel/`**: LAI/FPAR cross-variable diagnostics: scatter plots, zonal profiles, and maps.
- **`tables/`**: Tabular outputs for figures or supplementary material.
- **`figures/`**: Diagnostic and comparative figure outputs.
- **`global_relative_trends_paperfig/`**: Relative trend figures prepared for manuscript.
- **`trends_drop_vs_keep/`**: Comparative analysis of trends for dropped vs. retained regions.

------------------------------------------------------------------------

## **Notes**

- All results are derived from products in `output/<run_tag>/masked_0p25/`.
- Naming conventions encode variable, metric (yearmean/yearmax), mask, and run tag.
- Figures and NetCDFs in this directory are intended to be publication-ready.
- This folder contains **results and diagnostics only**; scripts live in `R/` and are documented separately.
- Historical analyses are preserved in `archive/` for reference and reproducibility.
