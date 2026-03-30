# Natural LAI / FPAR Processing Pipeline — Comprehensive Vignette

## Overview

This vignette walks through the **Natural LAI / FPAR Processing Pipeline**, a reproducible workflow for deriving global natural-vegetation leaf area index (LAI) and fraction of absorbed photosynthetically active radiation (FPAR) by explicitly masking out anthropogenic land-use signals.

The pipeline addresses a fundamental challenge in satellite vegetation monitoring: satellite products integrate both **climate-driven ecosystem responses** and **land-use change impacts** (e.g., cropland expansion, deforestation, urbanization). This workflow disentangles these signals by applying land-cover masks before aggregation and trend analysis.

------------------------------------------------------------------------

## Scientific Background

### The Problem: Mixed Signals in Satellite Vegetation Data

Satellite-observed vegetation changes can reflect:

1. **Natural ecosystem responses** to climate variability, rising CO₂, and nitrogen deposition
2. **Anthropogenic effects** from land-use change, agricultural intensification, and urbanization

For **attribution studies**, **climate model evaluation**, and **nature-based solution assessments**, these signals must be separated. Analyzing natural vegetation in isolation reveals the true climate signal without confounding land-use artifacts.

### The Solution: Land-Use Masking

This pipeline constructs independent, persistence-based land-use masks from two satellite land-cover products:

- **ESA-CCI Land Cover**: Time-varying categorical maps (1992–2020)
- **GLC-FCS30D**: High-resolution annual maps (1985–2022)

Masks encode cropland, urban, and pasture/rangeland extents at native satellite resolution, preserving spatial detail. Both are applied to monthly LAI/FPAR fields, yielding independent estimates of natural-vegetation trends for robust attribution.

------------------------------------------------------------------------

## Workflow Architecture

### Stage 1: Setup and Reference Grids

**Script**: `R/00_setup.R`

Initializes the analysis framework:

- Defines canonical global reference grids at 0.05° (native) and 0.25° (analysis) resolutions
- Generates cell-area rasters for explicit area weighting
- Defines areas of interest (AOIs) for regional diagnostics
- Writes configuration to `config/config.yml`

This stage must run **once per run tag** (masking configuration).

**Key outputs**:
- `src/ref_0p05.tif`, `src/ref_0p05.nc` — Reference grid at native resolution
- `src/ref_0p25.tif`, `src/ref_0p25.nc` — Reference grid at analysis resolution
- `src/area_0p05_km2.nc`, `src/area_0p25_km2.nc` — Cell-area rasters
- `config/config.yml` — Project metadata and AOI definitions

### Stage 2: Georeferencing Raw LAI / FPAR

**Script**: `R/01_georef_0p05.R`

Converts monthly LAI and FPAR NetCDF time series (from satellite observations) to aligned 0.05° global rasters:

- Reprojects to EPSG:4326 (WGS84)
- Snaps to the canonical 0.05° grid
- Ensures consistent global extent and temporal axis
- Fills gaps and flags invalid data

**Key outputs**:
- `data/georef/georef_lai_0p05/` — Monthly LAI at 0.05°
- `data/georef/georef_fpar_0p05/` — Monthly FPAR at 0.05°

These products serve as the **direct inputs to masking** and are treated as immutable intermediate data.

### Stage 3: Land-Cover Preprocessing

**Scripts**: `R/03_cci_frac_0p05.R`, `R/05_glc_stack_0p05.R`

#### ESA-CCI Fractional Products

`R/03_cci_frac_0p05.R` derives fractional cover layers from ESA-CCI land-cover categorical maps:

- Cropland fraction (% of grid cell)
- Urban fraction
- Grass fraction
- Bare/sparse fraction

Fractions are aggregated annually over the CCI time window (1992–2020).

**Key outputs**:
- `data/frac/cci_frac_0p05/` — Fractional cover at 0.05°

#### GLC-FCS30D Stacking

`R/05_glc_stack_0p05.R` regrids annual GLC-FCS30D categorical maps (30 m → 0.05°):

- Derives equivalent fractional products (cropland, urban, pasture)
- Supports time-varying mask construction
- Maintains consistency with CCI approach

**Key outputs**:
- `data/frac/glc_frac_0p05/` — Equivalent fractional products from GLC

### Stage 4: Mask Construction

**Scripts**: 
- `R/04_cci_mask_0p05.R` — ESA-CCI mask
- `R/06_glc_mask_0p05.R` — GLC-FCS30D mask
- `R/07_abiotic_static_from_cci_0p05.R` — Static abiotic mask
- `R/10_luh_pasture_overlap_025.R` — LUH2 pasture diagnostics

#### Natural-Vegetation Mask Construction

Two independent masks are built from CCI and GLC land-cover data:

1. **Definition**: A grid cell is flagged for masking (removal) if it contains:
   - Cropland (>threshold fraction)
   - Urban area (>threshold fraction)
   - Pasture/rangeland (from LUH2 pasture extent)
   - Abiotic exclusions (water, snow, ice)

2. **Persistence-based approach**: Multi-year averaging to distinguish persistent land use from interannual crop rotation variability

3. **Binary format**: 
   - **1 = drop** (anthropogenic)
   - **0 = keep** (natural)
   - **NA = undefined** (valid domain edge or incomplete data)

#### Static Abiotic Mask

`R/07_abiotic_static_from_cci_0p05.R` creates a permanent mask for non-vegetation areas:

- Water bodies (inland, coastal)
- Snow and ice (polar regions, high mountains)
- Bare rock and deserts (optional, climate-dependent)

This mask is applied uniformly across all years.

**Key outputs**:
- `output/{RUN_TAG}/masks/cci_mask_0p05.nc` — ESA-CCI derived mask
- `output/{RUN_TAG}/masks/glc_mask_0p05.nc` — GLC-FCS30D derived mask
- `output/{RUN_TAG}/masks/abiotic_static_0p05.nc` — Abiotic exclusions

### Stage 5: Mask Application and Temporal Filtering

**Script**: `R/11_apply_mask_0p05.R`

Applies selected masks to monthly LAI and FPAR:

- Loads monthly georeferenced fields
- Sets masked pixels to NA (not analyzed)
- Preserves temporal continuity for time-series analysis
- Writes masked monthly fields

Masking is **non-destructive**: original fields remain available for alternative masking strategies.

**Key outputs**:
- `output/{RUN_TAG}/masked_0p05/{LAI,FPAR}/masked_{VAR}_*.nc` — Masked monthly fields

### Stage 6: Area-Weighted Aggregation

**Script**: `R/12_agg_0p25.R`

Aggregates masked 0.05° fields to 0.25° using explicit area weighting:

- Computes weighted means accounting for partial cells at domain boundaries
- Preserves spatial heterogeneity while reducing computational load
- Produces analysis-ready time series

**Key outputs**:
- `output/{RUN_TAG}/masked_0p25/{LAI,FPAR}/` — Aggregated monthly fields at 0.25°
- `data/georef/georef_{LAI,FPAR}_0p25/` — Unmasked aggregates for reference

### Stage 7: Trend Analysis and Diagnostics

**Scripts**: `R/13_analyse_georef_0p25.R`, `R/14_compare_lai_fpar_trends.R`, etc.

Computes temporal trends and generates diagnostic outputs:

#### Annual Aggregation

Raw monthly fields are aggregated to annual metrics:

- **yearmean**: Annual mean LAI/FPAR
- **yearmax**: Annual maximum (peak growing season)
- **yearmin**: Annual minimum (dormant season)
- **yearamp**: Annual amplitude (difference between max and min)

#### Trend Estimation

Pixel-wise linear regression over the full time series:

- Slope: change per year (units/year)
- p-value: significance via Mann–Kendall test
- FDR: false-discovery rate correction for multiple testing

#### Comparative Analysis

- **Masked vs. Unmasked**: Direct pixel-by-pixel comparison of trends
- **CCI vs. GLC**: Sensitivity to mask source
- **Dropped regions**: Trends computed only on masked-out pixels to quantify anthropogenic effects

#### Zonal and Global Summaries

- Global mean time series (with uncertainty bounds)
- Zonal means by latitude band
- Regional aggregates (optional AOIs: India, China, Europe, etc.)
- Köppen–Geiger climate zone analysis

**Key outputs**:
- `analysis/results/global_mean_relative_trends/` — Global time series and trends
- `analysis/results/zonal_relative_trends/` — Zonal profiles
- `analysis/results/kg_trends/` — Climate-zone stratified analysis
- `output/{RUN_TAG}/eval/` — Pixel-level trend estimates and significance

------------------------------------------------------------------------

## Running the Pipeline

### Prerequisites

1. **Software dependencies**:
   - R (≥4.1) with packages: `terra`, `sf`, `tidyverse`, `ncdf4`, `parallel`
   - Command-line tools: `gdal_translate`, `cdo`
   - Bash shell

2. **External datasets** (not version controlled):
   - LAI/FPAR monthly NetCDFs (in `data-raw/LAI/`, `data-raw/FPAR/`)
   - ESA-CCI Land Cover (in `data-raw/ESACCI/`)
   - GLC-FCS30D annual maps (in `data-raw/GLC_FCS30D/`)
   - LUH2 v2h land-use harmonization data (in `data-raw/LUH2_v2h/`)

3. **Configuration**:
   - Edit `config/config.yml` to set project metadata, AOIs, and run tag

### Full Pipeline Execution

From the repository root:

```bash
cd R
make analysis
```

This sequentially executes:

1. Setup (once per run tag)
2. Georeferencing
3. Land-cover preprocessing
4. Mask construction
5. Mask application
6. Aggregation
7. Analysis

### Partial Execution

To run only data preparation:

```bash
make pipeline
```

To run with specific masks:

```bash
make pipeline MASKS=CCI   # ESA-CCI only
make pipeline MASKS=GLC   # GLC-FCS30D only
```

To regenerate quicklooks:

```bash
make ql
```

### Manual Script Execution

Individual scripts can be run independently:

```bash
# Setup
Rscript R/00_setup.R

# Georeferencing
Rscript R/01_georef_0p05.R

# Land-cover preprocessing
Rscript R/03_cci_frac_0p05.R
Rscript R/05_glc_stack_0p05.R

# Mask construction
Rscript R/04_cci_mask_0p05.R
Rscript R/06_glc_mask_0p05.R
Rscript R/07_abiotic_static_from_cci_0p05.R
Rscript R/09_luh_use_masks.R
Rscript R/10_luh_pasture_overlap_025.R

# Masking and aggregation
Rscript R/11_apply_mask_0p05.R
Rscript R/12_agg_0p25.R

# Analysis (subset of available scripts)
Rscript R/13_analyse_georef_0p25.R
Rscript R/14_compare_lai_fpar_trends.R
Rscript R/15_analyse_masked_trends.R
Rscript R/16_compare_masked_unmasked.R
Rscript R/17_analyse_trend_mask.R
```

------------------------------------------------------------------------

## Understanding the Outputs

### Directory Structure

```
output/{RUN_TAG}/
├── masked_0p05/
│   ├── LAI/
│   │   └── masked_LAI_*.nc          (monthly, masked)
│   └── FPAR/
│       └── masked_FPAR_*.nc         (monthly, masked)
├── masked_0p25/
│   ├── LAI/
│   │   └── masked_LAI_*.nc          (monthly, aggregated)
│   └── FPAR/
│       └── masked_FPAR_*.nc
├── masks/
│   ├── cci_mask_0p05.nc             (ESA-CCI mask)
│   ├── glc_mask_0p05.nc             (GLC-FCS30D mask)
│   └── abiotic_static_0p05.nc
└── eval/
    ├── trend_LAI_CCI/
    │   ├── LAI_yearmean_0p25.nc
    │   ├── LAI_yearmax_0p25.nc
    │   └── *_trend_mk_pval_0p25.nc  (significance)
    └── trend_FPAR_GLC/
        └── ...
```

### Key Variables in Trend Output

**For each variable, metric, and mask combination:**

- `{VAR}_{METRIC}_0p25.nc`: Time series at 0.25° resolution
- `{VAR}_{METRIC}_trend_slope_0p25.nc`: Linear trend (units/year)
- `{VAR}_{METRIC}_trend_pval_0p25.nc`: Mann–Kendall p-value
- `{VAR}_{METRIC}_trend_fdr_0p25.nc`: FDR-adjusted p-value

### Interpretation

1. **Positive slope + low p-value**: Significant greening
2. **Negative slope + low p-value**: Significant browning
3. **High p-value**: No significant trend
4. **Compare CCI vs. GLC**: Assess mask sensitivity
5. **Compare masked vs. unmasked**: Quantify land-use effect

------------------------------------------------------------------------

## Example Use Cases

### Use Case 1: Estimating Natural Greening Trends

**Goal**: Determine the contribution of climate to observed vegetation increases, excluding land-use effects.

**Steps**:

1. Load masked LAI yearmean trend from `output/tau_0.1/eval/trend_LAI_CCI/`
2. Mask pixels with p-value > 0.05 (non-significant)
3. Calculate global mean slope (units/year)
4. Compare with unmasked trend to quantify the effect of land-use masking

**Interpretation**: If masked trend is smaller than unmasked, land use has amplified the observed greening signal.

### Use Case 2: Regional Attribution in India

**Goal**: Assess natural vegetation dynamics in India while accounting for agricultural land.

**Steps**:

1. Load time series from `analysis/results/` filtered to India AOI
2. Compare yearmean and yearmax trends separately (seasonal dynamics)
3. Examine CCI vs. GLC mask sensitivity (crop rotation vs. permanent cropland)
4. Overlay with climate or CO₂ data to attribute drivers

### Use Case 3: Evaluating Land-Use Impacts

**Goal**: Quantify the contribution of anthropogenic land use to observed trends.

**Steps**:

1. Load trends for both masked (natural) and unmasked (all land) datasets
2. Compute difference: trend_unmasked − trend_masked
3. Spatial analysis: where is land use amplifying/dampening greening?
4. Temporal analysis: has the anthropogenic effect changed over time?

------------------------------------------------------------------------

## Customization and Extensions

### Changing Mask Thresholds

Edit the fraction thresholds in:

- `R/04_cci_mask_0p05.R`: CCI mask construction
- `R/06_glc_mask_0p05.R`: GLC mask construction

Example: to exclude pixels with **>5% urban cover** instead of >10%, modify:

```R
frac_urban_threshold <- 0.05  # was 0.10
```

Re-run mask construction and all downstream stages:

```bash
Rscript R/04_cci_mask_0p05.R
Rscript R/11_apply_mask_0p05.R
Rscript R/12_agg_0p25.R
# ... analysis scripts
```

### Adding Custom AOIs

Edit `config/config.yml`:

```yaml
aois:
  my_region:
    lon_min: 10.0
    lon_max: 15.0
    lat_min: 45.0
    lat_max: 50.0
```

Re-run setup:

```bash
Rscript R/00_setup.R
```

Custom AOIs will be available in analysis scripts and quicklooks.

### Extending Time Series

If new LAI/FPAR data becomes available:

1. Add monthly NetCDFs to `data-raw/{LAI,FPAR}/`
2. Update year range in `config/config.yml`
3. Re-run georeferencing and all downstream stages

------------------------------------------------------------------------

## Quality Control and Diagnostics

The pipeline generates multiple **quicklooks** for quality assurance:

### Monthly Products

- Georeferenced LAI/FPAR with global extent
- Masked vs. unmasked comparison at sample months
- Mask footprints (CCI vs. GLC)

### Annual Aggregates

- Yearmean, yearmax, yearmin heatmaps
- Mask consistency across years
- Abiotic exclusion coverage

### Trends

- Pixel-level trend maps (slope, p-value)
- Zonal mean time series with uncertainty bands
- Spatial distribution of significant trends

### Regional

- AOI-specific time series and trends
- Land-cover composition by region
- Anthropogenic effect magnitude

All quicklooks are in `output/{RUN_TAG}/quicklooks/` and `analysis/results/paper_figures/`.

------------------------------------------------------------------------

## References and Links

- **Repository**: https://github.com/maloan/natural_LAI_FPAR
- **ESA-CCI Land Cover**: https://www.esa-landcover-cci.org/
- **GLC-FCS30D**: https://gee-community-catalog.org/projects/glc_fcs/
- **LAI/FPAR observations**: [Jeong et al., 2024 - See data-raw/README.md]

------------------------------------------------------------------------

**Getting Help**

- **Questions on setup**: See [README.md](../README.md)
- **Script documentation**: See [R/README.md](../R/README.md)
- **Data structure**: See [data/README.md](../data/README.md)
- **Input data**: See [data-raw/README.md](../data-raw/README.md)
- **Reference grids**: See [src/README.md](../src/README.md)
- **Analysis outputs**: See [analysis/README.md](../analysis/README.md)
- **Trend workflows**: See [trends/README.md](../trends/README.md)
- **Reporting bugs**: Open an issue on GitHub

------------------------------------------------------------------------

**Last updated**: 2026-03-30  
