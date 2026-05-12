# Analysis Workflows (R/analysis/)

This folder contains R scripts for statistical summaries, visualizations, and
diagnostics of LAI/FPAR trends across different masking scenarios.

## Scientific methodology overview

### Data flow

```
trends/ (shell scripts)
  ↓
  Computes annual metrics, OLS trends, relative trends, MK p-values
  ↓
output/<TAU>/eval/trend_<VAR>_<MASK>/*.nc (trend rasters)
  ↓
R/analysis/*.R (this folder)
  ↓
Global/zonal/class aggregations + bootstrap CI + paper figures
  ↓
analysis/results/ (tables, figures)
```

### Statistical methods

**Global/zonal/class-level aggregations:**
1. Load trend rasters (slope in native units/year, or relative %/year)
2. Restrict to valid domain (nonmissing mask)
3. Weight by pixel area (0.25° = varying km² per latitude)
4. Compute area-weighted mean: Σ(trend × area) / Σ(area)
5. Bootstrap confidence interval: resample pixels with replacement (preserving weights)
6. Significance: mark with * if 95% CI does not cross zero

**Spatial autocorrelation:**
- Preserved through **zonal aggregation** (all pixels in a class treated as a cluster)
- Bootstrap weights respect cluster membership
- CIs account for spatial dependence without explicit adjustment

## Scripts by function

### Global summaries

**02_global_relative_trends_summary.R**
- Global mean relative greening (% per year) across all masking scenarios
- Outputs: CSV with mean ± 95% CI, significance flag
- Variables: LAI, FPAR (both yearmean and yearmax)

**03_global_absolute_trends_summary.R**
- Global mean absolute trends (native units/year)
- LAI: dimensionless units/year; FPAR: fAPAR units/year
- Outputs: CSV with formatted paper table + long-form overview
- Includes: effective area (km²), n_pixels, significance

**04_global_absolute_trends_timeseries.R**
- Time series trajectories and OLS trend fitting
- Outputs: time series CSV + fitted values + OLS statistics
- Figures: time series plot with trend lines and 95% CI bands (PNG/PDF)

### Zonal (latitude band) summaries

**05_zonal_seasonal_amplitude.R**
- Seasonal amplitude (yearmax - yearmin) by latitude band (1° bands)
- Outputs: CSV with zonal mean amplitude across scenarios
- Figures: line plot showing latitudinal variation

**06_zonal_relative_trends_all_masks.R**
- Relative trends by latitude band across all masking scenarios
- Outputs: CSV with zonal mean relative trend ± CI
- Figures: overlay plot showing scenario comparisons

### Stratified (land-cover / climate-zone) summaries

**12_landcover_trend_summary.R**
- Trends by IPCC-aggregated land-cover classes
- Workflow:
  1. Load yearly 0.25° LC majority maps (1992-2020)
  2. Merge subclasses (e.g., cropland mosaics → 10)
  3. Compute "stable" class as mode across years
  4. Bootstrap CI by class (accounts for spatial autocorrelation within class)
- Outputs:
  - CSV: full details (area, trend, CI, retained-area %, significance)
  - CSV: paper-grade summary (selected columns, rounded)
  - Figures: bar plot with error bars + retained-area labels (PNG/PDF)

**13_kg_trend_summary.R**
- Trends by Köppen-Geiger climate zones
- Workflow:
  1. Load trend rasters
  2. Lookup KG zone for each pixel (kgc::LookupCZ, chunked to manage memory)
  3. Aggregate to both 3-letter (full class) and 2-letter (main group)
  4. Bootstrap CI by zone
- Outputs:
  - CSV: full 3-letter classes
  - CSV: aggregated 2-letter groups
  - Figures: bar plots for both aggregation levels

### Diagnostics & support

**00_area_validdomain_after_nonvegetated.R**
- Compute valid-domain area mask (nonmissing, non-water, non-snow)
- Output: `area_0p25_validdomain_km2.nc` used in all trend aggregations

**01_masking_footprint_summary.R**
- Summarize masking impact (how much area removed by each mask type/tau)

**07_nonveg_snapshot_sensitivity.R**
- Sensitivity analysis: impact of different non-vegetated pixel thresholds

**08_lai_yearmean_trend_maps.R**, **09_lai_yearmax_trend_maps.R**
- Spatial trend maps (slope, significance) across all scenarios

**10_cropland_pasture_trend_diagnostics.R**
- Detailed agricultural trend analysis (cropland/pasture subregions)

**11_zonal_diagnostics_combined_figure.R**
- Multi-panel diagnostic figure combining global + zonal information

## Configuration and parameters

### Fixed parameters (in scripts)

| Parameter | Values | Meaning |
|-----------|--------|---------|
| `cci_taus` | τ=0.05, 0.1, 0.2 | Cloud cover threshold options (CCI mask) |
| `glc_run_tag` | τ=0.1 (default) | GLC mask tau setting |
| `band_deg` | 1 | Latitude band width (degrees) |
| `lc_year_start:end` | 1992-2020 | Land-cover years for stability calculation |
| `n_boot` | 1000 | Bootstrap samples (affects CI precision) |
| `conf` | 0.95 | Confidence level (95% CI) |
| `seed` | 42 | Fixed seed for reproducibility |

### Relative trend thresholds (from trends/)

| Variable | EPS | Meaning |
|----------|-----|---------|
| LAI | 0.05 | Min temporal mean for relative trend (avoids 0/0) |
| FPAR | 0.02 | Min temporal mean for relative trend |

### Output directories

```
analysis/results/
  ├── tables/
  │   └── trends/              # Global trend summaries (CSV)
  ├── figures/
  │   ├── timeseries/          # Time-series plots
  │   └── summaries/           # Zonal, LC, KG plots + paper-grade tables
  ├── lc_trends/
  │   └── <TAU>/               # Land-cover summaries per tau
  └── kg_trends/
      └── <TAU>/               # Köppen-Geiger summaries per tau
```

## Running scripts

### Single script
```r
Rscript R/analysis/02_global_relative_trends_summary.R
```

### With environment variable configuration
```bash
# Set relative trend threshold for FPAR
export EPS_FPAR=0.015
Rscript R/analysis/12_landcover_trend_summary.R --tau tau_0.1 --var FPAR

# Configure custom tau
export GLC_RUN_TAG=tau_0.2
Rscript R/analysis/02_global_relative_trends_summary.R
```

### Via Makefile (from repo root)
```bash
make analysis
```

## Output interpretation

### CSV tables

All CSV outputs include:
- **mean_est** (or reltrend_pct_per_year, abs_trend_per_year): point estimate
- **ci_lower, ci_upper**: 95% confidence interval bounds
- **sig_flag**: * = significant (CI does not cross zero), blank = not significant
- **n_pixels**: effective sample size (number of valid pixels)
- **area_km2**: total area in valid domain
- **n_classes** (in stratified outputs): number of subregions

### Figures

Standard styling:
- Color palette: grey20 (unmasked), teal/orange/purple (CCI tau sweep), blue (GLC)
- Error bars: 95% CI from bootstrap resampling
- Significance: implicit (if CI crosses zero, trend is non-significant)

## Reproducibility

- All random seeds are fixed (seed = 42)
- Bootstrap resampling is deterministic
- Paths relative to repo root via `here()` package
- Dependency versions documented in DESCRIPTION

## References for methods

- **Bootstrap**: Efron & Tibshirani (1993), *An Introduction to the Bootstrap*
- **Weighted bootstrap**: Barbe & Bertail (1995), *The Weighted Bootstrap*
- **Mann-Kendall**: Mann (1945), Kendall (1975)
- **Zonal statistics**: terra::zonal() documentation
