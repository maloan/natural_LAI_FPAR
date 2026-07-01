# Natural LAI / FPAR Pipeline Vignette

This vignette gives a short tour of the workflow used in this repository. The aim is to build LAI and FPAR products that reflect natural vegetation by masking out land-use dominated areas before aggregation and analysis. Satellite vegetation data mixes climate responses with land-use effects. If those signals are analyzed together, it is hard to tell what is driving a trend. This pipeline separates them by applying land-cover masks first.

The workflow does four things:

1.  prepares shared reference grids,
2.  georeferences monthly LAI and FPAR,
3.  builds and applies land-use masks,
4.  aggregates and analyzes the masked products.

Two independent masking branches are used:

- ESA-CCI/C3S land cover
- GLC_FCS30D land cover

Both are applied at 0.05 degree, then aggregated to 0.25 degree for analysis.

## Main steps

### 1. Set up the project

Script: [R/00_setup.R](../R/00_setup.R)

This creates the reference grids, area rasters, and project configuration. Run it once for each run tag.

### 2. Georeference LAI and FPAR

Script: [R/01_georef_0p05.R](../R/01_georef_0p05.R)

This aligns the monthly LAI and FPAR NetCDF files to the 0.05 degree grid.

### 3. Prepare land-cover inputs

Scripts:

- [R/02_cci_frac_0p05.R](../R/02_cci_frac_0p05.R)
- [R/04_glc_stack_0p05.R](../R/04_glc_stack_0p05.R)

These create the land-cover inputs used by the masking scripts.

### 4. Build masks

Scripts:

- [R/03_cci_mask_0p05.R](../R/03_cci_mask_0p05.R)
- [R/05_glc_mask_0p05.R](../R/05_glc_mask_0p05.R)
- [R/06_nonveg_static_from_cci_0p05.R](../R/06_nonveg_static_from_cci_0p05.R)
- [R/09_luh_pasture_overlap_0p25.R](../R/09_luh_pasture_overlap_0p25.R)

Mask convention:

- 1 = drop
- 0 = keep
- NA = undefined

### 5. Apply masks and aggregate

Scripts:

- [R/10_apply_mask_0p05.R](../R/10_apply_mask_0p05.R)
- [R/11_agg_0p25.R](../R/11_agg_0p25.R)
- [R/11_agg_0p5.R](../R/11_agg_0p5.R)

These produce the analysis-ready 0.25 degree products.

### 6. Run trend analysis

Trend scripts live in [trends/](../trends/).

Typical outputs include annual summaries, trend slopes, significance tests, masked vs unmasked comparisons, and zonal or regional summaries.

## Where outputs go

- [data/](../data/) holds intermediate products.
- [output/](../output/) holds masked data, evaluation products, and quicklooks.
- [analysis/](../analysis/) holds summaries and figures.

Outputs are grouped by run tag, such as `tau_0.1` or `tau_0.2`.

## Running it

From [R/](../R/):

``` bash
make pipeline
```

To include the analysis steps:

``` bash
make analysis
```

To use only one masking branch:

``` bash
make pipeline MASKS=CCI
make pipeline MASKS=GLC
make pipeline RUN_TAG=tau_0.1
```

To regenerate quicklooks:

``` bash
make ql
```

## Customization

Most project settings live in [config/config.yml](../config/config.yml). Common changes are run tags, year ranges, AOIs, mask thresholds, and class mappings.

If you change those settings, rerun setup and any downstream steps that depend on them.

## Further reading

- [README.md](../README.md)
- [R/README.md](../R/README.md)
- [data/README.md](../data/README.md)
- [data-raw/README.md](../data-raw/README.md)
- [src/README.md](../src/README.md)
- [analysis/README.md](../analysis/README.md)
- [trends/README.md](../trends/README.md)

Last updated: 2026-07-01
