# Core R Code

This directory contains the main components of the LAI/FPAR processing workflow.  
Scripts are organised by stage, and shared modules provide geometry, I/O, option
handling, and visualisation tools.

---

## Structure

### 0× — Initialisation and setup
- `00_setup.R`  
  Prepare the runtime environment, load core modules, and validate paths.

### 1× — Geolocation
- `01_georef_0p05.R`  
  Georeference raw LAI/FPAR fields to a 0.05° lon–lat grid.

### 2× — ESA-CCI fractional cover
- `02_cci_frac_0p05.R`  
  Derive annual fractional-cover layers (grass, bare, water, urban, etc.).

### 3× — GLC yearstack
- `03_glc_stack_0p05.R`  
  Build yearly stacks of Landsat GLC_FCS30D (categorical + fractional).

### 4× — Mask construction (CCI, GLC, LUH)
- `04_cci_mask_0p05.R`  
  CCI-based drop-mask using τ,k thresholds.
- `04_glc_mask_0p05.R`  
  GLC used≥N mask (temporal persistence).
- `04_luh_pasture_overlap_025.R`  
  LUH pasture–grass overlap mask (0.25° and 0.05° variants).

### 5× — Static abiotic masks
- `05_abiotic_static_from_cci_0p05.R`  
  Water/ice/bare threshold masks from CCI fractions.
- `05_abiotic_static_from_glc_0p05.R`  
  Same for GLC.

### 7× — Mask application
- `07_apply_mask_0p05.R`  
  Apply binary drop masks to monthly LAI/FPAR (0.05°), with quicklooks.
- `07_luh_use_masks.R`  
  Combine LUH-derived drop masks with other mask sets.

### 8× — Aggregation
- `08_agg_0p25.R`  
  Area-weighted 0.05° → 0.25° aggregation for LAI/FPAR.

---

## Shared modules

### Geometry and grid handling
- `geom.R`  
  Raster alignment, CF-time parsing, longitude rotation, area computation.

### Input/output
- `io.R`  
  GDAL write presets, provenance metadata, session logging, manifest helpers.

### Visualisation
- `viz.R`  
  Standardised global/AOI quicklooks for masks, LAI/FPAR, and fractional cover.

### Naming and configuration
- `names.R`  
  Strict, tokenised naming conventions.
- `options.R`  
  Environment-variable parsing and runtime toggles.

### General utilities
- `utils.R`  
  Path helpers, configuration loading, timing, GDAL options, small geometry tools.

### Mask statistics
- `stats.R`  
  Area-weighted means, drop fractions, and pairwise contingency/Jaccard metrics.

---

## Workflow summary

1. **Georeference raw satellite inputs** (0.05°).  
2. **Construct fractional layers and categorical stacks** (CCI, GLC).  
3. **Generate land-use masks** (CCI τ–k, GLC persistence, LUH pasture overlap).  
4. **Add abiotic static masks** (water/ice/bare).  
5. **Apply masks to LAI/FPAR monthly fields**.  
6. **Aggregate masked data to 0.25°** for model ingestion.  
7. **Evaluate** masks and outputs via area-weighted statistics and quicklooks.

---

