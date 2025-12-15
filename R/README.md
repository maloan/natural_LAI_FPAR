# **Core R Code**

This directory contains the principal components of the LAI/FPAR processing and analysis workflow. Scripts are organised by processing stage (00–17), supported by shared modules for geometry, I/O, configuration, statistics, and visualisation.

------------------------------------------------------------------------

## **Script overview**

### **00 — Initialisation**

-   **`00_setup.R`** Initialise configuration, verify paths, load shared modules, and prepare runtime options.

------------------------------------------------------------------------

### **01 — Georeferencing (0.05°)**

-   **`01_georef_0p05.R`** Harmonise raw LAI and FPAR inputs to the common 0.05° grid and write monthly georeferenced layers.

------------------------------------------------------------------------

### **02–03 — Fractional cover & GLC stacks**

-   **`02_apply_abiotic_only_0p05.R`** Apply only the abiotic static mask to 0.05° fields (stand-alone utility step).
-   **`03_cci_frac_0p05.R`** Derive ESA-CCI fractional cover layers (grass, urban, agricultural, etc.) at 0.05°.
-   **`05_glc_stack_0p05.R`** Build yearly Landsat GLC_FCS30D stacks (categorical persistence and fractional summaries).

------------------------------------------------------------------------

### **04–10 — Mask construction**

-   **`04_cci_mask_0p05.R`** Dynamic CCI mask based on τ, k thresholds for fractional-cover transitions.
-   **`06_glc_mask_0p05.R`** GLC persistence mask (e.g., used ≥ N years).
-   **`07_abiotic_static_from_cci_0p05.R`** CCI-derived static mask for water, and ice.
-   **`08_abiotic_static_from_glc_0p05.R`** GLC-derived abiotic mask (alternative source).
-   **`09_luh_use_masks.R`** Assemble LUH2-based land-use (pasture) layers.
-   **`10_luh_pasture_overlap_025.R`** Quantify LUH pasture × satellite grass overlap; derive drop masks at 0.25° and propagate to 0.05°.

------------------------------------------------------------------------

### **11–12 — Mask application & aggregation**

-   **`11_apply_mask_0p05.R`** Apply selected mask configuration (CCI/GLC, LUH overlap, abiotic static) to georeferenced LAI/FPAR at 0.05°.
-   **`12_agg_0p25.R`** Area-weighted aggregation of masked 0.05° layers to 0.25° resolution.

------------------------------------------------------------------------

### **13–17 — Analytical and diagnostic scripts**

-   **`13_analyse_georef_0p25.R`** Diagnostic summaries and maps for georeferenced (unmasked) 0.25° LAI/FPAR.
-   **`14_compare_lai_fpar_trends.R`** Compare LAI vs FPAR trend structures and spatial patterns.
-   **`15_analyse_masked_trends.R`** Trend analysis for masked 0.25° products.
-   **`16_compare_masked_unmasked.R`** Quantify differences between masked and unmasked trends.
-   **`17_analyse_trend_mask.R`** Examine trend sensitivity to specific masking schemes.

------------------------------------------------------------------------

## **Shared modules**

### Geometry & grid handling

-   **`geom.R`** Grid alignment, area calculations, coordinate utilities, and spatial helpers.

### Input/output

-   **`io.R`** NetCDF/TIFF read–write helpers, GDAL settings, provenance metadata.

### Visualisation

-   **`viz.R`** Global/AOI quicklooks, mask diagnostics, and standardised plotting utilities.

### Naming, configuration, options

-   **`names.R`** Controlled naming conventions for intermediate products.
-   **`options.R`** Parsing of environment variables and runtime switches.

### General utilities & statistics

-   **`utils.R`** Configuration readers, path utilities, timing wrappers, file management.
-   **`stats.R`** Area-weighted metrics, mask comparison statistics, Jaccard/overlap measures.

------------------------------------------------------------------------

## **Workflow summary**

1.  **Setup** environment configuration.
2.  **Georeference** LAI/FPAR monthly fields (0.05°).
3.  **Derive** fractional cover (CCI) and categorical persistence (GLC).
4.  **Construct** masks (CCI, GLC, LUH pasture overlap, abiotic).
5.  **Apply** masks to monthly LAI/FPAR using flexible combinations.
6.  **Aggregate** masked products to 0.25° for analysis and modelling.
7.  **Analyse** trends, masking effects, and LAI–FPAR consistency (scripts 13–17).

------------------------------------------------------------------------
