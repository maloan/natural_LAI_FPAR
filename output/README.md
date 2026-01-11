# **Final Outputs (`output/`)**

This directory contains all **generated products** from the
natural-vegetation LAI/FPAR workflow. Outputs are organised by **run
tag** (e.g. `tau_0.05`, `tau_0.1`, `tau_0.2`), where each run tag
corresponds to a specific masking configuration and parameter set.

All downstream analyses and figures draw exclusively from this
directory.

------------------------------------------------------------------------

## **Top-level structure**

```         
output/
├── georef_biotic/
└── tau_<run_tag>/
```

### **`georef_biotic/`**

Georeferenced, *unmasked* LAI and FPAR fields at 0.05°. These serve as
the common baseline for all masking and aggregation steps.

------------------------------------------------------------------------

## **Per-run-tag contents**

Each `tau_<run_tag>/` directory contains the following subdirectories.

### **1. `masked_0p05/`**

Masked monthly LAI and FPAR fields at native **0.05°** resolution.

-   Separate subfolders for `LAI/` and `FPAR/`

-   Masks applied according to the selected configuration (CCI or GLC),
    with optional:

    -   static abiotic exclusions,
    -   LUH pasture–grass overlap removal

-   Used as the direct input for spatial aggregation

------------------------------------------------------------------------

### **2. `masked_0p25/`**

Area-weighted aggregates on the **0.25°** grid.

-   Monthly LAI and FPAR fields aggregated from 0.05°
-   Includes an `eval/` subdirectory with summaries and diagnostics
-   Forms the basis for all trend, zonal, and global analyses

------------------------------------------------------------------------

### **3. `masks/`**

Binary mask layers used in the workflow (`1 = drop`, `0 = keep`):

-   **`mask_cci/`** – CCI used-land masks
-   **`mask_glc/`** – GLC used-land persistence masks
-   **`mask_luh/`** – LUH pasture / rangeland masks
-   **`mask_luh_overlap/`** – pasture–grass overlap masks
-   **`mask_abiotic/`** – static water, ice, and bare-ground masks

Masks are stored primarily at **0.05°** and replicated to **0.25°**
where needed.

------------------------------------------------------------------------

### **4. `eval/`**

Evaluation and trend products derived from masked 0.25° fields.

-   Trend maps and statistics for LAI and FPAR
-   Separate outputs for CCI- and GLC-based masking
-   Inputs to figures, tables, and manuscript analyses
