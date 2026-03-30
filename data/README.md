# **Intermediate Processed Data (`data/`)**

This directory contains **harmonized intermediate products** generated
by the preprocessing stages of the LAI/FPAR pipeline. These datasets are
analysis-ready, spatially consistent, and treated as **stable inputs**
for all downstream masking, aggregation, and evaluation steps.

All files here are created programmatically from `data-raw/` and are not
edited manually.

------------------------------------------------------------------------

## **Directory structure**

```         

data/
├── frac/
│   ├── cci_frac_0p05/
│   └── glc_frac_0p05/
└── georef/
    ├── georef_lai_0p05/
    └── georef_fpar_0p05/
```

------------------------------------------------------------------------

## **Contents**

### **1) Fractional land-cover products (`frac/`)**

Fractional cover layers at **0.05°**, derived from categorical
land-cover maps.

- **`cci_frac_0p05/`**: Cropland, urban, grass, and fused-class fractions derived from ESA-CCI / C3S land cover.
- **`glc_frac_0p05/`**: Equivalent products derived from GLC_FCS30D.

**Purpose**
- Construction of land-cover masks (threshold- and persistence-based)
- Sensitivity tests comparing CCI- vs GLC-derived land-use information
- Diagnostic plots and consistency checks

------------------------------------------------------------------------

### **2) Georeferenced LAI / FPAR (`georef/`)**

Monthly LAI and FPAR fields at **0.05°**, aligned to the project's
canonical grid.

- **`georef_lai_0p05/`**: LAI monthly rasters after reprojection, snapping, and extent correction.
- **`georef_fpar_0p05/`**: FPAR monthly rasters processed identically.

**Processing steps include**
- Reprojection to EPSG:4326
- Strict alignment to the 0.05° reference grid
- Consistent global extent
- Restoration of a continuous monthly time axis

These products are the **direct inputs** to masking
(`11_apply_mask_0p05.R`) and subsequent aggregation to 0.25°.

------------------------------------------------------------------------

## **Conventions**

- CRS, resolution, and extent conform exactly to the reference grids defined in `config/config.yml`
- Once written, files in `data/` are treated as immutable; regeneration requires rerunning the corresponding preprocessing scripts

------------------------------------------------------------------------

## **Role in the workflow**

`data/` serves as the **boundary between raw inputs and analysis**:
- Upstream scripts populate this directory
- Downstream scripts assume its presence and consistency
- No analysis results are written here (those belong in `output/` and `analysis/`)
