# Intermediate Processed Data (data)

This folder stores the cleaned, harmonized intermediate products used by the LAI/FPAR workflow.
Everything here is generated from data-raw and should not be edited by hand.

## Directory Layout

```text
data/
├── frac/
│   ├── cci_frac_0p05/
│   └── glc_frac_0p05/
└── georef/
    ├── georef_lai_0p05/
    └── georef_fpar_0p05/
```

## What Each Folder Contains

### Fractional Land-Cover Products (frac)

Fractional cover layers at 0.05 degrees derived from categorical land-cover maps.

- cci_frac_0p05/: Cropland, urban, grass, and fused fractions from ESA-CCI/C3S.
- glc_frac_0p05/: Equivalent fractional products from GLC_FCS30D.

These layers are used to build land-cover masks, compare CCI- and GLC-based masking behavior, and support diagnostics.

### Georeferenced LAI and FPAR (georef)

Monthly LAI and FPAR fields at 0.05 degrees, aligned to the project reference grid.

- georef_lai_0p05/: LAI rasters after reprojection, alignment, and extent correction.
- georef_fpar_0p05/: FPAR rasters processed the same way.

Processing includes:

- Reprojection to EPSG:4326.
- Strict alignment to the 0.05 degree reference grid.
- Consistent global extent.
- Continuous monthly time indexing.

These are direct inputs to masking and then to aggregation at coarser resolution.

## Conventions

- CRS, resolution, and extent must match the reference grids defined in config/config.yml.
- Files in data are treated as immutable outputs. If something changes upstream, regenerate them with the relevant scripts.

## Role in the workflow

The data folder is the boundary between preprocessing and analysis:

- Upstream preprocessing scripts write files here.
- Downstream masking and trend scripts read from here.
- Final analysis outputs belong in output and analysis, not in data.
