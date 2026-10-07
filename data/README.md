# Intermediate Processed Data

This folder contains generated intermediate products used by the LAI and FPAR processing pipeline.

Files in this directory are derived from `data-raw/` and should not be edited manually.

## Layout

```text
data/
├── frac/
│   ├── cci_frac_0p05/
│   └── glc_frac_0p05/
└── georef/
    ├── georef_lai_0p05/
    └── georef_fpar_0p05/
```

## Fractional land-cover products

`frac/` contains 0.05° fractional land-cover layers derived from categorical land-cover data.

- `cci_frac_0p05/` contains CCI-based fractions used for mask construction.
- `glc_frac_0p05/` contains GLC-based fractional products used for masking and diagnostics.

## Georeferenced LAI and FPAR

`georef/` contains monthly LAI and FPAR fields aligned to the common 0.05° project grid.

- `georef_lai_0p05/`
- `georef_fpar_0p05/`

These products are the direct inputs to masking and subsequent spatial aggregation.

## Conventions

All files use the project reference grids and are treated as generated outputs. If upstream inputs or preprocessing change, regenerate the affected products using the corresponding pipeline scripts.