# **Raw Datasets (`data-raw/`)**

This directory contains **external source datasets** used as inputs to
the natural-vegetation LAI/FPAR workflow. Files here are treated as
**immutable**: the pipeline reads from them but does not modify them in
place. All derived products are written to `data/` and `output/`.

## **Directory layout**

```         

data-raw/
├── ESACCI/
├── GLC_FCS30D/
├── LAI/
├── FPAR/
└── LUH2_v2h/
```

------------------------------------------------------------------------

## **Datasets**

### **1) ESA-CCI / C3S Land Cover (`ESACCI/`)**

Annual land-cover products used to derive cropland/urban fractions and
the CCI used-land mask.

**Contents** - Download scripts: `download_landcover_1992_2015.py`,
`download_landcover_2016_2022.py` - Annual maps (GeoTIFF) under
`ESACCI/ESACCI_1992-2020/`

**Source** - ESA Climate Change Initiative. (n.d.). *Land cover (CCI /
C3S)* [Data set]. Copernicus Climate Data Store.
<https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=download>

*Note:* This includes ESA-CCI v2.0.7 (1992–2015) and C3S v2.1.1
(2016–2022). The analysis window used by the pipeline is configured in
`config/config.yml`.

------------------------------------------------------------------------

### **2) GLC_FCS30D (`GLC_FCS30D/`)**

Annual land-cover maps used for GLC-based persistence masking (and
optionally fractional cover).

**Contents** - GeoTIFF exports at 0.05° (e.g.,
`GLC_FCS30D_0p05deg_1985.tif` … `GLC_FCS30D_0p05deg_2022.tif`) - Export
recipe used in this project: `R/000_google_engine_glc_code.txt`

**Source** - Liu, L., Zhang, X., & Zhao, T. (2023). *GLC_FCS30D: The
first global 30-m land-cover dynamic monitoring product with fine
classification system from 1985 to 2022* [Data set]. Zenodo.
<https://doi.org/10.5281/zenodo.8239305>\
- GEE Community Catalog. (n.d.). *GLC_FCS30D project page*.
<https://gee-community-catalog.org/projects/glc_fcs/>

------------------------------------------------------------------------

### **3) LAI and fAPAR monthly inputs (`LAI/`, `FPAR/`)**

Monthly NetCDF time series used as the observational basis before
georeferencing, masking, and aggregation.

**Contents** - `LAI/lai_1982-2024/` (monthly NetCDF) -
`FPAR/fpar_1982-2024/` (monthly NetCDF)

**Source** - Jeong, S., Ryu, Y., Gentine, P., Lian, X., Fang, J., Li,
X., Dechant, B., Kong, J., Choi, W., Jiang, C., Keenan, T. F., Harrison,
S. P., & Prentice, I. C. (2024). Persistent global greening over the
last four decades using novel long-term vegetation index data with
enhanced temporal consistency. *Remote Sensing of Environment, 311*,
114282. <https://doi.org/10.1016/j.rse.2024.114282>

------------------------------------------------------------------------

### **4) LUH2 v2h (`LUH2_v2h/`)**

Land-use state variables and transitions used for pasture/rangeland
fractions and overlap masking.

**Contents** - `states.nc`, `transitions.nc`, `management.nc`,
`staticData_quarterdeg.nc` - `LUH2_v2h_README.pdf`

**Source** - Land-Use Harmonization Project. (n.d.). *LUH2 data and
documentation* [Data set]. <https://luh.umd.edu/data.shtml>

------------------------------------------------------------------------

## **How these inputs are used (high level)**

These datasets provide inputs for:

-   Land-cover fractions and used-land masks (CCI, GLC),

-   LUH2 pasture/rangeland constraints and overlap diagnostics,

-   LAI/FPAR georeferencing and time-series construction.

Derived products are written to `data/` (intermediate) and `output/`
(final + evaluation).

------------------------------------------------------------------------

## **Versioning and reproducibility**

-   `data-raw/` contents are not tracked by Git.
-   Reproducibility depends on upstream versions and the local file
    inventory.
-   The pipeline’s paths, years, and class mappings are defined in
    `config/config.yml`. \`\`\`
