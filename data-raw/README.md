# Raw Datasets

This folder contains local copies of the external source datasets used by the processing pipeline.

```text
data-raw/
├── ESACCI/
├── GLC_FCS30D/
├── LAI/
├── FPAR/
└── LUH2_v2h/
```

## ESA-CCI / C3S land cover

Used to derive fractional land-cover layers, CCI-based managed-land masks, and the static non-vegetated mask.

Data are stored under `ESACCI/ESACCI_1992-2022/`. The workflow reads the `lccs_class` variable from the annual NetCDF files.

Versions used:

- 1992--2015: v2.0.7cds
- 2016--2022: v2.1.1

Source:

- https://cds.climate.copernicus.eu/datasets/satellite-land-cover
- DOI: https://doi.org/10.24381/cds.006f2c9a

Download scripts are provided in the `ESACCI/` directory.

## GLC_FCS30D v2

Used for the GLC-based masking branch and for fractional grass-cover estimates used with LUH2 pasture data.

The original 30 m archives are stored under `GLC_FCS30D/archives/`. Derived 0.05° products include categorical modal-class maps and fractional grass-cover maps.

Source:

- Liu, L., Zhang, X., and Li, Z. (2025). GLC_FCS30D (version 2). Zenodo. https://doi.org/10.5281/zenodo.15063683
- Zhang et al. (2024), *Earth System Science Data*, 16, 1353--1381. https://doi.org/10.5194/essd-16-1353-2024

## LAI and FPAR

Monthly LAI and FPAR inputs cover 1982--2024 and are stored under:

```text
LAI/lai_1982-2024/
FPAR/fpar_1982-2024/
```

Source:

- https://www.environment.snu.ac.kr/data/longterm-lai
- Jeong et al. (2024), *Remote Sensing of Environment*, 311, 114282. https://doi.org/10.1016/j.rse.2024.114282
- Jeong et al. (2026), *Sustained global greening driven by continuous CO2 fertilization*, Research Square preprint, version 1. https://doi.org/10.21203/rs.3.rs-8621169/v1

## LUH2 v2h

LUH2 land-use state variables are used to derive the pasture-overlap mask.

Main input:

```text
LUH2_v2h/states.nc
```

Source:

- https://luh.umd.edu/data.shtml
- Hurtt et al. (2020), Harmonization of global land use change and management for the period 850–2100 (LUH2) for CMIP6. *Geoscientific Model Development*, 13, 5425–5464. https://doi.org/10.5194/gmd-13-5425-2020

