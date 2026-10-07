# Raw Datasets

This folder contains external source datasets used by the processing pipeline. Raw data are not tracked by Git.

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

The original 30 m archives are stored locally under `GLC_FCS30D/archives/` and are not tracked by Git. Derived 0.05° products include categorical modal-class maps and fractional grass-cover maps.

Source:

- https://doi.org/10.5281/zenodo.15063683
- Zhang et al. (2024), *Earth System Science Data*, 16, 1353--1381. https://doi.org/10.5194/essd-16-1353-2024

## LAI and FPAR

Monthly LAI and FPAR inputs cover 1982--2024 and are stored under:

```text
LAI/lai_1982-2024/
FPAR/fpar_1982-2024/
```

Source:

- https://www.environment.snu.ac.kr/data/longterm-lai
- Jeong et al. (2024), *Remote Sensing of Environment*, 311, 114282
- Jeong et al. (2026), *Sustained global greening driven by continuous CO2 fertilization*, Research Square preprint

## LUH2 v2h

LUH2 land-use state variables are used to derive the pasture-overlap mask.

Main input:

```text
LUH2_v2h/states.nc
```

Source:

- https://luh.umd.edu/data.shtml
- Hurtt et al. (2020), LUH2 for CMIP6
- Historical and future LUH2 input4MIPs datasets are available through ESGF

## Reproducibility

Raw datasets are not included in the repository. Dataset paths, analysis years, thresholds, and class mappings are defined in:

```text
config/config_<RUN_TAG>.yml
```