# Raw Datasets (data-raw)

This folder contains external source datasets used by the workflow. The pipeline reads from data-raw and writes derived products to data and output.

## Folder layout

```text
data-raw/
├── ESACCI/
├── GLC_FCS30D/
├── LAI/
├── FPAR/
└── LUH2_v2h/
```

## Datasets

### ESA-CCI / C3S land cover (ESACCI)

Used to build fractional land-cover layers and CCI-based masks.

Included content:

- Download scripts: `download_landcover_1992_2015.py` and
  `download_landcover_2016_2022.py`.
- Downloaded ZIP archives and extracted annual NetCDF maps in
  `ESACCI/ESACCI_1992-2022/`.

The preprocessing reads the `lccs_class` variable directly from the extracted NetCDF files. The expected versions are v2.0.7cds for 1992--2015 and v2.1.1 for 2016--2022.

Data source:
- https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=download
- Copernicus Climate Change Service, Climate Data Store, (2019): Land cover classification gridded maps from 1992 to present derived from satellite observation. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: 10.24381/cds.006f2c9a (Accessed on 01.07.2026)

### GLC_FCS30D (GLC_FCS30D)

Used for the GLC-based masking branch.

Included content:

- The 36 original GLC_FCS30D v2 ZIP archives are stored in `data-raw/GLC_FCS30D/archives/` (approximately 135.8 GB in total). This local raw-data directory is excluded from version control.
- Global categorical 0.05 degree GeoTIFF maps named `GLC_FCS30D_mode_0p05_<year>.tif`.
- A global fractional-grass stack named `GLC_FCS30D_grass_fraction_0p05_1985_2022.tif`.


Data source:
- https://doi.org/10.5281/zenodo.15063683
- Zhang, X., Zhao, T., Xu, H., Liu, W., Wang, J., Chen, X., and Liu, L.: GLC_FCS30D: the first global 30 m land-cover dynamics monitoring product with a fine classification system for the period from 1985 to 2022 generated using dense-time-series Landsat imagery and the continuous change-detection method, Earth Syst. Sci. Data, 16, 1353–1381, https://doi.org/10.5194/essd-16-1353-2024, 2024.
- Liangyun Liu, Xiao Zhang, & Zhehua Li. (2025). GLC_FCS30D: the first global 30-m land-cover dynamic monitoring product with fine classification system from 1985 to 2022 (Version v2) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.15063683

### LAI and FPAR monthly inputs (LAI, FPAR)

Monthly NetCDF time series used before georeferencing, masking, and aggregation.

Included content:

- LAI/lai_1982-2024/
- FPAR/fpar_1982-2024/

Data source:

- https://www.environment.snu.ac.kr/data/longterm-lai
- Jeong, S., Ryu, Y., Gentine, P., Lian, X., Fang, J., Li, X., Dechant, B., Kong, J., Choi, W., Jiang, C., Keenan, T. F., Harrison, S. P., & Prentice, I. C. (2024). Persistent global greening over the last four decades using novel long-term vegetation index data with enhanced temporal consistency. Remote Sensing of Environment, 311, 114282. [HTML] (2024.09)​
- Sungchan Jeong, Youngryel Ryu, Pierre Gentine et al. Sustained global greening driven by continuous CO2 fertilization, 21 January 2026, PREPRINT (Version 1) available at Research Square [HTML]​

### LUH2 v2h (LUH2_v2h)

Land-use state variables used for the pasture--grass overlap mask.

Included content:

- states.nc
- supporting LUH2 documentation files.

Data source:

- https://luh.umd.edu/data.shtml

- Hurtt, G. C., L. Chini, R. Sahajpal, S. Frolking, B. L. Bodirsky, K. Calvin, J. C. Doelman, J. Fisk, S. Fujimori, K. K. Goldewijk, T. Hasegawa, P. Havlik, A. Heinimann, F. Humpenöder, J. Jungclaus, Jed Kaplan, J. Kennedy, T. Kristzin, D. Lawrence, P. Lawrence, L. Ma, O. Mertz, J. Pongratz, A. Popp, B. Poulter, K. Riahi, E. Shevliakova, E. Stehfest, P. Thornton, F. N. Tubiello, D. P. van Vuuren, X. Zhang (2020). Harmonization of Global Land-Use Change and Management for the Period 850-2100 (LUH2) for CMIP6.Geoscientifc Model Development Discussions. https://doi.org/10.5194/gmd-2019-360
- Hurtt, G. C., Chini, L., Sahajpal, R., Frolking, S., Bodirsky, B. L., Calvin, K., Doelman, J., Fisk, J., Fujimori, S., Goldewijk, K. K., Hasegawa, T., Havlik, P., Heinimann, A., Humpenöder, F., Jungclaus, J., Kaplan, J., Krisztin, T., Lawrence, D., Lawrence, P., Mertz, O., Pongratz, J., Popp, A., Riahi, K., Shevliakova, E., Stehfest, E., Thornton, P., van Vuuren, D., Zhang, X. (2019). Harmonization of Global Land Use Change and Management for the Period 2015-2300. Version 20190529. Earth System Grid Federation. https://doi.org/10.22033/ESGF/input4MIPs.10468
- Hurtt, G. C., Chini, L., Sahajpal, R., Frolking, S., Bodirsky, B. L., Calvin, K., Doelman, J., Fisk, J., Fujimori, S., Goldewijk, K. K., Hasegawa, T., Havlik, P., Heinimann, A., Humpenöder, F., Jungclaus, J., Kaplan, J., Krisztin, T., Lawrence, D., Lawrence, P., Mertz, O., Pongratz, J., Popp, A., Riahi, K., Shevliakova, E., Stehfest, E., Thornton, P., van Vuuren, D., Zhang, X. (2019). Harmonization of Global Land Use Change and Management for the Period 850-2015. Version 20190529. Earth System Grid Federation. https://doi.org/10.22033/ESGF/input4MIPs.10454
 
---
## Reproducibility notes

- data-raw contents are not tracked by git.
- Paths, year windows, and class mappings are defined in `config/config_<RUN_TAG>.yml`.
