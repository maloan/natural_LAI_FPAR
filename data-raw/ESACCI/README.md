# ESA-CCI and C3S annual land cover

The workflow reads the 1992--2022 annual 300 m land-cover ZIP archives from `ESACCI_1992-2022/`. 

- `download_landcover_1992_2015.py` downloads ESA-CCI v2.0.7cds.
- `download_landcover_2016_2022.py` downloads C3S v2.1.1.
- `R/12_make_lc025_fractions.R` temporarily extracts one annual NetCDF at a time, validates it, and writes the versioned 22-layer 0.25-degree fractions.

The data use the 22-class global UN-LCCS legend. The workflow preserves these 22 primary classes during spatial aggregation, including codes 30 and 40 as separate agricultural mosaic classes. It combines the fractions into the 19 analysis groups defined in `R/helpers/landcover_classification.R` before selecting the dominant class.

Dataset record and DOI:
https://doi.org/10.24381/cds.006f2c9a
