#!/bin/bash
# This script downloads the data set corresponding to Munoz-Sabater et al. 2021 
# (https://doi.org/10.5194/essd-13-4349-2021)

## E.g. ESA CCI Land Cover data: https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=overview
# Explanations for new computer: https://cds.climate.copernicus.eu/how-to-api
# see currently running requests: https://cds.climate.copernicus.eu/requests?tab=running

#**Citation:**
# Copernicus Climate Change Service, Climate Data Store, (2019): Land cover classification gridded maps from 1992 to present derived from satellite observation. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: [10.24381/cds.006f2c9a](https://doi.org/10.24381/cds.006f2c9a) (Accessed on 13-May-2025)
# Defourny, P., Lamarche, C., Brockmann, C., Boettcher, M., Bontemps, S., De Maet, T., Duveiller, G. L. Harper, K., Hartley A., Kirches, G., Moreau, I., Peylin, P., Ottlé, C., Radoux J., Van Bogaert, E., Ramoino, F., Albergel, C., and Arino, O.: Observed annual global land-use change from 1992 to 2020 three times more dynamic than reported by inventory-based statistics, in preparation, 2023.
# Harper, K. L., Lamarche, C., Hartley, A., Peylin, P., Ottlé, C., Bastrikov, V., San Martín, R., Bohnenstengel, S. I., Kirches, G., Boettcher, M., Shevchuk, R., Brockmann, C., and Defourny, P.: A 29-year time series of annual 300 m resolution plant-functional-type maps for climate models, Earth Syst. Sci. Data, 15, 1465–1499, https://doi.org/10.5194/essd-15-1465-2023, 2023. 
#**Attribution:**
#> \[Generated using/Contains modified] Copernicus Climate Change Service information \[2020]. Neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains.


##################################################################################################
# SETUP in interactive session: (2025-02-11)
# # with conda (since we need to install the library cfgrib outside of python):
# ### source: https://github.com/geco-bern/system_administration/blob/main/README.md#python-and-cuda
# ### cd ~
# ### wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# ### bash Miniconda3-latest-Linux-x86_64.sh
#
# # Then create conda environment and install needed libraries and Python packages
# conda create -n ecmwf_download
# conda activate ecmwf_download
# conda install -c conda-forge curl cfgrib h5netcdf dask xarray time cdsapi setuptools
        # conda install openssl=3.2.0 # to avoid openssl error from xarray

# OPTIONAL: # Setup R project to use that conda environment
# OPTIONAL: # go to Tools => Project Options => Python => Select conda environment "ecmwf_download"
# OPTIONAL: # reticulate::conda_list()
# OPTIONAL: # reticulate::use_condaenv(condaenv = "ecmwf_download")
# OPTIONAL: # ### reticulate::conda_install(envname = "ecmwf_download", c("curl", "cfgrib", "h5netcdf", "dask", "xarray", "time"))

# Check from within Python what is used:
# import os; print(os.system("conda env list"))
##################################################################################################


##################################################################################################
# RUN:
# in RStudio: should be taken care of by using a Project-based session and setting the conda environmnet in the Project settings
# in Terminal: ssh dash
#              cd /data/archive/landcover_defourny_2023/
#              tmux
#              conda activate ecmwf_download
#              python3 -c "$pycode" 2>&1 | tee "download_log_$(date +%Y-%m-%d_%Hh%M).txt"

# Define Python code in a variable
pycode=$(cat <<EOF
import cdsapi

client = cdsapi.Client()

dataset = "satellite-land-cover"
request = {
    "variable": "all",
    "year": ["1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015"],
    "version": ["v2_1_1"]
}
target = "land_cover_1992-2015.zip"

client.retrieve(dataset, request, target)
EOF
)

# Download zip file
cd /storage/research/giub_geco/data_2/landcover_defourny_1992-2015/
# conda run -n ecmwf_download --live-stream python3 download_landcover_defourny_2023.py 2>&1 | tee "download_log_$(date +%Y-%m-%d_%Hh%M).txt"
conda run -n ecmwf_download --live-stream python3 -c "$pycode" 2>&1 | tee "download_log_$(date +%Y-%m-%d_%Hh%M).txt"

# Extract zip file
unzip land_cover_1992-2015.zip -d data && rm land_cover_1992-2015.zip


# Run map2tidy
# TODO activate below code...
# R
# library(map2tidy)
# map2tidy::map2tidy(nclist = "data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc", outdir="data_tidy")
# map2tidy::map2tidy(nclist = "data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc", outdir="data_tidy", lonnam="lon", latnam="lat", varnam="lccs_class")
# map2tidy::map2tidy(nclist = "data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc", outdir="data_tidy", fileprefix = "C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1", lonnam="lon", latnam="lat", varnam="lccs_class", do_chunks=TRUE)

# map2tidy::map2tidy(
#   nclist     = "data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc", 
#   outdir     = "data_tidy", 
#   fileprefix = "C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1", 
#   lonnam     = "lon", 
#   latnam     = "lat", 
#   varnam     = c("lccs_class", "processed_flag","current_pixel_state","observation_count","change_count"),
#   do_chunks  = TRUE, filter_lon_between_degrees=c(1, 1.005))


# map2tidy::map2tidy(
#   nclist     = "data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc", 
#   outdir     = "data_tidy", 
#   fileprefix = "C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1", 
#   lonnam     = "lon", 
#   latnam     = "lat", 
#   varnam     = c("lccs_class"),
#   do_chunks  = TRUE, 
#   ncores     = 10) 
#   filter_lon_between_degrees = c(1, 1.005))



##################################################################################################
