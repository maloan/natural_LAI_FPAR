cd ../output/tau_0.2/masked_0p25/FPAR/masked_FPAR_GLC
ls FPAR_masked_*.tif | head
rm *.nc
mkdir -p tmp_nc

for tif in FPAR_masked_*.tif; do
  # Extract YYYYMM
  base=$(basename "$tif")
  yyyymm=${base#FPAR_masked_}
  yyyymm=${yyyymm%_0p25.tif}   # e.g. 198201

  year=${yyyymm:0:4}
  mon=${yyyymm:4:2}

  nc_tmp="tmp_nc/FPAR_${yyyymm}.nc"
  echo "Converting $tif → $nc_tmp"

  # 1) GeoTIFF → NetCDF
  gdal_translate -of NETCDF "$tif" "$nc_tmp"

  # 2) Set date for CDO (one day per month is fine)
  cdo -O setdate,${year}-${mon}-01 "$nc_tmp" "tmp_nc/FPAR_${yyyymm}_dated.nc"
done

cd tmp_nc
cdo -O mergetime FPAR_*_dated.nc ../../FPAR_masked_monthly_0p25.nc
cd ../..

rm -r tmp_nc

# Yearly mean FPAR
cdo -O yearmean FPAR_masked_monthly_0p25.nc FPAR_yearmean_0p25.nc

cdo -info FPAR_yearmean_0p25.nc

# Yearly max FPAR
cdo -O yearmax FPAR_masked_monthly_0p25.nc FPAR_yearmax_0p25.nc

cdo -info FPAR_yearmax_0p25.nc
# Trend on annual mean
cdo -O trend FPAR_yearmean_0p25.nc \
    FPAR_yearmean_trend_intercept_0p25.nc \
    FPAR_yearmean_trend_slope_peryear_0p25.nc

# Trend on annual max
cdo -O trend FPAR_yearmax_0p25.nc \
    FPAR_yearmax_trend_intercept_0p25.nc \
    FPAR_yearmax_trend_slope_peryear_0p25.nc

# Yearly minimum
cdo -O yearmin FPAR_masked_monthly_0p25.nc FPAR_yearmin_0p25.nc

# Trend on minimum
cdo -O trend FPAR_yearmin_0p25.nc \
    FPAR_yearmin_trend_intercept_0p25.nc \
    FPAR_yearmin_trend_slope_peryear_0p25.nc

# Amplitude = max - min
cdo -O sub FPAR_yearmax_0p25.nc FPAR_yearmin_0p25.nc FPAR_yearamp_0p25.nc

# Trend on amplitude
cdo -O trend FPAR_yearamp_0p25.nc \
    FPAR_yearamp_trend_intercept_0p25.nc \
    FPAR_yearamp_trend_slope_peryear_0p25.nc

mkdir -p ../../../../eval/trend_FPAR_GLC
mv *.nc ../../../../eval/trend_FPAR_GLC/.

cd ../../../../../../
