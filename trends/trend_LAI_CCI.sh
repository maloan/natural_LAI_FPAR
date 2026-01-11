cd ../output/tau_0.2/masked_0p25/LAI/masked_LAI_CCI
ls LAI_masked_*.tif | head
rm *.nc
mkdir -p tmp_nc

for tif in LAI_masked_*.tif; do
  # Extract YYYYMM
  base=$(basename "$tif")
  yyyymm=${base#LAI_masked_}
  yyyymm=${yyyymm%_0p25.tif}   # e.g. 198201

  year=${yyyymm:0:4}
  mon=${yyyymm:4:2}

  nc_tmp="tmp_nc/LAI_${yyyymm}.nc"
  echo "Converting $tif → $nc_tmp"

  # 1) GeoTIFF → NetCDF
  gdal_translate -of NETCDF "$tif" "$nc_tmp"

  # 2) Set date for CDO (one day per month is fine)
  cdo -O setdate,${year}-${mon}-01 "$nc_tmp" "tmp_nc/LAI_${yyyymm}_dated.nc"
done

cd tmp_nc
cdo -O mergetime LAI_*_dated.nc ../LAI_masked_monthly_0p25.nc
cd ..

rm -r tmp_nc

# Yearly mean LAI
cdo -O yearmean LAI_masked_monthly_0p25.nc LAI_yearmean_0p25.nc

cdo -info LAI_yearmean_0p25.nc

# Yearly max LAI
cdo -O yearmax LAI_masked_monthly_0p25.nc LAI_yearmax_0p25.nc

cdo -info LAI_yearmax_0p25.nc
# Trend on annual mean
cdo -O trend LAI_yearmean_0p25.nc \
    LAI_yearmean_trend_intercept_0p25.nc \
    LAI_yearmean_trend_slope_peryear_0p25.nc

# Trend on annual max
cdo -O trend LAI_yearmax_0p25.nc \
    LAI_yearmax_trend_intercept_0p25.nc \
    LAI_yearmax_trend_slope_peryear_0p25.nc

# Yearly minimum
cdo -O yearmin LAI_masked_monthly_0p25.nc LAI_yearmin_0p25.nc

# Trend on minimum
cdo -O trend LAI_yearmin_0p25.nc \
    LAI_yearmin_trend_intercept_0p25.nc \
    LAI_yearmin_trend_slope_peryear_0p25.nc

# Amplitude = max - min
cdo -O sub LAI_yearmax_0p25.nc LAI_yearmin_0p25.nc LAI_yearamp_0p25.nc

# Trend on amplitude
cdo -O trend LAI_yearamp_0p25.nc \
    LAI_yearamp_trend_intercept_0p25.nc \
    LAI_yearamp_trend_slope_peryear_0p25.nc

mkdir -p ../../../eval/trend_LAI_CCI
mv *.nc ../../../eval/trend_LAI_CCI/.

cd ../../../../../trends
