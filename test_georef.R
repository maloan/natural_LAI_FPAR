library(terra)
library(scico)
plot(rast("output/georef_biotic/georef_biotic_LAI_0p05/LAI_202405_0p05_masked_abiotic.tif"))


f_lai <- "~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202412.nc"
lai    <- rast("~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202412.nc")

plot(
  rast("~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202112.nc")[[1]],
  col  = scico(100, palette = "batlow"),
  main = "SNU_fAPAR_v1_202312.nc"
)

f_lai <- "~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202411.nc"
lai    <- rast(f_lai[[1]])

plot(
  rast("~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202412.nc")[[1]],
  col  = scico(100, palette = "batlow"),  # try also: "roma", "batlow", "vik"
  main = "SNU_fAPAR_v1_202411.nc"
)

f_lai <- "~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202410.nc"
lai    <- rast(f_lai[[1]])

plot(
  lai,
  col  = scico(100, palette = "batlow"),  # try also: "roma", "batlow", "vik"
  main = "SNU_fAPAR_v1_202410.nc"
)

f_lai <- "~/GitHub/natural_LAI_FPAR/data-raw/FPAR/fpar_1982-2024/SNU_fAPAR_v1_202409.nc"
lai    <- rast(f_lai[[1]])

plot(
  lai,
  col  = scico(100, palette = "batlow"),  # try also: "roma", "batlow", "vik"
  main = "SNU_fAPAR_v1_202410.nc"
)

f_lai <- "~/GitHub/natural_LAI_FPAR/data-raw/LAI/lai_1982-2024/SNU_LAI_v1_202411.nc"
lai    <- rast(f_lai[[1]])

plot(
  lai,
  col  = scico(100, palette = "batlow"),  # try also: "roma", "batlow", "vik"
  main = "SNU_LAI_v1_202412.nc"
)
