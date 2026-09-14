# =============================================================================
# 12_chapter2_inputs_0p5.R — Build Chapter 2 fAPAR and mask inputs at 0.5°
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(terra)
})

source(here("R", "helpers", "io.R"))

terraOptions(progress = 1, memfrac = 0.25)

ref005 <- rast(here("src", "ref_0p05.nc"))
area005 <- rast(here("src", "area_0p05_km2.nc"))
area005 <- align_to_template(area005, ref005, method = "bilinear")

out_dir <- here("output", "chapter2")
tmp_dir <- file.path(out_dir, "tmp_fpar_unmasked_0p5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# Unified CCI alpha=0.1 and pasture mask. Values are 1=drop and 0=keep.
cci_file <- here(
  "output", "alpha_0.1", "masks", "mask_cci",
  "mask_used_frac_fused_alpha0p10_k3_1992-2022_0p05.tif"
)
pasture_file <- here(
  "output", "alpha_0.1", "masks", "mask_luh_overlap",
  "mask_luh_overlap_CCI_Gmin0p10_Pmin0p10_beta0p50_1992-2015_0p05_rep.tif"
)

if (!file.exists(cci_file) || !file.exists(pasture_file)) {
  stop("Missing CCI or pasture mask.")
}

cci005 <- align_to_template(rast(cci_file), ref005, method = "near")
pasture005 <- align_to_template(rast(pasture_file), ref005, method = "near")

cci005 <- ifel(is.na(cci005), 0L, cci005 >= 1)
pasture005 <- ifel(is.na(pasture005), 0L, pasture005 >= 1)
mask005 <- ifel(cci005 == 1 | pasture005 == 1, 1L, 0L)

# Preserve the excluded-area fraction before creating the conservative binary mask.
mask_area050 <- aggregate(mask005 * area005, fact = 10, fun = "sum", na.rm = TRUE)
total_area050 <- aggregate(area005, fact = 10, fun = "sum", na.rm = TRUE)
excluded_fraction050 <- mask_area050 / total_area050
excluded_fraction050 <- clamp(excluded_fraction050, 0, 1)
names(excluded_fraction050) <- "excluded_fraction"

mask050 <- ifel(excluded_fraction050 > 0, 1L, 0L)
names(mask050) <- "mask"

mask_file <- file.path(out_dir, "mask_CCI_alpha0p10_pasture_any_0p5.nc")
fraction_file <- file.path(
  out_dir,
  "mask_CCI_alpha0p10_pasture_excluded_fraction_0p5.tif"
)

writeCDF(
  mask050,
  mask_file,
  overwrite = TRUE,
  varname = "mask",
  longname = "CCI alpha 0.1 or pasture exclusion mask; 1=drop, 0=keep",
  unit = "1",
  prec = "short",
  missval = -9999
)
writeRaster(
  excluded_fraction050,
  fraction_file,
  overwrite = TRUE,
  wopt = wopt_f32(FALSE)
)

# Area-weight unmasked monthly fAPAR from 0.05° to 0.5°.
months <- format(
  seq(as.Date("1982-01-01"), as.Date("2024-12-01"), by = "month"),
  "%Y%m"
)
input_dir <- here("data", "georef", "georef_fpar_0p05")
input_files <- file.path(input_dir, sprintf("FPAR_%s_0p05.tif", months))
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing unmasked fAPAR inputs:\n", paste(missing_files, collapse = "\n"))
}

monthly_files <- file.path(tmp_dir, sprintf("FPAR_unmasked_%s_0p5.tif", months))

for (i in seq_along(input_files)) {
  if (file.exists(monthly_files[i])) {
    next
  }

  message("Aggregating ", months[i], " (", i, "/", length(months), ")")
  fpar005 <- rast(input_files[i])
  compareGeom(fpar005, ref005, stopOnError = TRUE)

  numerator <- aggregate(
    fpar005 * area005,
    fact = 10,
    fun = "sum",
    na.rm = TRUE
  )
  denominator <- aggregate(
    ifel(is.finite(fpar005), area005, 0),
    fact = 10,
    fun = "sum",
    na.rm = TRUE
  )
  fpar050 <- ifel(denominator > 0, numerator / denominator, NA)
  names(fpar050) <- "fpar"

  writeRaster(
    fpar050,
    monthly_files[i],
    overwrite = TRUE,
    wopt = wopt_f32(FALSE)
  )
}

fpar_stack <- rast(monthly_files)
names(fpar_stack) <- paste0("fpar_", months)
time(fpar_stack) <- as.Date(paste0(months, "01"), format = "%Y%m%d")

fpar_file <- file.path(out_dir, "fpar_unmasked_0p5_monthly_1982-2024.nc")
writeCDF(
  fpar_stack,
  fpar_file,
  overwrite = TRUE,
  varname = "fpar",
  longname = "Unmasked fraction of absorbed photosynthetically active radiation",
  unit = "1",
  prec = "float",
  compression = 4,
  missval = -9999
)

check <- rast(fpar_file)
if (nlyr(check) != 516L || !all(res(check) == 0.5)) {
  stop("The Chapter 2 fAPAR NetCDF failed its geometry or time-layer check.")
}

value_range <- global(check, range, na.rm = TRUE)
if (min(value_range[, 1]) < 0 || max(value_range[, 2]) > 1) {
  stop("The Chapter 2 fAPAR NetCDF contains values outside 0--1.")
}

unlink(tmp_dir, recursive = TRUE)

message("Wrote ", fpar_file)
message("Wrote ", mask_file)
message("Wrote ", fraction_file)
