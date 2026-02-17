## =============================================================================
# 07_abiotic_static_from_cci.R — Build abiotic mask for ONE year (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.6)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
cci_dir <- cfg$paths$cci_dir

TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE <- as.numeric(Sys.getenv("TAU_ICE", "0.05"))
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "FALSE"))
YEAR <- as.integer(Sys.getenv("YEAR", "2007"))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tauW_tok <- gsub("\\.", "p", sprintf("%.2f", TAU_WATER))
tauI_tok <- gsub("\\.", "p", sprintf("%.2f", TAU_ICE))

out_tif <- file.path(out_dir, sprintf(
  "mask_abiotic_CCI_%d_tauW%s_tauI%s_0p05.tif", YEAR, tauW_tok, tauI_tok
))

if (SKIP_EXISTING && file.exists(out_tif)) {
  message("✓ Abiotic mask already exists — skipping: ", out_tif)
  quit(save = "no")
}

ESACCI <- cfg$esa_cci$classes
vals_water <- as.integer(unlist(ESACCI$water, use.names = FALSE))
vals_ice <- as.integer(unlist(ESACCI$snow_ice, use.names = FALSE))
nodata_vals <- unique(c(as.integer(
  unlist(ESACCI$nodata, use.names = FALSE)
), 255L))

# ---  input for YEAR  -----------------------------
files <- list.files(cci_dir, pattern = "\\.tif$", full.names = TRUE)

get_year <- function(p) {
  m <- regexpr("(19|20)\\d{2}", basename(p), perl = TRUE)
  if (m[1] < 0) {
    return(NA_integer_)
  }
  as.integer(substr(basename(p), m[1], m[1] + attr(m, "match.length") - 1))
}
rank_src <- function(p) if (grepl("^C3S", basename(p))) 2L else 1L

yrs <- vapply(files, get_year, integer(1))
cand <- files[!is.na(yrs) & yrs == YEAR]

cand_rank <- vapply(cand, rank_src, integer(1))
in_file <- cand[which.max(cand_rank)]

message("→ Processing ", basename(in_file), " (YEAR=", YEAR, ")")

r <- rast(in_file)
if (is.na(crs(r))) {
  crs(r) <- crs(ref005)
}

r <- terra::subst(r, nodata_vals, NA)

# water / ice fractions at 0.05° (single year)
pW <- resample(classify(
  r, cbind(vals_water, 1),
  others = 0
), ref005, method = "average")
pI <- resample(classify(
  r, cbind(vals_ice, 1),
  others = 0
), ref005, method = "average")

water_drop <- ifel(pW >= TAU_WATER, 1L, 0L)
ice_drop <- ifel(pI >= TAU_ICE, 1L, 0L)
both_drop <- ifel((water_drop == 1L) & (ice_drop == 1L), 1L, 0L)
abi_mask <- ifel((water_drop == 1L) | (ice_drop == 1L), 1L, 0L)

# write components
writeRaster(c(water_drop, ice_drop, both_drop, abi_mask),
  file.path(out_dir, sprintf("abiotic_components_%d_0p05.tif", YEAR)),
  overwrite = TRUE, wopt = wopt_byte(FALSE, na = 255L)
)

abi_mask <- ifel((pW >= TAU_WATER) | (pI >= TAU_ICE), 1L, 0L)
names(abi_mask) <- "abiotic_drop"

writeRaster(
  abi_mask,
  out_tif,
  overwrite = TRUE,
  wopt = wopt_byte(FALSE, na = 255L)
)

p_drop <- tryCatch(
  global(
    abi_mask, "mean",
    na.rm = TRUE
  )[1, 1],
  error = function(e) NA_real_
)

gc()
