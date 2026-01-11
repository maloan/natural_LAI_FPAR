## =============================================================================
# 07_abiotic_static_from_cci.R — Build abiotic mask for ONE year (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
# source(here("R","helpers", "viz.R"))
#
cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.6)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
cci_dir <- cfg$paths$cci_dir

TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE <- as.numeric(Sys.getenv("TAU_ICE", "0.05"))
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
YEAR <- as.integer(Sys.getenv("YEAR", "2007"))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

tauW_tok <- gsub("\\.", "p", sprintf("%.2f", TAU_WATER))
tauI_tok <- gsub("\\.", "p", sprintf("%.2f", TAU_ICE))

out_tif <- file.path(out_dir, sprintf(
  "mask_abiotic_CCI_%d_tauW%s_tauI%s_0p05.tif", YEAR, tauW_tok, tauI_tok
))

if (SKIP_EXISTING && file.exists(out_tif)) {
  message("✓ Abiotic mask already exists — skipping: ", out_tif)
  invisible(gc())
  quit(save = "no")
}

ESACCI <- cfg$esa_cci$classes
vals_water <- as.integer(unlist(ESACCI$water, use.names = FALSE))
vals_ice <- as.integer(unlist(ESACCI$snow_ice, use.names = FALSE))
nodata_vals <- unique(c(as.integer(
  unlist(ESACCI$nodata, use.names = FALSE)
), 255L))

# --- pick input for YEAR (prefer C3S if multiple) -----------------------------
files <- list.files(cci_dir, pattern = "\\.tif$", full.names = TRUE)
if (!length(files)) stop("No CCI GeoTIFFs found in: ", cci_dir, call. = FALSE)

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
if (!length(cand)) {
  stop("No CCI file found for YEAR=", YEAR,
    " in: ", cci_dir,
    call. = FALSE
  )
}

cand_rank <- vapply(cand, rank_src, integer(1))
in_file <- cand[which.max(cand_rank)]

message("→ Processing ", basename(in_file), " (YEAR=", YEAR, ")")

r <- rast(in_file)
if (is.na(crs(r))) crs(r) <- crs(ref005)

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

abi_mask <- ifel((pW >= TAU_WATER) | (pI >= TAU_ICE), 1L, 0L)
names(abi_mask) <- "abiotic_drop"

writeRaster(
  abi_mask,
  out_tif,
  overwrite = TRUE,
  wopt = wopt_byte(FALSE, na = 255L)
)

quicklook_mask_all_aois(
  mask = abi_mask,
  title = sprintf(
    "Abiotic mask — %d (pW≥%.2f ∨ pI≥%.2f)", YEAR, TAU_WATER, TAU_ICE
  ),
  tag = sprintf("abiotic_%d", YEAR),
  cfg = cfg,
  ql_root = ql_dir,
  down = 4L,
  include_global = TRUE,
  drop_global_key = FALSE
)

p_drop <- tryCatch(
  global(
    abi_mask, "mean",
    na.rm = TRUE
  )[1, 1],
  error = function(e) NA_real_
)

cat(sprintf(
  "\nWritten: %s\nYear: %d\nDrop thresholds — water=%.2f,
  ice=%.2f\nGlobal drop fraction: %.4f\nQuicklooks: %s\n",
  out_tif, YEAR, TAU_WATER, TAU_ICE, p_drop, ql_dir
))

gc()
