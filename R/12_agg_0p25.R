## =============================================================================
# 12_agg_0p25.R — Area-weighted aggregation of masked LAI/FPAR (0.05° → 0.25°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

# --- repo helpers --------------------------------------------------------------
source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# ------------------------------------------------------------------------------
# Env
# ------------------------------------------------------------------------------
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "FALSE"))
OVERWRITE <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

VAR <- toupper(Sys.getenv("VAR", "FPAR")) # LAI|FPAR
MASK <- toupper(Sys.getenv("MASK", "CCI")) # CCI|GLC

# ------------------------------------------------------------------------------
# Refs and weights
# ------------------------------------------------------------------------------
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

if (!compareGeom(area005, ref005, stopOnError = FALSE)) {
  area005 <- resample(area005, ref005, method = "bilinear")
}

# ------------------------------------------------------------------------------
# dirs
# ------------------------------------------------------------------------------
key_in <- sprintf("masked_%s_%s_0p05_dir", tolower(VAR), tolower(MASK))
key_out <- sprintf("masked_%s_%s_0p25_dir", tolower(VAR), tolower(MASK))

in_dir <- cfg$paths[[key_in]]
out_dir <- cfg$paths[[key_out]]

stopifnot(is.character(in_dir), length(in_dir) == 1, nzchar(in_dir))
stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))

patt <- sprintf("^%s_.*_\\d{6}_0p05_masked\\.tif$", VAR)
ql_title <- VAR

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
qdir <- file.path(out_dir, "quicklooks")
dir.create(qdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------
stopifnot(dir.exists(in_dir))
files <- sort(list.files(in_dir, pattern = patt, full.names = TRUE))

# ------------------------------------------------------------------------------
# Quicklook helper
# ------------------------------------------------------------------------------
quicklook_025 <- function(r025, ym, title = ql_title) {
  out_png <- file.path(qdir, sprintf("quicklook_%s_0p25_%s.png", title, ym))
  png(out_png, width = 1400, height = 700, res = 120)
  op <- par(oma = c(0, 0, 2.2, 0), mar = c(3, 3, 3, 6))
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  terra::plot(
    rr,
    main   = sprintf("%s 0.25° %s", title, ym),
    col    = pal_green(64),
    colNA  = col_na,
    axes   = TRUE,
    legend = TRUE
  )
  .add_overlays(rr)
  mtext("Longitude (°E)", 1, line = 2)
  mtext("Latitude (°N)", 2, line = 2)
  mtext("Area-weighted aggregation to 0.25°", 3, outer = TRUE, cex = 1.05)

  invisible(out_png)
}

# ------------------------------------------------------------------------------
# Aggregation loop
# ------------------------------------------------------------------------------
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_masked_%s_0p25.tif", VAR, ym))
  do_write <- OVERWRITE || !file.exists(out) || !SKIP_EXISTING
  do_ql <- (substr(ym, 5, 6) %in% c("01", "07")) &&
    (REMAKE_QL || !file.exists(file.path(
      qdir, sprintf("quicklook_%s_0p25_%s.png", ql_title, ym)
    )))

  if (!do_write && !do_ql) {
    next
  }

  if (do_write) {
    r <- rast(f)

    if (!compareGeom(r, ref005, stopOnError = FALSE)) {
      r <- resample(r, ref005, method = "bilinear")
    }

    # numerator/denominator (area-weighted mean, handling NA)
    num <- aggregate(r * area005, fact = 5, fun = "sum", na.rm = TRUE)
    den <- aggregate((!is.na(r)) * area005, fact = 5, fun = "sum", na.rm = TRUE)

    r025 <- ifel(den == 0, NA, num / den)

    if (!compareGeom(r025, ref025, stopOnError = FALSE)) {
      r025 <- resample(r025, ref025, method = "near")
    }
    wopt <- wopt_f32(opts$SPEED_OVER_SIZE)$wopt
    writeRaster(
      r025,
      out,
      overwrite = TRUE,
      wopt = wopt
    )
  } else {
    r025 <- rast(out)
    if (!compareGeom(r025, ref025, stopOnError = FALSE)) {
      r025 <- resample(r025, ref025, method = "near")
    }
  }

  if (do_ql) {
    quicklook_025(r025, ym, title = ql_title)
  }

  gc()
}

gc()
message("Aggregated monthly ", VAR, " written to: ", out_dir)
