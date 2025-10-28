

## =============================================================================
# 08_agg_0p25.R — Area-weighted aggregation of masked LAI/FPAR from 0.05° → 0.25°

suppressPackageStartupMessages({
  library(terra)
  library(stringr)
  library(glue)
})

# --- config & refs -------------------------------------------------------------
# Root dir
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
),
winslash = "/",
mustWork = FALSE)
# Other source files
source(file.path(ROOT, "R", "00_utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "FALSE"))
OVERWRITE     <- as.logical(Sys.getenv("OVERWRITE", "TRUE"))
REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))

# --- choose variable -----------------------------------------------------------
VAR  <- toupper(Sys.getenv("VAR", "FPAR"))   # FPAR | FPAR
MASK <- toupper(Sys.getenv("MASK", "GLC"))   # CCI | GLC

# --- paths, patterns, visualization ranges ------------------------------------
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)

if (VAR == "LAI") {
  in_dir  <- if (MASK == "CCI") {
    cfg$paths$masked_lai_cci_005_dir
  } else {
    cfg$paths$masked_lai_glc_005_dir
  }
  out_dir <- if (MASK == "CCI") {
    cfg$paths$masked_lai_cci_025_dir
  } else {
    cfg$paths$masked_lai_glc_025_dir
  }
  patt     <- "^LAI_.*_\\d{6}_0p05_masked\\.tif$"
  out_name <- "LAI_masked_{ym}_0p25.tif"
  ql_title <- "LAI"
  zlim     <- c(cfg$variables$lai$clamp$min, cfg$variables$lai$clamp$max)
} else {
  in_dir  <- if (MASK == "CCI") {
    cfg$paths$masked_fpar_cci_005_dir
  } else {
    cfg$paths$masked_fpar_glc_005_dir
  }
  out_dir <- if (MASK == "CCI") {
    cfg$paths$masked_fpar_cci_025_dir
  } else {
    cfg$paths$masked_fpar_glc_025_dir
  }
  patt     <- "^FPAR_.*_\\d{6}_0p05_masked\\.tif$"
  out_name <- "FPAR_masked_{ym}_0p25.tif"
  ql_title <- "FPAR"
  zlim     <- c(cfg$variables$fpar$clamp$min,
                cfg$variables$fpar$clamp$max)
}

qdir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(qdir, TRUE, showWarnings = FALSE)

# --- latitude weights on the 0.05° grid ---------------------------------------
area005 <- rast(cfg$grids$grid_005$area_raster)

# --- discover inputs -----------------------------------------------------------
files <- sort(list.files(in_dir, patt, full.names = TRUE))
if (!length(files))
  stop("No ", VAR, " inputs found in: ", in_dir)

# --- pretty quicklook (global, overlays) --------------------------------------
quicklook <- function(r025,
                      ym,
                      down = 1L,
                      title = ql_title,
                      zlim = zlim) {
  rr <- if (down > 1L) {
    aggregate(r025, down, mean, na.rm = TRUE)
  } else {
    r025
  }
  png(
    file.path(qdir, sprintf("quicklook_%s_0p25_%s.png", title, ym)),
    width = 1400,
    height = 700,
    res = 120
  )
  op <- par(oma = c(0, 0, 2.2, 0), mar = c(3, 3, 3, 6))
  terra::plot(
    rr,
    main   = sprintf("%s 0.25° %s", title, ym),
    col    = pal_green(64),
    colNA  = col_na,
    zlim   = zlim,
    axes   = TRUE,
    legend = TRUE
  )
  .add_overlays(rr)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)", side = 2, line = 2)
  mtext(
    "Area-weighted aggregation to 0.25°",
    side = 3,
    outer = TRUE,
    cex = 1.05
  )
  par(op)
  dev.off()
}

# --- aggregation loop ----------------------------------------------------------
for (f in files) {
  ym  <- str_extract(basename(f), "\\d{6}")
  out <- file.path(out_dir, glue(out_name))

  if (!(SKIP_EXISTING && file.exists(out)) || OVERWRITE) {
    r <- rast(f)  # 0.05° masked monthly field

    # Snap to 0.05° template if needed (continuous → bilinear is fine)
    if (!compareGeom(r, ref005, stopOnError = FALSE)) {
      r <- resample(r, ref005, method = "bilinear")
    }
    # Area-weighted mean over 5×5 blocks
    num <- aggregate(
      r * area005,
      fact = 5,
      fun = function(x, ...)
        sum(x, na.rm = TRUE)
    )
    den <- aggregate((!is.na(r)) * area005,
                     fact = 5,
                     fun = function(x, ...)
                       sum(x, na.rm = TRUE)
    )
    r025 <- ifel(den == 0, NA, num / den)

    # Snap to canonical 0.25° template (index-safe)
    if (!compareGeom(r025, ref025, stopOnError = FALSE)) {
      r025 <- resample(r025, ref025, method = "near")
    }

    writeRaster(
      r025,
      out,
      overwrite = TRUE,
      gdal   = gdal_wopt("FLT4S")$gdal,
      NAflag = -9999
    )
  } else {
    r025 <- rast(out)  # load existing for QL
  }

  # Lightweight QA thumbs (Jan/Jul)
  if (substr(ym, 5, 6) %in% c("01", "07")) {
    ql_png <- file.path(qdir, sprintf("quicklook_%s_0p25_%s.png", ql_title, ym))
    if (REMAKE_QL || !(SKIP_EXISTING && file.exists(ql_png))) {
      quicklook(r025,
                ym,
                down = 1L,
                title = ql_title,
                zlim = zlim)
    }
  }
}

message("Aggregated monthly ", VAR, " written to: ", out_dir)
