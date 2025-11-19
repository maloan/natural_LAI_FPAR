## =============================================================================
# 04_luh_pasture_overlap_025.R — Build LUH pasture-overlap mask (0.25° + 0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(glue)
  library(stringr)
  library(here)
})

# --- config & refs -------------------------------------------------------------
ROOT <- here()

source(here("R", "utils.R"))
source(here("R", "io.R"))
source(here("R", "geom.R"))
source(here("R", "viz.R"))
source(here("R", "options.R"))

cfg  <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

ref005  <- rast(cfg$grids$grid_005$ref_raster)
ref025  <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

# --- params (env) --------------------------------------------------------------
GRASS_SOURCE <- toupper(Sys.getenv("GRASS_SOURCE", "CCI"))
REMAKE_QL    <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))

G_MIN  <- as.numeric(Sys.getenv("G_MIN",  "0.1"))
P_MIN  <- as.numeric(Sys.getenv("P_MIN",  "0.1"))
ALPHA  <- as.numeric(Sys.getenv("ALPHA",  "0.5"))

Y0 <- env_get_int("LUH_AVG_START", cfg$project$years$cci_start)
Y1 <- env_get_int("LUH_AVG_END",   cfg$project$years$cci_end)

stopifnot(is.finite(G_MIN),
          is.finite(P_MIN),
          is.finite(ALPHA),
          is.finite(Y0),
          is.finite(Y1))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
ql_dir  <- file.path(out_dir, "quicklooks")

dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

fname <- glue(
  "mask_luh_overlap_{GRASS_SOURCE}_Gmin{tok(G_MIN)}_",
  "Pmin{tok(P_MIN)}_alpha{tok(ALPHA)}_{Y0}-{Y1}_0p25.tif"
)

if (!exists("pal_grass", inherits = TRUE)) {
  pal_grass <- hcl.colors(64, "Greens", rev = TRUE)
}

# --- 1) grass fraction @0.05° --------------------------------------------------
grass_005 <- switch(GRASS_SOURCE,

                    "CCI" = {
                      frac_dir <- cfg$paths$cci_out_dir
                      f <- list.files(frac_dir, "ESACCI_frac_\\d{4}_0p05\\.tif$", full.names = TRUE)
                      if (!length(f)) stop("No CCI fraction files in: ", frac_dir)
                      yrs <- as.integer(str_extract(basename(f), "\\d{4}"))
                      keep <- which(yrs >= Y0 & yrs <= Y1)
                      if (!length(keep)) stop("No CCI fraction years in window ", Y0, "-", Y1)
                      stk <- rast(lapply(f[keep], function(x) rast(x)[["frac_grass"]]))
                      mean(stk, na.rm = TRUE)
                    },

                    "GLC_TEMP" = {
                      p <- file.path(cfg$paths$glc_out_dir, "glc_cat_yearstack_0p05.tif")
                      if (!file.exists(p)) stop("GLC yearstack not found: ", p)
                      s <- rast(p)
                      y <- suppressWarnings(as.integer(substr(names(s), 2, 5)))
                      keep <- which(y >= Y0 & y <= Y1)
                      if (!length(keep)) stop("No GLC years in window ", Y0, "-", Y1)

                      vec_int <- function(x) as.integer(unlist(x, use.names = FALSE))
                      grass_vals <- vec_int(cfg$glc$classes$grassland)

                      is_grass <- classify(s[[keep]], cbind(grass_vals, 1), others = 0)
                      app(is_grass, mean, na.rm = TRUE)
                    },

                    stop("Unknown GRASS_SOURCE='", GRASS_SOURCE, "'. Use CCI | GLC_TEMP.")
)

if (!compareGeom(grass_005, ref005, stopOnError = FALSE))
  grass_005 <- resample(grass_005, ref005, method = "bilinear")

names(grass_005) <- "grass"

# --- 2) aggregate grass 0.05° → 0.25° -----------------------------------------
num <- aggregate(grass_005 * area005, fact = 5, fun = function(x) sum(x, na.rm = TRUE))
den <- aggregate((!is.na(grass_005)) * area005, fact = 5, fun = function(x) sum(x, na.rm = TRUE))

grass_025 <- ifel(den == 0, NA, num / den)
grass_025 <- clamp(grass_025, 0, 1)

if (!compareGeom(grass_025, ref025, stopOnError = FALSE))
  grass_025 <- resample(grass_025, ref025, method = "near")

names(grass_025) <- "grass_025"

# --- 3) LUH pasture @0.25° (mean over window) ---------------------------------
luh_nc <- cfg$luh2$states_nc
if (!file.exists(luh_nc)) stop("LUH file not found: ", luh_nc)

v_pas <- cfg$luh2$variables$pasture

sds_names <- try(names(sds(luh_nc)), silent = TRUE)
if (inherits(sds_names, "try-error")) sds_names <- character(0)
if (!(v_pas %in% sds_names))
  stop("LUH subdataset '", v_pas, "' not found. Available: ",
       paste(sds_names, collapse=", "))

pas <- rast(luh_nc, subds = v_pas)
ty  <- suppressWarnings(as.integer(time(pas)))
if (all(is.na(ty))) stop("LUH time axis missing/NA for ", v_pas)

keep <- which(ty >= Y0 & ty <= Y1)
if (!length(keep)) stop("No LUH timesteps in window ", Y0, "-", Y1)

pasture_025 <- mean(pas[[keep]], na.rm = TRUE)

if (!compareGeom(pasture_025, ref025, stopOnError = FALSE))
  pasture_025 <- resample(pasture_025, ref025, method = "bilinear")

pasture_025 <- clamp(pasture_025, 0, 1)
names(pasture_025) <- "pasture"

# --- 4) decision at 0.25° ------------------------------------------------------
eps <- 1e-9
denom <- ifel(is.na(grass_025), NA, grass_025 + eps)
ratio <- clamp(pasture_025 / denom, 0, 1)

drop_025 <- (grass_025 >= G_MIN) &
  (pasture_025 >= P_MIN) &
  (ratio >= ALPHA)

mask_025 <- ifel(drop_025, 1L, 0L)

# --- 5) write 0.25° + 0.05° replica ------------------------------------------
wopt <- wopt_byte(speed = TRUE, na = 255L)

out025 <- file.path(out_dir, fname)

writeRaster(
  mask_025,
  out025,
  overwrite = TRUE,
  datatype = wopt$datatype,
  gdal     = wopt$gdal,
  NAflag   = wopt$NAflag
)

mask_005 <- disagg(mask_025, fact = 5, method = "near")

if (!compareGeom(mask_005, ref005, stopOnError = FALSE))
  mask_005 <- resample(mask_005, ref005, method = "near")

out005 <- file.path(out_dir,
                    sub("_0p25\\.tif$", "_0p05_rep.tif", basename(out025))
)

writeRaster(
  mask_005,
  out005,
  overwrite = TRUE,
  datatype = wopt$datatype,
  gdal     = wopt$gdal,
  NAflag   = wopt$NAflag
)

# --- 6) quicklooks -------------------------------------------------------------
ql_png <- file.path(ql_dir, sub("\\.tif$", ".png", basename(out025)))

if (REMAKE_QL || !file.exists(ql_png)) {

  png(ql_png, 1600, 750, res = 120)
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 4, 4), oma=c(2,0,2,0))

  ok <- TRUE

  tryCatch(plot(grass_025, main="Grass fraction (0.25°)",
                col=pal_grass, zlim=c(0,1)),
           error=function(e){ ok <<- FALSE })

  tryCatch(plot(pasture_025, main="Pasture fraction (0.25°)",
                col=pal_grass, zlim=c(0,1)),
           error=function(e){ ok <<- FALSE })

  tryCatch(plot(mask_025,
                main="LUH overlap mask (1=drop)",
                col=c("#f0f0f0","#d73027"),
                breaks=c(-0.5,0.5,1.5),
                legend=FALSE, axes=TRUE, box=TRUE),
           error=function(e){ ok <<- FALSE })

  if (ok) {
    legend("bottomleft",
           fill=c("#f0f0f0","#d73027"),
           legend=c("0 keep","1 drop"),
           bty="n")
  }

  if (ok) {
    mtext(
      sprintf(
        "LUH Pasture Overlap Mask (α=%s, G≥%s, P≥%s), Source=%s, %d–%d",
        tok(ALPHA), tok(G_MIN), tok(P_MIN),
        GRASS_SOURCE, Y0, Y1
      ),
      outer = TRUE, cex = 1.4, line = 0
    )
  }

  par(op)
  dev.off()
}

# --- AOI quicklooks ------------------------------------------------------------
aoi_root <- file.path(ql_dir, "aois")
dir.create(aoi_root, recursive = TRUE, showWarnings = FALSE)

aois <- cfg$aois

for (nm in names(aois)) {

  a <- aois[[nm]]
  ext_aoi <- ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)

  g_aoi <- try(crop(grass_025, ext_aoi, snap="near"), silent=TRUE)
  p_aoi <- try(crop(pasture_025, ext_aoi, snap="near"), silent=TRUE)
  m_aoi <- try(crop(mask_025, ext_aoi, snap="near"), silent=TRUE)

  if (inherits(g_aoi,"try-error") ||
      inherits(p_aoi,"try-error") ||
      inherits(m_aoi,"try-error")) next

  aoi_dir <- file.path(aoi_root, nm)
  dir.create(aoi_dir, recursive = TRUE, showWarnings = FALSE)

  ql_aoi <- file.path(
    aoi_dir,
    sprintf("luh_overlap_%s_%s_%s_%s.png",
            GRASS_SOURCE, tok(G_MIN), tok(P_MIN), nm)
  )

  if (REMAKE_QL || !file.exists(ql_aoi)) {

    png(ql_aoi, 1400, 600, res = 120)
    op <- par(mfrow=c(1,3), mar=c(3,3,3,4), oma=c(2,0,2,0))

    ok_aoi <- TRUE

    tryCatch(plot(g_aoi, main="Grass (0.25°)",
                  col=pal_grass, zlim=c(0,1)),
             error=function(e){ ok_aoi <<- FALSE })

    tryCatch(plot(p_aoi, main="Pasture (0.25°)",
                  col=pal_grass, zlim=c(0,1)),
             error=function(e){ ok_aoi <<- FALSE })

    tryCatch(plot(m_aoi,
                  main="Mask (1=drop)",
                  col=c("#f0f0f0","#d73027"),
                  breaks=c(-0.5,0.5,1.5),
                  legend=FALSE, axes=TRUE, box=TRUE),
             error=function(e){ ok_aoi <<- FALSE })

    if (ok_aoi) {
      legend("bottomleft",
             fill=c("#f0f0f0","#d73027"),
             legend=c("0 keep","1 drop"),
             bty="n")
    }

    mtext(sprintf("LUH Pasture-Overlap Mask — AOI: %s — %d–%d",
                  nm, Y0, Y1),
          outer=TRUE, line=0, cex=1.3)

    par(op)
    dev.off()
  }
}

gc()
cat(glue("
Wrote:
  - {out025}
  - {out005}
Rule: drop if grass≥{G_MIN} & pasture≥{P_MIN} & pasture/grass≥{ALPHA}
Source grass={GRASS_SOURCE}; window={Y0}-{Y1}; semantics: 1=drop, 0=keep, 255=NA
"))
