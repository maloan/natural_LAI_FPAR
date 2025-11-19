#!/usr/bin/env Rscript
## =============================================================================
# 23_master_trend_plots.R — Full diagnostics for LAI/FPAR masked/unmasked trends
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(Kendall)
  library(here)
})

# ---------------------------------------------------------------------
# Root and output using here()
# ---------------------------------------------------------------------
ROOT <- here()
OUTP <- here("analysis", "trend_plots")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Coastline
# ---------------------------------------------------------------------
coast_path <- here("src", "ne_110m_coastline.gpkg")
coast <- if (file.exists(coast_path)) read_sf(coast_path) else NULL

# Palettes
pal_seq <- function(n=100) hcl.colors(n, "batlow")
pal_div <- function(n=100) hcl.colors(n, "Blue-Red 3")

# ---------------------------------------------------------------------
# Helper: global plot
# ---------------------------------------------------------------------
plot_global_map <- function(r, title, outfile, lim = NULL, diverging = FALSE) {
  png(outfile, width = 1800, height = 900, res = 150)
  par(mai = c(0.6,0.8,0.7,0.8), mgp = c(2,0.6,0), tck = -0.01)

  if (is.null(lim)) lim <- range(values(r), na.rm = TRUE)

  pal <- if (diverging) pal_div(100) else pal_seq(100)

  plot(r, col = pal, range = lim, axes = FALSE, main = title)
  if (!is.null(coast)) plot(coast$geom, add = TRUE, lwd = 0.3)

  dev.off()
}

# ---------------------------------------------------------------------
# Find all slope files
# ---------------------------------------------------------------------
files <- list.files(
  ROOT,
  pattern = "trend_slope_yearmean_decade.*\\.nc$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) stop("No trend slope files found.")

# ---------------------------------------------------------------------
# Parse metadata from filename
# ---------------------------------------------------------------------
parse_info <- function(f){
  var  <- if (grepl("FPAR", f)) "FPAR" else "LAI"
  mask <- "UNMASKED"
  if (grepl("CCI", f)) mask <- "CCI"
  if (grepl("GLC", f)) mask <- "GLC"
  if (grepl("georef", f)) mask <- "GEOREF"

  tau <- "NA"
  mt <- regmatches(f, regexpr("tau_0\\.[0-9]+", f))
  if (length(mt) == 1) tau <- mt

  list(var = var, mask = mask, tau = tau)
}

# ---------------------------------------------------------------------
# Load rasters
# ---------------------------------------------------------------------
slopes <- lapply(files, function(f){
  info <- parse_info(f)
  r <- rast(f)
  names(r) <- paste(info$var, info$mask, info$tau, sep = "_")
  list(info = info, file = f, r = r)
})

# ---------------------------------------------------------------------
# 1. Global maps
# ---------------------------------------------------------------------
for (s in slopes){
  info <- s$info
  outdir <- file.path(OUTP, info$var, info$mask, info$tau)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  outpng <- file.path(outdir,
                      sprintf("%s_%s_%s_yearmean_trend.png",
                              info$var, info$mask, info$tau))

  plot_global_map(
    s$r,
    sprintf("%s — %s — %s\nTrend (yearmean per decade)",
            info$var, info$mask, info$tau),
    outpng,
    lim = c(-0.05, 0.05),
    diverging = TRUE
  )
}

# ---------------------------------------------------------------------
# 2. Masked vs unmasked
# ---------------------------------------------------------------------
unmasked <- slopes[sapply(slopes, function(s) s$info$mask == "UNMASKED")]

for (um in unmasked){
  var <- um$info$var

  masked <- slopes[sapply(slopes, function(s)
    s$info$var == var && s$info$mask != "UNMASKED")]

  for (mk in masked){
    diff_r <- mk$r - um$r

    outdir <- file.path(OUTP, "diff_masked_unmasked", var)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    outpng <- file.path(
      outdir,
      sprintf("%s_%s_%s_minus_unmasked.png", var, mk$info$mask, mk$info$tau)
    )

    plot_global_map(
      diff_r,
      sprintf("%s: Masked (%s %s) − Unmasked\nTrend difference per decade",
              var, mk$info$mask, mk$info$tau),
      outpng,
      lim = c(-0.05, 0.05),
      diverging = TRUE
    )
  }
}

# ---------------------------------------------------------------------
# 3. CCI τ comparisons
# ---------------------------------------------------------------------
cci <- slopes[sapply(slopes, function(s) s$info$mask == "CCI")]

for (var in c("LAI","FPAR")){
  ccivar <- cci[sapply(cci, function(s) s$info$var == var)]
  taus <- unique(sapply(ccivar, function(s) s$info$tau))
  if (length(taus) < 2) next

  for (i in seq_along(ccivar)){
    for (j in seq_along(ccivar)){
      if (i >= j) next

      r1 <- ccivar[[i]]$r
      r2 <- ccivar[[j]]$r
      t1 <- ccivar[[i]]$info$tau
      t2 <- ccivar[[j]]$info$tau

      diff_r <- r1 - r2

      outdir <- file.path(OUTP, "cci_tau_comparison", var)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

      outpng <- file.path(
        outdir,
        sprintf("%s_CCI_%s_minus_%s.png", var, t1, t2)
      )

      plot_global_map(
        diff_r,
        sprintf("%s CCI: %s - %s\nTrend difference (per decade)", var, t1, t2),
        outpng,
        lim = c(-0.05, 0.05),
        diverging = TRUE
      )
    }
  }
}

# ---------------------------------------------------------------------
# 4. GLC vs CCI
# ---------------------------------------------------------------------
glc <- slopes[sapply(slopes, function(s) s$info$mask == "GLC")]
cci <- slopes[sapply(slopes, function(s) s$info$mask == "CCI")]

for (var in c("FPAR","LAI")){
  vglc <- glc[sapply(glc, function(s) s$info$var == var)]
  vcci <- cci[sapply(cci, function(s) s$info$var == var)]

  for (g in vglc){
    t <- g$info$tau
    ccimatch <- vcci[sapply(vcci, function(s) s$info$tau == t)]
    if (length(ccimatch) == 0) next

    diff_r <- g$r - ccimatch[[1]]$r

    outdir <- file.path(OUTP, "glc_vs_cci", var)
    dir.create(outdir, recursive = TRUE)

    outpng <- file.path(outdir, sprintf("%s_GLC_vs_CCI_tau_%s.png", var, t))

    plot_global_map(
      diff_r,
      sprintf("%s: GLC - CCI (τ=%s)\nTrend difference per decade", var, t),
      outpng,
      lim = c(-0.05, 0.05),
      diverging = TRUE
    )
  }
}

# ---------------------------------------------------------------------
# 5. Scatter plots
# ---------------------------------------------------------------------
df_scatter <- function(r1, r2){
  d <- cbind(values(r1), values(r2))
  d <- d[complete.cases(d), ]
  colnames(d) <- c("x","y")
  as.data.frame(d)
}

plot_scatter <- function(d, title, outfile){
  p <- ggplot(d, aes(x=x, y=y)) +
    geom_point(alpha=0.25, size=0.5, color="#3366AA") +
    geom_abline(intercept=0, slope=1, color="red") +
    coord_fixed() +
    labs(x="Variant 1", y="Variant 2", title=title) +
    theme_minimal(base_size=14)

  ggsave(outfile, p, width=8, height=6)
}

um <- unmasked[[1]]
masked <- slopes[sapply(slopes, function(s) s$info$mask != "UNMASKED")]

for (mk in masked){
  if (mk$info$var != um$info$var) next

  d <- df_scatter(um$r, mk$r)

  outdir <- file.path(OUTP, "scatter", mk$info$var)
  dir.create(outdir, recursive = TRUE)

  outpng <- file.path(
    outdir,
    sprintf("scatter_%s_masked_vs_unmasked_%s.png",
            mk$info$var, mk$info$tau)
  )

  plot_scatter(
    d,
    sprintf("%s masked (%s) vs unmasked", mk$info$var, mk$info$tau),
    outpng
  )
}

# ---------------------------------------------------------------------
# 6. Latitude-binned summaries
# ---------------------------------------------------------------------
lat_summary <- function(r){
  xy <- crds(r, df=TRUE)
  vals <- values(r)
  d <- cbind(xy, slope=vals)
  d <- d[complete.cases(d),]
  d$latbin <- cut(d$y, breaks=seq(-90,90,5), include.lowest=TRUE)
  d %>% group_by(latbin) %>% summarize(mean_slope = mean(slope))
}

for (s in slopes){
  info <- s$info
  df <- lat_summary(s$r)

  outdir <- file.path(OUTP, "lat_summary", info$var, info$mask)
  dir.create(outdir, recursive = TRUE)

  outpng <- file.path(
    outdir,
    sprintf("latitudinal_slope_%s_%s_%s.png",
            info$var, info$mask, info$tau)
  )

  p <- ggplot(df, aes(x=latbin, y=mean_slope)) +
    geom_bar(stat="identity", fill=hcl.colors(5, "batlow")[3]) +
    coord_flip() +
    labs(
      title = sprintf("%s %s %s: Latitudinal Mean Trend",
                      info$var, info$mask, info$tau),
      x="Latitude band",
      y="Slope per decade"
    ) +
    theme_minimal(base_size=13)

  ggsave(outpng, p, width=8, height=8)
}

# ---------------------------------------------------------------------
# 7. AOI time series (placeholder)
# ---------------------------------------------------------------------
aois <- list(
  amazon = list(lon_min=-80, lon_max=-45, lat_min=-20, lat_max=6),
  europe = list(lon_min=-25, lon_max=45, lat_min=33, lat_max=72)
)

extract_aoi_timeseries <- function(filepath, aoi){
  r <- rast(filepath)
  e <- ext(aoi$lon_min,aoi$lon_max,aoi$lat_min,aoi$lat_max)
  rs <- crop(r,e)
  ts <- colMeans(values(rs), na.rm=TRUE)
  ts
}

# (section can be filled later)

# ---------------------------------------------------------------------
# 8. Significance maps (MK)
# ---------------------------------------------------------------------
mk_significance <- function(monthly_nc){
  x <- rast(monthly_nc)
  vals <- values(x)
  pvals <- apply(vals, 1, function(ts){
    if (all(is.na(ts))) return(NA)
    res <- try(MannKendall(ts), silent=TRUE)
    if (inherits(res, "try-error")) return(NA)
    res$sl
  })
  r <- rast(x[[1]])
  values(r) <- pvals
  r
}

# Example:
# pval_r <- mk_significance("monthly_0p25.nc")

# ---------------------------------------------------------------------
# 9. Agreement maps
# ---------------------------------------------------------------------
agreement_map <- function(r1, r2){
  a <- sign(values(r1))
  b <- sign(values(r2))
  out <- ifelse(is.na(a)|is.na(b), NA,
                ifelse(a == b, 1, -1))
  rr <- rast(r1); values(rr) <- out; rr
}

# Example:
# ag <- agreement_map(cci_r, glc_r)

cat("All diagnostics computed in:", OUTP, "\n")
