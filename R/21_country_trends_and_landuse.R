## =============================================================================
## 21_country_trends_and_landuse.R
## Country-level LAI / fAPAR trends and land-use fractions (0.25°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(here)
  library(readr)
  library(rnaturalearth)
  library(exactextractr)
  library(lmtest)
  library(sandwich)
})
ROOT <- here::here()

DIR_UNM <- file.path(ROOT, "analysis", "unmasked", "0p25")

# Cropland / urban fraction rasters at 0.25°
# (0–1 fractions; binary masks also fine)
F_CROP  <- file.path(ROOT, "data/landuse/cropland_frac_0p25.tif")
F_URBAN <- file.path(ROOT, "data/landuse/urban_frac_0p25.tif")

OUTDIR <- file.path(ROOT, "analysis", "country_trends")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")

trend_test_hac <- function(df, y = "value", x = "year") {
  fit <- lm(reformulate(x, y), data = df)
  ct  <- lmtest::coeftest(fit, vcov. = sandwich::NeweyWest(fit))
  tibble(
    slope = unname(coef(fit)[2]),
    p     = ct[2, 4],
    r2    = summary(fit)$r.squared
  )
}
countries <- rnaturalearth::ne_countries(
  scale = 110,
  returnclass = "sf"
) |>
  st_transform("EPSG:4326") |>
  select(iso_a3, admin)

crop_r  <- rast(F_CROP)
urban_r <- rast(F_URBAN)

df_crop <- exactextractr::exact_extract(crop_r, countries, "mean") |>
  tibble(frac_cropland = .)

df_urban <- exactextractr::exact_extract(urban_r, countries, "mean") |>
  tibble(frac_urban = .)

df_landuse <- bind_cols(
  countries |> st_drop_geometry(),
  df_crop,
  df_urban
)

results <- list()

for (var in VARS) {
  for (met in METRICS) {

    message("Processing: ", var, " / ", met)

    f_nc <- file.path(DIR_UNM, sprintf("%s_georef_%s_0p25.nc", var, met))
    if (!file.exists(f_nc)) next

    r <- rast(f_nc)
    years <- 1982:(1982 + nlyr(r) - 1L)

    # ---- global reference trend ---------------------------------------------
    glob <- terra::global(r, "mean", na.rm = TRUE) |>
      as_tibble() |>
      rename(value = 1) |>
      mutate(year = years)

    base <- mean(glob$value[glob$year >= 1982 & glob$year <= 2000], na.rm = TRUE)
    if (!is.finite(base) || abs(base) < 1e-8) next

    glob <- glob |> mutate(value = 100 * value / base)

    M <- exactextractr::exact_extract(r, countries, "mean")
    colnames(M) <- years

    df_cty <- as.data.frame(M) |>
      bind_cols(countries |> st_drop_geometry())

    df_long <- df_cty |>
      pivot_longer(
        cols = as.character(years),
        names_to = "year",
        values_to = "value"
      ) |>
      mutate(year = as.integer(year))

    df_tr <- df_long |>
      group_by(iso_a3, admin) |>
      group_modify(~ {
        if (sum(is.finite(.x$value)) < 15) {
          return(tibble(slope = NA_real_, p = NA_real_, r2 = NA_real_))
        }
        base <- mean(.x$value[.x$year >= 1982 & .x$year <= 2000], na.rm = TRUE)
        if (!is.finite(base) || abs(base) < 1e-8) {
          return(tibble(slope = NA_real_, p = NA_real_, r2 = NA_real_))
        }
        .x <- mutate(.x, value = 100 * value / base)
        trend_test_hac(.x)
      }) |>
      ungroup()

    results[[length(results) + 1]] <- df_tr
  }
}

df_trends <- bind_rows(results)

df_out <- df_trends |>
  left_join(df_landuse, by = c("iso_a3", "admin"))

write_csv(
  df_out,
  file.path(OUTDIR, "country_trends_lai_fpar.csv")
)

message("Saved country trends to:\n  ", OUTDIR)
