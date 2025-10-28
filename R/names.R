# ===============================
# File: ref/names.R
# ===============================
name_stack_cat <- function(dataset = c("cci", "landsat")) {
  sprintf("%s_cat_yearstack_0p05.tif", match.arg(dataset))
}
name_frac <- function(kind, year) {
  sprintf("cci_frac_%s_%d_0p05.tif", kind, year)
}
name_cat_year <- function(dataset, year) {
  sprintf("%s_cat_%d_0p05.tif", dataset, year)
}
