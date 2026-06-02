# Functions for processing eVIIRS phenology metrics.
#
# These read the per-metric, per-year 375 m eVIIRS GeoTIFFs, decode the
# producer NoData/water codes and scaling, reproject each year to the Daymet
# 1 km LCC grid, and collapse across years. Day-of-year (timing) metrics are
# handled with circular encoding during spatial aggregation to avoid the
# Dec 31 / Jan 1 seam, which is a real concern for a CONUS-wide product
# (winter-active vegetation in the South).
#
# Assumes `daymet_grid` and `crs_daymet` are available (via source('Functions/init.R')).

library(terra)


# Metric metadata -------------------------------------------------------------

#' Lookup table of eVIIRS phenology metric properties.
#'
#' One row per metric we model (the 6 primitives + TIN). Encodes the valid
#' range, the codes to set to NA, the scaling needed to recover real units,
#' and whether the metric is a day-of-year ("timing") quantity that needs
#' circular handling.
#'
#' @return A tibble with columns: metric, type, valid_min, valid_max,
#'   na_codes (list), scale_fun (list of functions), is_circular, period.
#' @export
phenology_metric_info <- function() {
  # Scaling: NDVI metrics use (value - 100) * 0.01 to recover unscaled NDVI.
  # Timing metrics (SOST, EOST, MAXT) are already in day-of-year units.
  # TIN and AMP are unitless indices (no rescale needed for modelling).
  ndvi_scale <- function(x) (x - 100) * 0.01
  identity_scale <- function(x) x

  tribble(
    ~metric, ~type,     ~valid_min, ~valid_max, ~na_codes,             ~scale_fun,     ~is_circular, ~period,
    "SOST",  "timing",  -150,        365,        list(c(-1000, 1000)),  list(identity_scale), TRUE,  365,
    "EOST",  "timing",     1,        450,        list(c(-1000, 1000)),  list(identity_scale), TRUE,  365,
    "MAXT",  "timing",     1,        365,        list(c(-1000, 1000)),  list(identity_scale), FALSE, 365,
    "SOSN",  "ndvi",     101,        200,        list(c(100, 255)),     list(ndvi_scale),     FALSE, NA,
    "EOSN",  "ndvi",     101,        200,        list(c(100, 255)),     list(ndvi_scale),     FALSE, NA,
    "MAXN",  "ndvi",     101,        200,        list(c(100, 255)),     list(ndvi_scale),     FALSE, NA,
    "TIN",   "index",      1,        200,        list(c(0, 255)),       list(identity_scale), FALSE, NA
  )
}


# File discovery --------------------------------------------------------------

#' Build a tibble of phenology GeoTIFF paths by metric and year.
#'
#' Globs the raw phenology directory rather than constructing paths, because the
#' source folder names contain an inconsistent space ("eVIIRS _Phenology").
#'
#' @param dir Directory containing the per-year subfolders of .tif files.
#' @param metrics Character vector of metric codes to keep (default: the 7 we model).
#' @return A tibble with columns metric, year, path.
#' @export
find_phenology_files <- function(dir,
                                 metrics = c("SOST", "EOST", "MAXT",
                                             "SOSN", "EOSN", "MAXN", "TIN")) {
  tifs <- list.files(dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
  stopifnot("No .tif files found" = length(tifs) > 0)

  files <- tibble(path = tifs) |>
    mutate(
      fname  = basename(path),
      metric = str_extract(fname, "_CONUS_([A-Z]+)_\\d{4}\\.tif$", group = 1),
      year   = as.integer(str_extract(fname, "(\\d{4})\\.tif$", group = 1))
    ) |>
    filter(metric %in% metrics) |>
    select(metric, year, path) |>
    arrange(metric, year)

  # Sanity: each requested metric should appear for each year present.
  n_years <- n_distinct(files$year)
  counts <- count(files, metric)
  if (any(counts$n != n_years)) {
    warning("Some metrics are missing years:\n",
            paste(capture.output(print(counts)), collapse = "\n"))
  }
  files
}


# Decoding --------------------------------------------------------------------

#' Decode a raw phenology raster: set NoData/water codes to NA and clamp.
#'
#' Sets the producer codes to NA, then sets anything outside the documented
#' valid range to NA as well (defensive: catches unexpected fill values). Does
#' NOT rescale here -- rescaling happens after aggregation so the NA handling
#' and circular logic operate on the native integer encoding.
#'
#' @param r A SpatRaster (single layer) of raw metric values.
#' @param info One row of `phenology_metric_info()` for this metric.
#' @return A SpatRaster with invalid cells set to NA (still in raw units).
#' @export
decode_phenology <- function(r, info) {
  na_codes <- info$na_codes[[1]][[1]]
  r <- terra::classify(r, cbind(na_codes, NA))
  # Clamp to valid range -> NA outside it (values = FALSE keeps in-range only).
  r <- terra::clamp(r, lower = info$valid_min, upper = info$valid_max, values = FALSE)
  r
}


# Spatial aggregation (375 m -> Daymet 1 km) ----------------------------------

#' Project one decoded year of a metric to the Daymet 1 km grid.
#'
#' Continuous metrics are projected with method = "average". Timing metrics
#' flagged circular (SOST, EOST) are encoded to sin/cos on their period,
#' each component projected with "average", then decoded back to day-of-year.
#' This makes the spatial averaging seam-safe (e.g. DOY 360 and DOY 5 average
#' to ~+2/-3, not ~180). MAXT is bounded 1-365 and projected linearly.
#'
#' Also returns a "valid fraction" layer: the proportion of contributing 375 m
#' cells that were non-NA, obtained by projecting a 0/1 validity mask with
#' "average".
#'
#' @param r Decoded SpatRaster (single layer, raw units, NAs set).
#' @param info One row of `phenology_metric_info()`.
#' @param grid Target grid (SpatRaster template), e.g. `daymet_grid`.
#' @return A SpatRaster with 2 layers: <metric> (on grid) and valid_frac.
#' @export
project_phenology_year <- function(r, info, grid) {
  # Validity mask -> fraction of valid fine pixels per coarse cell.
  valid <- terra::project(terra::ifel(is.na(r), 0, 1), grid, method = "average")
  names(valid) <- "valid_frac"

  if (isTRUE(info$is_circular)) {
    period <- info$period
    ang <- 2 * pi * r / period
    sin_p <- terra::project(sin(ang), grid, method = "average")
    cos_p <- terra::project(cos(ang), grid, method = "average")
    # Circular mean direction -> back to [0, period). Note: for SOST this maps
    # the negative (prior-year) starts onto the same circle as ~360, which is
    # exactly the intended wrap. Downstream, interpret values near `period` as
    # late-December starts.
    ang_mean <- terra::atan2(sin_p, cos_p)
    out <- (ang_mean / (2 * pi)) * period
    out <- terra::ifel(out < 0, out + period, out)
    names(out) <- info$metric
  } else {
    out <- terra::project(r, grid, method = "average")
    names(out) <- info$metric
  }
  c(out, valid)
}


# Cross-year collapse ---------------------------------------------------------

#' Collapse a multi-year stack of one metric to a single layer.
#'
#' Computes BOTH a linear mean and a linear median across years, plus their
#' absolute difference (a diagnostic for seam/outlier contamination). For
#' circular timing metrics, also computes a circular mean across years, so the
#' linear-vs-circular gap can be inspected. The caller decides which to keep.
#'
#' @param stk SpatRaster of the metric across years (one layer per year, raw units).
#' @param info One row of `phenology_metric_info()`.
#' @return A SpatRaster with layers: mean, median, mean_median_absdiff, and
#'   (if circular) circ_mean and circ_linear_absdiff.
#' @export
collapse_years <- function(stk, info) {
  m_mean   <- terra::app(stk, mean,   na.rm = TRUE)
  m_median <- terra::app(stk, median, na.rm = TRUE)
  names(m_mean) <- "mean"; names(m_median) <- "median"
  absdiff <- abs(m_mean - m_median); names(absdiff) <- "mean_median_absdiff"
  out <- c(m_mean, m_median, absdiff)

  if (isTRUE(info$is_circular)) {
    period <- info$period
    ang <- 2 * pi * stk / period
    s <- terra::app(sin(ang), mean, na.rm = TRUE)
    cc <- terra::app(cos(ang), mean, na.rm = TRUE)
    circ <- (terra::atan2(s, cc) / (2 * pi)) * period
    circ <- terra::ifel(circ < 0, circ + period, circ)
    names(circ) <- "circ_mean"
    cl_diff <- abs(m_mean - circ); names(cl_diff) <- "circ_linear_absdiff"
    out <- c(out, circ, cl_diff)
  }
  out
}


# Rescaling -------------------------------------------------------------------

#' Apply the metric's scaling to recover real units (e.g. NDVI from byte range).
#'
#' @param r SpatRaster in raw units.
#' @param info One row of `phenology_metric_info()`.
#' @return SpatRaster in real units.
#' @export
rescale_phenology <- function(r, info) {
  f <- info$scale_fun[[1]][[1]]
  f(r)
}
