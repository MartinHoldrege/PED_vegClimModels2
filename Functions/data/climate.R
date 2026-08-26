# //////////////////////////////////////////////////////////////////////////
# climate_functions.R
#
# Helper functions for calculating climate metrics from monthly raster stacks
# (tmin, tmax, prcp). Used to derive end-of-century climate "normals" from MACA
# projections so they can serve as inputs to the vegetation-climate models.
#
# Metric definitions mirror the Daymet training pipeline (06_GettingWeatherData.R),
# with the deliberate exceptions noted inline:
#   - durationFrostFreeDays uses a robust month -> day-of-year lookup instead of
#     the fragile string-date construction in the training code.
#   - weighted_annual_mean() defaults to the mathematically correct denominator
#     sum(weights)
#
# //////////////////////////////////////////////////////////////////////////

# Day-of-month weights used for monthly -> annual averaging. February is given
# 28.25 days, matching the training pipeline. Order is Jan..Dec.
.month_weights <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


#' Day-weighted monthly-to-annual mean
#'
#' Computes an annual mean from 12 monthly values, weighting each month by its
#' length. By default the denominator is `sum(weights)` (365.5), giving a
#' properly normalised weighted mean.
#'
#'
#' @param x12 Numeric vector of 12 monthly values (Jan..Dec). May contain NA.
#' @param weights Numeric vector of 12 month-length weights.
#' @param denom Denominator for the weighted sum. Defaults to `sum(weights)`.
#' @return A single numeric value (the day-weighted annual mean).
weighted_annual_mean <- function(x12,
                                 weights = .month_weights,
                                 denom = sum(weights)) {
  sum(x12 * weights) / denom
}


#' Precipitation seasonality (coefficient of variation)
#'
#' @param prcp12 Numeric vector of 12 monthly precipitation totals (Jan..Dec).
#' @return sd / mean of monthly precipitation. NA/NaN if mean precip is 0;
#'   the caller substitutes 2 for such pixels (matching training).
precip_seasonality <- function(prcp12) {
  stats::sd(prcp12) / mean(prcp12)
}


#' Correlation between monthly precipitation and monthly maximum temperature
#'
#' @param prcp12 Numeric vector of 12 monthly precipitation totals (Jan..Dec).
#' @param tmax12 Numeric vector of 12 monthly mean-tmax values (Jan..Dec).
#' @return Pearson correlation. NA if either input has zero variance; the
#'   caller substitutes -0.25 for such pixels (matching training).
precip_temp_corr <- function(prcp12, tmax12) {
  stats::cor(x = prcp12, y = tmax12)
}


#' Isothermality
#'
#' Mean monthly temperature range (tmax - tmin) divided by the annual range
#' (max tmax - min tmin), as a percentage.
#'
#' @param tmin12 Numeric vector of 12 monthly tmin values (Jan..Dec).
#' @param tmax12 Numeric vector of 12 monthly tmax values (Jan..Dec).
#' @return Isothermality (percent).
isothermality <- function(tmin12, tmax12) {
  ann_range <- max(tmax12) - min(tmin12)
  mean(tmax12 - tmin12) / ann_range * 100
}


#' First month with night-time temperatures above freezing
#'
#' @param tmin12 Numeric vector of 12 monthly tmin values (Jan..Dec).
#' @return Month number (1-12) of the first month with tmin > 0, or NA if no
#'   month is above freezing. Caller substitutes 8 for NA (matching training).
first_above_freezing <- function(tmin12) {
  which(tmin12 > 0)[1]
}


#' Last month with night-time temperatures above freezing
#'
#' @param tmin12 Numeric vector of 12 monthly tmin values (Jan..Dec).
#' @return Month number (1-12) of the last month with tmin > 0, or NA if no
#'   month is above freezing.
last_above_freezing <- function(tmin12) {
  above <- which(tmin12 > 0)
  if (length(above) > 0) max(above) else NA_real_
}


# Day-of-year of the first day of each month (non-leap year). Index by month.
.first_doy <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
# Day-of-year of the last day of each month (non-leap year). Index by month.
.last_doy  <- c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)


#' Duration of the frost-free period, in days
#'
#' Defined as (last day of the last above-freezing month) minus (first day of
#' the first above-freezing month). 
#'
#' @param above_month First above-freezing month (1-12) or NA.
#' @param last_above_month Last above-freezing month (1-12) or NA.
#' @return Frost-free duration in days, or NA if either month is NA. Caller
#'   substitutes 0 for NA (matching training intent).
frost_free_days <- function(above_month, last_above_month) {
  if (is.na(above_month) || is.na(last_above_month)) return(NA_real_)
  .last_doy[last_above_month] - .first_doy[above_month]
}

# saturation vapor pressure
svp <- function(x) {
  # constants for SVP calculation
  #calculate SVP according to Williams et al NatCC 2012 supplementary material -  units Pa (Pascals)
  # https://static-content.springer.com/esm/art%3A10.1038%2Fnclimate1693/MediaObjects/41558_2013_BFnclimate1693_MOESM272_ESM.pdf
  # this is Horner's method
  a0<-6.107799961
  a1<-0.4436518521
  a2<-0.01428945805
  a3<-0.0002650648471
  a4<-0.000003031240396
  a5<-0.00000002034080948
  a6<-0.00000000006136820929
  svp_hpa <- (a0+ x*(a1+ x *(a2+ x *(a3+ x *(a4	+ x *(a5	+ x *a6)))))) # eq S1
  svp_hpa*100 # convert hPa to Pa
}

#' Calculate all per-year annual climate metrics for one year
#'
#' Takes three 12-layer monthly stacks (tmin, tmax, prcp, all Jan..Dec, in
#' degrees C and mm) on a common grid and returns a multi-layer SpatRaster of
#' the annual metrics used downstream. Per-year NA substitutions are applied
#' here so that across-year aggregation sees no NAs:
#'   precip_Seasonality -> 2, PrecipTempCorr -> -0.25,
#'   aboveFreezing_month -> 8, durationFrostFreeDays -> 0.
#'
#' @param tmin12 SpatRaster, 12 layers of monthly tmin (Jan..Dec), deg C.
#' @param tmax12 SpatRaster, 12 layers of monthly tmax (Jan..Dec), deg C.
#' @param prcp12 SpatRaster, 12 layers of monthly precip (Jan..Dec), mm.
#' @param vp12 SpatRaster, 12 layers of monthly vapor pressure, Pa (optional)
#' @param denom Denominator passed to `weighted_annual_mean()` for the
#'   day-weighted temperature/VPD means. Defaults to the correct `sum(weights)`.
#' @return A SpatRaster with one layer per annual metric.
calc_annual_metrics <- function(tmin12, tmax12, prcp12,
                                vp12 = NULL,
                                denom = sum(.month_weights)) {

  # ---- monthly mean temperature and monthly VPD (12-layer stacks) ----
  tmean12 <- (tmax12 + tmin12) / 2
  # svp() is polynomial arithmetic that works directly on SpatRasters,
  # avoiding the per-pixel R function overhead of app().
  
  # if no vapor pressure data (e.g. for testing projections
  # for maca models where we just have tmin/tmax)
  if(is.null(vp12)) {
    vpd12 <- svp(tmax12) - svp(tmin12)
  } else {
    # for daymet data
    vpd12 <- svp(tmax12) - vp12
  }
  

  # ---- simple monthly reductions ----
  totalAnnPrecip   <- sum(prcp12)
  T_warmestMonth   <- max(tmax12)
  T_coldestMonth   <- min(tmin12)
  precip_wettest   <- max(prcp12)
  precip_driest    <- min(prcp12)

  # ---- day-weighted annual means (direct raster math, avoids app() overhead) ----
  weighted_raster_mean <- function(r12, weights = .month_weights, d = denom) {
    out <- r12[[1]] * weights[1]
    for (i in 2:12) out <- out + r12[[i]] * weights[i]
    out / d
  }
  Tmin_annAvg  <- weighted_raster_mean(tmin12)
  Tmax_annAvg  <- weighted_raster_mean(tmax12)
  tmean_annAvg <- weighted_raster_mean(tmean12)
  annVPD_mean  <- weighted_raster_mean(vpd12)

  # ---- VPD extremes ----
  annVPD_max <- max(vpd12)
  annVPD_min <- min(vpd12)

  r_not_na <- !is.na(Tmin_annAvg)
  
  # ---- precipitation seasonality (CV), NA -> 2 ----
  precip_Seasonality <- terra::app(prcp12, fun = precip_seasonality)
  precip_Seasonality <- terra::ifel(is.na(precip_Seasonality) & r_not_na, 2, precip_Seasonality)

  # ---- precip-temp correlation (NA -> -0.25) ----
  # app over a 24-layer stack: first 12 = prcp, last 12 = tmax.
  PrecipTempCorr <- terra::app(c(prcp12, tmax12),
                               fun = function(x) precip_temp_corr(x[1:12], x[13:24]))
  PrecipTempCorr <- terra::ifel(is.na(PrecipTempCorr) & r_not_na, -0.25, PrecipTempCorr)

  # ---- isothermality (direct raster math) ----
  isotherm <- mean(tmax12 - tmin12) / (max(tmax12) - min(tmin12)) * 100

  # ---- thaw timing and frost-free duration ----
  aboveFreezing_month     <- terra::app(tmin12, fun = first_above_freezing)
  lastAboveFreezing_month <- terra::app(tmin12, fun = last_above_freezing)

  durationFrostFreeDays <- terra::lapp(
    c(aboveFreezing_month, lastAboveFreezing_month),
    fun = function(a, l) mapply(frost_free_days, a, l)
  )
  durationFrostFreeDays <- terra::ifel(is.na(durationFrostFreeDays) & r_not_na, 0, durationFrostFreeDays)

  # aboveFreezing_month NA -> 8 (after it has been used for frost-free days)
  aboveFreezing_month <- terra::ifel(is.na(aboveFreezing_month) & r_not_na, 8, aboveFreezing_month)

  # ---- monthly water deficit and wet degree days ----
  # awd = tmean*2 - prcp (Thornthwaite-style approximation, ported as-is).
  awd12 <- tmean12 * 2 - prcp12
  annWaterDeficit <- terra::app(awd12, fun = function(x) sum(x[x > 0]))

  # wet degree days: months where tmean*2 < prcp contribute tmean*30, else 0.
  awdd12 <- terra::ifel(tmean12 * 2 < prcp12, tmean12 * 30, 0)
  annWetDegDays <- terra::app(awdd12, fun = function(x) sum(x[x > 0]))

  out <- c(totalAnnPrecip, T_warmestMonth, T_coldestMonth,
           Tmin_annAvg, Tmax_annAvg, tmean_annAvg,
           precip_wettest, precip_driest, precip_Seasonality,
           PrecipTempCorr, aboveFreezing_month, isotherm,
           annWaterDeficit, annWetDegDays,
           annVPD_mean, annVPD_max, annVPD_min,
           durationFrostFreeDays)

  names(out) <- c("totalAnnPrecip", "T_warmestMonth", "T_coldestMonth",
                  "tmin_annAvg", "tmax_annAvg", "tmean",
                  "precip_wettestMonth", "precip_driestMonth", "precip_Seasonality",
                  "PrecipTempCorr", "aboveFreezing_month", "isothermality",
                  "annWaterDeficit", "annWetDegDays",
                  "annVPD_mean", "annVPD_max", "annVPD_min",
                  "durationFrostFreeDays")
  out
}


# Most metrics: mean across years. A few use percentiles, matching training.

# Across-year reductions --------------------------------------------------
# Two implementations of the same three reductions: one for multi-layer
# SpatRasters (year = layer), one for numeric vectors (year = element).

.reduction_raster <- list(
  mean = function(r) terra::app(r, fun = mean, na.rm = TRUE),
  q95  = function(r) terra::app(r, fun = function(x) stats::quantile(x, 0.95, na.rm = TRUE)),
  q05  = function(r) terra::app(r, fun = function(x) stats::quantile(x, 0.05, na.rm = TRUE))
)


.reduction_vector <- list(
  mean = function(x) mean(x, na.rm = TRUE),
  q95  = function(x) unname(stats::quantile(x, 0.95, na.rm = TRUE)),
  q05  = function(x) unname(stats::quantile(x, 0.05, na.rm = TRUE))
)


# Map: per-year metric name -> list(output _CLIM name, reduction name,
# anomaly type). Anomaly type controls how the trailing-window anomaly is
# expressed: "abs" for quantities on an interpretable absolute scale
# (temperatures, VPD, correlations, day counts), "pct" for quantities where a
# proportional change is more meaningful (precipitation, deficits).
climate_reductions <- list(
  list("tmin_annAvg",            "tmin_meanAnnAvg_CLIM",                   "mean", "abs"),
  list("tmax_annAvg",            "tmax_meanAnnAvg_CLIM",                   "mean", "abs"),
  list("tmean",                  "tmean_meanAnnAvg_CLIM",                  "mean", "abs"),
  list("totalAnnPrecip",         "prcp_meanAnnTotal_CLIM",                 "mean", "pct"),
  list("T_warmestMonth",         "T_warmestMonth_meanAnnAvg_CLIM",         "mean", "abs"),
  list("T_coldestMonth",         "T_coldestMonth_meanAnnAvg_CLIM",         "mean", "abs"),
  list("precip_wettestMonth",    "precip_wettestMonth_meanAnnAvg_CLIM",    "mean", "pct"),
  list("precip_driestMonth",     "precip_driestMonth_meanAnnAvg_CLIM",     "mean", "pct"),
  list("precip_Seasonality",     "precip_Seasonality_meanAnnAvg_CLIM",     "mean", "pct"),
  list("PrecipTempCorr",         "PrecipTempCorr_meanAnnAvg_CLIM",         "mean", "abs"),
  list("aboveFreezing_month",    "aboveFreezing_month_meanAnnAvg_CLIM",    "mean", "abs"),
  list("isothermality",          "isothermality_meanAnnAvg_CLIM",          "mean", "abs"),
  list("annWaterDeficit",        "annWaterDeficit_meanAnnAvg_CLIM",        "mean", "pct"),
  list("annWetDegDays",          "annWetDegDays_meanAnnAvg_CLIM",          "mean", "pct"),
  list("annVPD_mean",            "annVPD_mean_meanAnnAvg_CLIM",            "mean", "abs"),
  list("annVPD_max",             "annVPD_max_meanAnnAvg_CLIM",             "mean", "abs"),
  list("annVPD_min",             "annVPD_min_meanAnnAvg_CLIM",             "mean", "abs"),
  list("annVPD_max",             "annVPD_max_95percentile_CLIM",           "q95",  "abs"),
  list("annWaterDeficit",        "annWaterDeficit_95percentile_CLIM",      "q95",  "pct"),
  list("annWetDegDays",          "annWetDegDays_5percentile_CLIM",         "q05",  "pct"),
  list("durationFrostFreeDays",  "durationFrostFreeDays_5percentile_CLIM", "q05",  "abs"),
  list("durationFrostFreeDays",  "durationFrostFreeDays_meanAnnAvg_CLIM",  "mean", "abs")
)

#' Anomaly of a trailing window relative to the long-term normal
#'
#' The anomaly is recent minus normal, so a recent period that is warmer or
#' wetter than normal gives a positive anomaly. Note this is the opposite sign
#' to some older code
#'
#' Percent anomalies are undefined when the normal is zero. Where both the
#' normal and the recent window are zero (e.g. driest-month precipitation in
#' the desert Southwest) the anomaly is set to 0: the recent period is not
#' different from normal. Where the normal is zero but the recent window is
#' not, the anomaly is left NA rather than being forced to a finite value.
#'
#' @param dat Data frame containing the `_CLIM` and short-window columns.
#' @param short_suffix Suffix of the short-window columns (e.g. "_3yr").
#' @return Data frame of anomaly columns, named `<metric><short_suffix>Anom`.
calc_anomalies <- function(dat, short_suffix = "_3yr") {
  clim_names  <- purrr::map_chr(climate_reductions, 2)
  anom_types  <- purrr::map_chr(climate_reductions, 4)
  short_names <- stringr::str_replace(clim_names, "_CLIM$", short_suffix)
  anom_names  <- stringr::str_replace(clim_names, "_CLIM$",
                                      paste0(short_suffix, "Anom"))
  
  missing <- setdiff(c(clim_names, short_names), names(dat))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  out <- purrr::pmap(list(clim_names, short_names, anom_types),
                     function(cl, sh, type) {
                       normal <- dat[[cl]]
                       recent <- dat[[sh]]
                       if (type == "abs") return(recent - normal)
                       dplyr::case_when(
                         normal == 0 & recent == 0 ~ 0,
                         normal == 0               ~ NA_real_,
                         TRUE                      ~ (recent - normal) / normal
                       )
                     })
  purrr::set_names(as.data.frame(out), anom_names)
}


# Row-wise helpers --------------------------------------------------------
# Vectorised equivalents of the scalar helpers, for n x 12 matrices. Kept
# separate (rather than mapping the scalar versions over rows) because these run over
# millions of rows; each is mathematically identical to its scalar counterpart.

#' Row-wise maximum / minimum of a matrix
#' @param m Numeric matrix.
#' @return Numeric vector of length nrow(m).
.row_max <- function(m) purrr::reduce(as.data.frame(m), pmax)
.row_min <- function(m) purrr::reduce(as.data.frame(m), pmin)

#' Row-wise sample standard deviation
#' @param m Numeric matrix.
#' @return Numeric vector of length nrow(m).
.row_sd <- function(m) {
  centred <- m - rowMeans(m)
  sqrt(rowSums(centred^2) / (ncol(m) - 1))
}

#' Row-wise Pearson correlation between two matrices
#'
#' Correlates row i of `x` against row i of `y`, using stats::cor() on each
#' pair of rows. Returns NA where either row has zero variance (cor() warns in
#' that case; the warning is suppressed because the NA is expected and is
#' substituted by the caller).
#'
#' @param x,y Numeric matrices of identical dimension.
#' @return Numeric vector of length nrow(x).
.row_cor <- function(x, y) {
  stopifnot(identical(dim(x), dim(y)))
  suppressWarnings(
    purrr::map2_dbl(asplit(x, 1), asplit(y, 1),
                    function(xi, yi) stats::cor(x = xi, y = yi))
  )
}


#' Calculate all per-year annual climate metrics for point data
#'
#' Data-frame analogue of `calc_annual_metrics()`. Takes 12-column matrices of
#' monthly values (Jan..Dec, one row per site-year) and returns a data frame of
#' the same annual metrics, with the same NA substitutions.
#'
#' @param tmin12,tmax12,prcp12 Numeric matrices, n x 12: monthly tmin (deg C),
#'   tmax (deg C) and precipitation total (mm).
#' @param vp12 Numeric matrix, n x 12, of monthly vapour pressure (Pa), or NULL
#'   to approximate the dew point with tmin (matching the MACA branch).
#' @param denom Denominator for the day-weighted annual means.
#' @return Data frame with one row per input row and one column per metric,
#'   using the same column names as `calc_annual_metrics()`.
calc_annual_metrics_df <- function(tmin12, tmax12, prcp12, vp12 = NULL,
                                   denom = sum(.month_weights)) {
  stopifnot(is.matrix(tmin12), ncol(tmin12) == 12,
            identical(dim(tmin12), dim(tmax12)),
            identical(dim(tmin12), dim(prcp12)))
  if (!is.null(vp12)) stopifnot(identical(dim(tmin12), dim(vp12)))
  
  tmean12 <- (tmax12 + tmin12) / 2
  vpd12 <- if (is.null(vp12)) svp(tmax12) - svp(tmin12) else svp(tmax12) - vp12
  
  # Day-weighted annual mean: matrix multiply by the month weights.
  wmean <- function(m) as.vector(m %*% .month_weights) / denom
  
  tmin_annAvg <- wmean(tmin12)
  r_not_na    <- !is.na(tmin_annAvg)
  
  # Precip seasonality (CV), NA -> 2.
  precip_Seasonality <- .row_sd(prcp12) / rowMeans(prcp12)
  precip_Seasonality[is.na(precip_Seasonality) & r_not_na] <- 2
  
  # Precip-tmax correlation, NA -> -0.25.
  PrecipTempCorr <- .row_cor(prcp12, tmax12)
  PrecipTempCorr[is.na(PrecipTempCorr) & r_not_na] <- -0.25
  
  # Thaw timing. max.col() gives the first/last TRUE; rows with no month above
  # freezing get NA.
  above    <- tmin12 > 0
  has_any  <- rowSums(above) > 0
  first_ab <- max.col(above, ties.method = "first")
  last_ab  <- 13 - max.col(above[, 12:1, drop = FALSE], ties.method = "first")
  first_ab[!has_any] <- NA_integer_
  last_ab[!has_any]  <- NA_integer_
  
  # Frost-free duration, NA -> 0. aboveFreezing_month NA -> 8, applied after.
  durationFrostFreeDays <- .last_doy[last_ab] - .first_doy[first_ab]
  durationFrostFreeDays[is.na(durationFrostFreeDays) & r_not_na] <- 0
  aboveFreezing_month <- first_ab
  aboveFreezing_month[is.na(aboveFreezing_month) & r_not_na] <- 8
  
  # Monthly water deficit and wet degree days; sum positive months only.
  awd12  <- tmean12 * 2 - prcp12
  awdd12 <- ifelse(tmean12 * 2 < prcp12, tmean12 * 30, 0)
  
  data.frame(
    totalAnnPrecip        = rowSums(prcp12),
    T_warmestMonth        = .row_max(tmax12),
    T_coldestMonth        = .row_min(tmin12),
    tmin_annAvg           = tmin_annAvg,
    tmax_annAvg           = wmean(tmax12),
    tmean                 = wmean(tmean12),
    precip_wettestMonth   = .row_max(prcp12),
    precip_driestMonth    = .row_min(prcp12),
    precip_Seasonality    = precip_Seasonality,
    PrecipTempCorr        = PrecipTempCorr,
    aboveFreezing_month   = aboveFreezing_month,
    isothermality         = rowMeans(tmax12 - tmin12) /
      (.row_max(tmax12) - .row_min(tmin12)) * 100,
    annWaterDeficit       = rowSums(pmax(awd12, 0)),
    annWetDegDays         = rowSums(pmax(awdd12, 0)),
    annVPD_mean           = wmean(vpd12),
    annVPD_max            = .row_max(vpd12),
    annVPD_min            = .row_min(vpd12),
    durationFrostFreeDays = durationFrostFreeDays
  )
}


#' Reduce per-year point metrics over a trailing window
#'
#' For each requested site-year, applies the `climate_reductions` across the
#' `n_years` years *preceding* that year (the observation year itself is
#' excluded). Windows are truncated where the record does not extend far enough
#' back; `min_years` sets how short a window is still acceptable.
#'
#' @param annual Data frame of per-year metrics, with columns `cell`, `year`,
#'   and one column per metric named in `climate_reductions`.
#' @param targets Data frame of the site-years to summarise, with columns
#'   `cell` and `year`.
#' @param n_years Length of the trailing window (e.g. 30 or 3).
#' @param out_suffix Suffix for the output columns, replacing the `_CLIM` in
#'   `climate_reductions` (e.g. "_CLIM" or "_3yr").
#' @param min_years Minimum number of years required; site-years with fewer
#'   available get NA. Defaults to `n_years` (no truncation allowed).
#' @return Data frame with `cell`, `year`, `n_years_used`, and one column per
#'   reduction.
roll_point_normals <- function(annual, targets, n_years, out_suffix,
                               min_years = n_years) {
  metrics   <- purrr::map_chr(climate_reductions, 1)
  out_names <- purrr::map_chr(climate_reductions, 2)
  out_names <- stringr::str_replace(out_names, "_CLIM$", out_suffix)
  red_names <- purrr::map_chr(climate_reductions, 3)
  
  target_years <- sort(unique(targets$year))
  
  purrr::map_dfr(target_years, function(yr) {
    cells_yr <- targets$cell[targets$year == yr]
    window   <- annual[annual$year >= yr - n_years & annual$year <= yr - 1 &
                         annual$cell %in% cells_yr, , drop = FALSE]
    if (nrow(window) == 0) return(NULL)
    
    split_by_cell <- split(window, window$cell)
    
    purrr::map_dfr(split_by_cell, function(d) {
      vals <- purrr::map2(metrics, red_names, function(m, red) {
        .reduction_vector[[red]](d[[m]])
      })
      out <- purrr::set_names(as.data.frame(vals), out_names)
      out$cell <- d$cell[1]
      out$year <- yr
      out$n_years_used <- nrow(d)
      out
    })
  }) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(out_names),
                                ~ ifelse(n_years_used >= min_years, .x, NA_real_)))
}
