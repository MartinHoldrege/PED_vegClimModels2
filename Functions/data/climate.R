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
#     sum(weights); the training code used a fixed denominator of 12 (see note
#     in that function).
#
# Per-pixel cross-month reductions are implemented with terra::app() applied to
# 12-layer stacks, so each function receives the 12 monthly values for a pixel as
# a plain numeric vector. This is slower than fully vectorised raster math but
# keeps the metric logic transparent and easy to check against the source.
# //////////////////////////////////////////////////////////////////////////

# Day-of-month weights used for monthly -> annual averaging. February is given
# 28.5 days, matching the training pipeline. Order is Jan..Dec.
.month_weights <- c(31, 28.5, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


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
  #calculate SVP according to Williams et al NatCC 2012 supplementary material -  units haPa
  # https://static-content.springer.com/esm/art%3A10.1038%2Fnclimate1693/MediaObjects/41558_2013_BFnclimate1693_MOESM272_ESM.pdf
  a0<-6.107799961
  a1<-0.4436518521
  a2<-0.01428945805
  a3<-0.0002650648471
  a4<-0.000003031240396
  a5<-0.00000002034080948
  a6<-0.00000000006136820929
  svp_hapa <- (a0+ x*(a1+ x *(a2+ x *(a3+ x *(a4	+ x *(a5	+ x *a6)))))) # eq S1
  svp_hapa
}

vpd <- function(tmean, tmin) {
  t_dewpoint <- tmin # approximation
  svp_mean <- svp(tmean)
  svp_dew_approx <- svp(t_dewpoint) # approximate actual vapor pressure
  vpd <- (svp_mean - svp_dew_approx)
  vpd
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
#' @param denom Denominator passed to `weighted_annual_mean()` for the
#'   day-weighted temperature/VPD means. Defaults to the correct `sum(weights)`.
#' @return A SpatRaster with one layer per annual metric.
calc_annual_metrics <- function(tmin12, tmax12, prcp12,
                                denom = sum(.month_weights)) {

  # ---- monthly mean temperature and monthly VPD (12-layer stacks) ----
  tmean12 <- (tmax12 + tmin12) / 2
  # svp() is polynomial arithmetic that works directly on SpatRasters,
  # avoiding the per-pixel R function overhead of app().
  vpd12 <- svp(tmean12) - svp(tmin12)

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

  # ---- precipitation seasonality (CV), NA -> 2 ----
  precip_Seasonality <- terra::app(prcp12, fun = precip_seasonality)
  precip_Seasonality <- terra::ifel(is.na(precip_Seasonality), 2, precip_Seasonality)

  # ---- precip-temp correlation (NA -> -0.25) ----
  # app over a 24-layer stack: first 12 = prcp, last 12 = tmax.
  PrecipTempCorr <- terra::app(c(prcp12, tmax12),
                               fun = function(x) precip_temp_corr(x[1:12], x[13:24]))
  PrecipTempCorr <- terra::ifel(is.na(PrecipTempCorr), -0.25, PrecipTempCorr)

  # ---- isothermality (direct raster math) ----
  isotherm <- mean(tmax12 - tmin12) / (max(tmax12) - min(tmin12)) * 100

  # ---- thaw timing and frost-free duration ----
  aboveFreezing_month     <- terra::app(tmin12, fun = first_above_freezing)
  lastAboveFreezing_month <- terra::app(tmin12, fun = last_above_freezing)

  durationFrostFreeDays <- terra::lapp(
    c(aboveFreezing_month, lastAboveFreezing_month),
    fun = function(a, l) mapply(frost_free_days, a, l)
  )
  durationFrostFreeDays <- terra::ifel(is.na(durationFrostFreeDays), 0, durationFrostFreeDays)

  # aboveFreezing_month NA -> 8 (after it has been used for frost-free days)
  aboveFreezing_month <- terra::ifel(is.na(aboveFreezing_month), 8, aboveFreezing_month)

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
