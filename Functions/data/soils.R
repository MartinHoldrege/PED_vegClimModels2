
#' Average over new intervals from point predictions
#' [function written by Daniel Schlaepfer]
#'
#' Weighted average of point predictions within new intervals using
#' the trapezoidal rule (see Hengl et al. 2017 PLoS ONE).
#'
#' Point values at new interval ends are estimated with linear interpolation
#' if not available. Extrapolation outside available point range is determined
#' by `extrapolation`.
#'
#' @param d Soil depth points
#' @param fd Soil properties at `d`
#' @param dnew New soil depth intervals (including uppermost and deepest depths)
#' @param extrapolation Method for extrapolation outside available point range,
#' see details.
#'
#' @section Details:
#' `extrapolation` method `"closest"` uses the value at the closest data
#' extreme;
#' `extrapolation` method `"naturalSpline"` builds a natural cubic spline to
#' extrapolate values outside available point range.
trapezoidSoilLayers <- function(
    d,
    fd,
    dnew,
    extraplation = c("closest", "naturalSpline")
) {
  # 1 / (b - a) * 1 / 2 * sum(k = 1; N - 1) {(x[k + 1] - x[k]) * (f(x[k + 1]) + f(x[k]))}
  extraplation <- match.arg(extraplation)
  res <- rep(NA_real_, length(dnew) - 1L)
  
  if (all(is.na(dnew)) || all(is.na(fd))) return(res)
  
  dd <- sort(unique(c(d, dnew)), na.last = NA) # remove NAs
  
  fdd <- if (sum(!is.na(fd)) > 1L) {
    # Need at least 2 values to interpolate
    stats::approx(
      x = d,
      y = fd,
      xout = dd,
      rule = switch(extraplation, closest = 2L, naturalSpline = 1L)
    )[["y"]]
  } else {
    tmp <- rep(NA_real_, length(dd))
    ids <- match(dd, d, nomatch = 0L)
    tmp[ids > 0L] <- fd[ids]
    tmp
  }
  
  if (extraplation == "naturalSpline") {
    isna <- which(is.na(fdd))
    if (length(isna)) {
      fdd[isna] <- stats::spline(
        x = d, y = fd, xout = dd[isna], method = "natural"
      )[["y"]]
    }
  }
  
  for (kl in seq_along(res)) {
    ids <- which(dnew[[kl]] <= dd & dd <= dnew[[kl + 1L]])
    nids <- length(ids)
    if (nids == 0L) next
    if (nids == 1L) {
      res[[kl]] <- fdd[[ids]]
    } else {
      tmps <- rep(0, nids - 1L)
      for (ks in seq_along(tmps)) {
        k <- ids[[ks]]
        tmps[[ks]] <- (dd[[k + 1L]] - dd[[k]]) * (fdd[[k + 1L]] + fdd[[k]])
      }
      res[[kl]] <- 1 / (dnew[[kl + 1L]] - dnew[[kl]]) * 1 / 2 *
        sum(tmps, na.rm = TRUE)
    }
  }
  
  res
}
