

plot_obs_vs_pred <- function(data, observed = 'observed', predicted = 'predicted',
                             size = 1e5, title = NULL) {
  
  size <- if(size < nrow(data)) size else nrow(data)
  data2 <- dplyr::sample_n(data, size = size)
  g1 <- ggplot(data2, aes(.data[[observed]],
                          .data[[predicted]])) +
    geom_point(alpha = 0.05, size = 0.3) +
    geom_smooth(se = FALSE) +
    geom_abline(slope = 1, intercept = 0, color = 'red')
  
  g2 <- g1 +
    scale_y_continuous(transform='log1p') +
    scale_x_continuous(transform ='log1p') +
    labs(x = 'Observed (log scale)',
         y = 'Predicted (log scale)')
  
  g1b <- g1 +
    labs(x = 'Observed', y = 'Predicted')
  
  g1b + g2 + patchwork::plot_annotation(title = title)
}


plot_residual <- function(data, observed = 'observed', predicted = 'predicted',
                             size = 1e5) {
  
  size <- if(!is.null(size) && size < nrow(data)) size else nrow(data)
  
  data2 <- data
  data2$residual <- log1p(data2[[observed]]) - log1p(data2[[predicted]])
  data3 <- dplyr::sample_n(data2, size = size)
  g1 <- ggplot(data3, aes(.data[[predicted]],
                          .data$residual))+
    geom_point(alpha = 0.05, size = 0.3) +
    geom_smooth(se = FALSE) +
    geom_abline(slope = 0, intercept = 0, color = 'red',
                linetype = 2) +
    labs(x = 'Predicted',
         y = 'residual (log1pscale)',
         subtitle = 'Residual vs predicted')
  
  # using full sized data set for histogram (not sampled) 
  g2 <- ggplot(data2, aes(.data$residual))+
    geom_histogram(bins = 100)+
    labs(x = 'residual (log1p scale)',
         subtitle = 'Residual distribution')
  
  g1 + g2
}

# log likelihood


#' Plot log-likelihood against observed and predicted biomass
#'
#' Creates a two-panel plot showing per-observation log-likelihood against
#' observed (left) and predicted (right) values. Optionally applies `log1p`
#' scaling to the x-axis.
#'
#' @param data Data frame with observed and predicted columns.
#' @param fit Fitted cwexp model object (passed to `cwexp_loglik`).
#' @param observed Character; name of observed response column. (as created by predict)
#' @param predicted Character; name of predicted column.
#' @param response Character; name of the response column in `data` for
#'   `cwexp_loglik`. Default same as `observed`.
#' @param shift Numeric; shift parameter for log-likelihood. Default 1.
#' @param log1p_x Logical; if TRUE, apply `log1p` to the x-axis. Default FALSE.
#' @param size Integer or NULL; subsample size for the scatterplots. If NULL
#'   or greater than `nrow(data)`, use all rows.
#'
#' @return A patchwork object (two ggplots side by side).
#' @export
plot_loglik <- function(data, fit,
                        response = NULL,
                        shift = NULL,
                        log1p_x = FALSE,
                        size = 1e5) {
  
  data2 <- data
  
  if (is.null(response)) {
    stopifnot(length(fit$spec$formula) == 3)
    response <-  fit$spec$formula[2] |> as.character()
  }
  stopifnot(response %in% names(data))
  
  dll <- fit$spec$dll
  if(is.null(shift) & !is.null(dll)) {
    shift <- if(dll == "cwexp_lognormal_en_tmb") {
      0
    } else {
      1
    }
    
  } else if (is.null(shift)){
    shift <- 1
  }
  data2$.ll <- cwexp_loglik(fit, data, response = response, shift = shift)
  data2$.predicted <- predict(fit, data, type = "mu")
  size <- if (!is.null(size) && size < nrow(data2)) size else nrow(data2)
  data3 <- dplyr::sample_n(data2, size = size)
  
  x_lab_suffix <- if (log1p_x) " (log1p scale)" else ""
  
  # observed vs log-likelihood
  g1 <- ggplot2::ggplot(data3, ggplot2::aes(
    x = if (log1p_x) log1p(.data[[response]]) else .data[[response]],
    y = .data$.ll
  )) +
    ggplot2::geom_point(alpha = 0.05, size = 0.3) +
    ggplot2::geom_smooth(se = FALSE) +
    ggplot2::labs(
      x = paste0("Observed", x_lab_suffix),
      y = "Log-likelihood",
      subtitle = "Observed vs log-likelihood"
    )
  
  # predicted vs log-likelihood
  g2 <- ggplot2::ggplot(data3, ggplot2::aes(
    x = if (log1p_x) log1p(.data[['.predicted']]) else .data[['.predicted']],
    y = .data$.ll
  )) +
    ggplot2::geom_point(alpha = 0.05, size = 0.3) +
    ggplot2::geom_smooth(se = FALSE) +
    ggplot2::labs(
      x = paste0("Predicted", x_lab_suffix),
      y = "Log-likelihood",
      subtitle = "Predicted vs log-likelihood"
    )
  
  g1 + g2
}
