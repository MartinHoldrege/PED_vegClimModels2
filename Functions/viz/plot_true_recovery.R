# functions for plotting how well true results
# are recovered by a model, (for when using simulated data
# where 'truth' is known)


#' Plot per-group mu recovery: true vs predicted (i.e. mean true biomass
#' a given pft cover and given covariate values)
#'
#' Creates a faceted scatterplot (or hex density plot) with one panel per
#' functional group, showing true per-group mu on the x-axis and predicted
#' per-group mu on the y-axis.
#'
#' @param truth Object with true parameters (`cwexp_dummy` or similar).
#' @param fit A single fitted model object.
#' @param newdata Optional data frame. If NULL, uses `truth$data`.
#' @param weighted Logical. If TRUE, plots cover-weighted group contributions.
#' @param density Logical. If TRUE, uses hex bin density instead of points.
#' @param bins Integer; number of hex bins per axis (used when `density = TRUE`).
#' @param title Optional plot title.
#' @param point_alpha Numeric; transparency for points (used when `density = FALSE`).
#' @param point_size Numeric; size for points (used when `density = FALSE`).
#'
#' @return A ggplot object.
#' @export
plot_group_recovery <- function(truth, fit, newdata = NULL,
                                weighted = TRUE,
                                density = FALSE,         
                                bins = 30,               
                                title = "Per-PFT mu recovery",
                                point_alpha = 0.3,
                                point_size = 0.5) {
  data <- if (!is.null(newdata)) newdata else truth$data
  
  mu_true_grp <- predict_by_group(truth, newdata = data, weighted = weighted)
  mu_hat_grp  <- predict_by_group(fit, newdata = data, weighted = weighted)
  
  group_names <- colnames(mu_true_grp)
  
  df <- dplyr::bind_rows(lapply(seq_along(group_names), function(g) {
    dplyr::tibble(
      group   = group_names[g],
      true_mu = mu_true_grp[, g],
      pred_mu = mu_hat_grp[, g]
    )
  }))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$true_mu, y = .data$pred_mu))
  
  if (density) {                                                   
    p <- p + ggplot2::geom_hex(bins = bins)                       
  } else {                                                        
    p <- p + ggplot2::geom_point(alpha = point_alpha, size = point_size)
  }                                                                
  
  p +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::facet_wrap(~ group, scales = "free") +
    ggplot2::labs(
      x     = "True per-PFT mu",
      y     = "Predicted per-PFT mu",
      title = title
    )
}
