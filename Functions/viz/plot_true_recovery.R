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
#' @param newdata Optional data frame. If NULL, uses `truth$data`. Must contain
#'   `{group}Bio` columns if `use_observed = "force"`.
#' @param weighted Logical. If TRUE, plots cover-weighted group contributions
#'   (only used when `use_observed = "never"`).
#' @param use_observed One of `"force"` (use `{group}Bio` columns from data,
#'   error if any absent) or `"never"` (predict from truth parameters).
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
                                use_observed = "never",                
                                density = FALSE,
                                bins = 30,
                                title = "Per-PFT biomass recovery",
                                point_alpha = 0.3,
                                point_size = 0.5) {
  stopifnot(use_observed %in% c("force", "never"))                    
  
  data        <- if (!is.null(newdata)) newdata else truth$data
  group_names <- names(truth$par$alpha)
  
  mu_true_grp <- if (use_observed == "force") {                       
    bio_cols <- paste0(str_replace(group_names, 'Cov$', ''), "Bio")                            
    missing  <- bio_cols[!bio_cols %in% names(data)]                  
    if (length(missing) > 0) {                                        
      stop("use_observed = 'force' but columns not found: ",          
           paste(missing, collapse = ", "))                           
    }                                                                 
    mat <- as.matrix(data[bio_cols])             
    colnames(mat) <- group_names       
    if(!weighted) { # dividing by cover to get 'biomass at 100% cover'
      stopifnot(stringr::str_detect(group_names, 'Cov'),
                group_names %in% names(data))
      mat <- mat/as.matrix(data[group_names])
    }
    mat                                                               
  } else {                                                            
    predict_by_group(truth, newdata = data, weighted = weighted)
  }                                                                   
  
  mu_hat_grp <- predict_by_group(fit, newdata = data, weighted = weighted)
  
  df <- dplyr::bind_rows(lapply(seq_along(group_names), function(g) {
    dplyr::tibble(
      group   = stringr::str_replace(group_names[g], 'Cov', ''),
      true_mu = mu_true_grp[, g],
      pred_mu = mu_hat_grp[, g]
    )
  }))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$true_mu, y = .data$pred_mu))
  
  if (density) {
    p <- p + ggplot2::geom_hex(bins = bins)+
      ggplot2::scale_fill_viridis_c(trans = "log", name = "count")
  } else {
    p <- p + ggplot2::geom_point(alpha = point_alpha, size = point_size)
  }
  
  x_lab <- if (use_observed == "force") "Observed (mean) per-PFT biomass" else "True per-PFT mu"
  y_lab <- "Predicted per-PFT mu"
  if(!weighted) {
   x_lab <- paste(x_lab, '(at 100% cover)')
   y_lab <- paste(y_lab, '(at 100% cover)')
  } 
  
    p +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed") +
    ggplot2::facet_wrap(~ group, scales = "free") +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = title
    )
}


#' Plot alpha recovery across simulation replicates
#'
#' @param results Tibble with columns PFT, alpha_est, alpha_true (from bias
#'   experiment).
#' @param title Character; plot title.
#'
#' @return A ggplot object.
plot_alpha_recovery <- function(results, 
                                title = "Alpha recovery across replicate simulations of data") {
  truth_df <- results |>
    dplyr::distinct(PFT, alpha_true)
  
  med <- results |> 
    dplyr::group_by(PFT) |> 
    summarize(median = median(alpha_est))
  ggplot(results, aes(x = 1, y = alpha_est)) +
    geom_boxplot(fill = "steelblue", alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
    geom_hline(aes(yintercept = alpha_true),
               linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_text(data = truth_df,
              aes(x = 0.8, y = -Inf, label = paste("true =", round(alpha_true, 1))),
              vjust = -0.5, size = 3, color = "red") +
    geom_text(data = med,
              aes(x = 0.8, y = -Inf, label = paste("median =\n", round(median, 1))),
              vjust = -2, size = 3, color = "blue") +
    facet_wrap(~ PFT, scales = "free_y") +
    labs(x = NULL, y = "Estimated alpha", title = title,
         caption = "Dashed red line = true alpha") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

