

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
