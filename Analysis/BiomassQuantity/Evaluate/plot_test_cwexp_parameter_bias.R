# plotting recovery of alpha values

source("Functions/init.R")

# params ------------------------------------------------------------------

vd <- 's07'
vp <- 'p01'
vm <- 'm03'

# reading in data ---------------------------------------------------------

suffix <- paste0(vd, '-', vp, '-',vm)
# created in"Analysis/BiomassQuantity/Fit/test_cwexp_parameter_bias.R"

p <- file.path(paths$large,
          "Data_processed/BiomassQuantityData/Fit",
          paste0("param_bias_experiment_", suffix, ".rds"))

obj <- readRDS(p)


# create plots ------------------------------------------------------------

ggplot(obj$results, aes(x = alpha_true, y = alpha_est)) +
  geom_point()
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
g <- plot_alpha_recovery(obj$results)
g
p_out <- paste0("Figures/BiomassQuantity/simulated/alpha_recovery_", 
                suffix, '.png')
ggsave(filename = p_out,
       plot = g,
       height = 7, width = 5)
