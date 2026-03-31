# plotting recovery of alpha values

source("Functions/init.R")
source_functions()
# params ------------------------------------------------------------------

vd <- 's07'
vp <- 'p01'
vm <- 'm08'

use_lambda <- TRUE 
# reading in data ---------------------------------------------------------

suffix <- paste0(vd, '-', vp, '-',vm)

suffix2 <- if(use_lambda){
  paste0(suffix, '_reg')
} else{
  suffix
}
# created in"Analysis/BiomassQuantity/Fit/test_cwexp_parameter_bias.R"

p <- file.path(paths$large,
          "Data_processed/BiomassQuantityData/Fit",
          paste0("param_bias_experiment_", suffix2, ".rds"))

obj <- readRDS(p)


# process -----------------------------------------------------------------
if(use_lambda) {
  lambda <- obj$spec$lambda
} else {
  lambda <- 0
}

results <- obj$results

B_true <- obj$true_par$B


results_B <- results |> 
  mutate(B = map2(B, PFT, \(x, pft) {
    est <- x[ , pft]
    true <- B_true[, pft]
    stopifnot(names(true) == names(est))
    tibble(B_true = true, B_est = est)
  })) |> 
  select(rep, PFT, B) |> 
  unnest(cols = B)
                  

# create plots ------------------------------------------------------------

# alpha recovery
g <- plot_alpha_recovery(obj$results) +
  labs(caption = paste0('Models fit w/ lambda =', signif(lambda, digits = 2),
                        '\n', suffix2))

p_out <- paste0("Figures/BiomassQuantity/simulated/alpha_recovery_", 
                suffix, '.png')
ggsave(filename = p_out,
       plot = g,
       height = 7, width = 5)

# B recovery
g2 <- ggplot(results_B, aes(B_true, B_est, group = B_true)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~PFT, scales = 'free') +
  geom_abline(slope = 1) +
  labs(subtitle = 'Beta recovery across replicate simulations of data',
       caption = paste0('Models fit w/ lambda =', signif(lambda, digits = 2),
                       '\n', suffix2))
p_out2 <- paste0("Figures/BiomassQuantity/simulated/beta_recovery_", 
                suffix, '.png')
ggsave(filename = p_out2,
       plot = g2,
       height = 4, width = 5)
