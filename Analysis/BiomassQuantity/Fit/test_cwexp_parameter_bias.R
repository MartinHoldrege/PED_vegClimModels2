# Purpose:
# Test whether the cwexp MLE produces biased estimates of alpha
# for small-cover PFTs. Fits the correctly-specified model to 20
# datasets generated from the same true parameters with different
# noise draws. Lambda = 0 (pure MLE, no regularization).
#
# Author: Martin Holdrege

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

vd <- 's06'
vp <- 'p01'
vm <- 'm01'

m_spec <- model_specs[[vm]]
p_spec <- purer_specs[[vp]]

n_reps <- 20
# setup -------------------------------------------------------------------


# --- load the base data (real cover + climate, pre-normalized) ---
p3 <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated',
                paste0('simBiomass_', vd, '.rds'))

d <- readRDS(p3)
pred_vars <- d$spec$pred_vars
cover_cols <- d$spec$cover_cols
formula_full <- d$spec$formula
dat_base <- d$data

# making sure that if we label output
# w/ the model version, that is correct (i.e)
# the model version has the same params as the data generating
# process
stopifnot(
  m_spec$pred_vars == pred_vars,
  m_spec$pfts == str_replace(cover_cols, 'Cov', '')
)


B <- d$par$B
intercepts <- d$par$alpha
names(intercepts) <- str_replace(names(intercepts), 'Cov$', '')
sigma_obs <- d$par$sigma_obs

# putting into format for sim_bio
coefs <- lapply(1:nrow(B), function(i) {
  coef <- B[i, ]
  names(coef) <- str_replace(names(coef), 'Cov$', '')
  list(var = rownames(B)[[i]],
       coef = coef)
})

m_spec$dll_path
dll <- cwexp_tmb_compile(m_spec$dll_path, quiet = TRUE)
# run experiment ----------------------------------------------------------

results <- vector("list", n_reps)

for (i in seq_len(n_reps)) {
  cat("Rep", i, "of", n_reps, "...")

  # generate new data with different noise (same X, C, true params)
  set.seed(i * 100)
  sim_i <- sim_bio(
    data = select(dat_base, -matches('Bio$'), -totalMu),
    coefs = coefs,
    intercepts = intercepts,
    pred_vars = pred_vars,
    sigma_pft = 0,
    sigma_obs = sigma_obs,
    sigma_region = 0,
    normalize = FALSE
  )

  dat_i <- sim_i$data

  data_i_purer <- select_purer_by_region(
    dat = dat_i,
    cover_vars = cover_cols,
    region_col = "region",
    q = p_spec$q,
    min_raw_cover = p_spec$min_raw_cover,
    seed = 123 # same seed so same purer selection
  )
  
  # fit at lambda = 0 (pure MLE)
  fit_i <- cwexp_fit_tmb(
      data = data_i_purer$data,
      formula = formula_full,
      cover_cols = cover_cols,
      dll = dll,
      penalty = "elastic_net",
      en_alpha = 0.5,
      lambda = 0,
      include_report = FALSE,
      control = list(iter.max = 1000, eval.max = 1000)
    )
    
  results[[i]] <- tibble::tibble(
    rep = i,
    convergence = fit_i$tmb$opt$convergence,
    sigma = fit_i$par$sigma,
    PFT = names(fit_i$par$alpha),
    alpha_est = unname(fit_i$par$alpha),
    alpha_true = unname(intercepts),
    B = list(fit_i$par$B)
  )
}

# summarize ---------------------------------------------------------------

res_df <- dplyr::bind_rows(results)

# save --------------------------------------------------------------------

saveRDS(list(results = res_df,
             true_par = d$par, 
             spec = d$spec),
        file.path(paths$large,
                  "Data_processed/BiomassQuantityData/Fit",
                  paste0("param_bias_experiment_", vd, '-', vp, '-',
                         vm, ".rds")))
