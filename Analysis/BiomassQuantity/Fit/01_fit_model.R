# Purpose:
# Fit cwexp model to biomass data via single-level CV to select lambda,
# then refit global model at selected lambda.
# For starters this is for testing on simulated data.
#
# Author: Martin Holdrege

# dependencies ------------------------------------------------------------

source('Functions/init.R')
source('Functions/read_write.R')
source('Functions/data/simulate_data.R')
source('Functions/models/cwexp_helpers.R')
source('Functions/models/cwexp_lognormal_tmb.R')
source('Functions/models/cwexp_lambda_path.R')
source('Functions/resampling/cv_metrics.R')
source('Functions/resampling/cv_main.R')
source('Functions/grouping/environmental_clustering.R')
source('Functions/grouping/make_folds.R')
source('Functions/sampling/select_purer.R')
library(future)
library(future.apply)

# params ------------------------------------------------------------------

test_run <- FALSE
pfts <- const$pfts
cover_cols <- paste0(pfts, "Cov")

# predictor variables (main effects only for now)
pred_vars <- c("tmean_CLIM",
               "precip_CLIM",
               "PrecipTempCorr_CLIM",
               "VPD_mean",
               "sand")

inter <- list(c("tmean_CLIM", "precip_CLIM"))

inter_terms <- purrr::map_chr(inter, ~ paste(.x, collapse = ":"))
all_terms <- c(pred_vars, inter_terms)
formula_full <- as.formula(
  paste("totalBio ~", paste(all_terms, collapse = " + "))
)

# --- config: all tunable settings in one place ----------------------------
config <- list(
  mv = 'm01',
  # model specification
  model = list(
    formula = formula_full,
    cover_cols = cover_cols,
    pred_vars = pred_vars
  ),
  # purer pixel selection
  purer = list(
    region_col = "region",
    q = 0.9,
    min_raw_cover = 0.05,
    max_n_per_region = NULL,  # NULL = no cap
    seed = 42
  ),
  # environmental clustering (for CV folds)
  clustering = list(
    vars = c("tmean_CLIM", "precip_CLIM", "VPD_mean","PrecipTempCorr_CLIM",
             "sand"),
    k = 20,
    seed = 1
  ),
  # cross-validation
  
  cv = list(
    n_folds = 4,
    en_alpha = 0.5,
    n_lambda = 15,
    include_zero = TRUE,
    select_rule = "1se",
    select_metric = "mae_log",
    fold_seed = 1
  ),
  # data
  data = list(
    version = if(opt$use_simulated) opt$vs else opt$vd,
    use_simulated = opt$use_simulated
  )
)

# read in data ------------------------------------------------------------

# reading output from DataPrep/07_make_analysis_ready.R
dat1 <- read_analysis_ready(opt)

if(test_run) {
  dat1 <- dplyr::sample_n(dat1, 1000)
}
# compile TMB model -------------------------------------------------------

dll_path <- "src/cwexp_lognormal_en_tmb.cpp"
dll_en <- cwexp_tmb_compile(dll_path, quiet = TRUE)

# find purer pixels -------------------------------------------------------

purer <- select_purer_by_region(
  dat = dat1,
  cover_vars = config$model$cover_cols,
  region_col = config$purer$region_col,
  q = config$purer$q,
  min_raw_cover = config$purer$min_raw_cover,
  max_n_per_region = config$purer$max_n_per_region,
  seed = config$purer$seed
)

dat_train <- purer$data
config$data$n_train <- nrow(dat_train)
cat("Purer pixel selection:\n")
cat("  Input rows:", nrow(dat1), "\n")
cat("  Selected rows:", nrow(dat_train), "\n")
cat("  Thresholds:\n")
print(purer$thresholds)
cat("\n")

# create environmental clusters and CV folds ------------------------------

clust <- make_env_clusters(
  data = dat_train,
  vars = config$clustering$vars,
  k = config$clustering$k,
  seed = config$clustering$seed
)
dat_train$env_cluster <- clust$env_cluster

folds <- make_cluster_folds(
  env_cluster = dat_train$env_cluster,
  n_folds = config$cv$n_folds,
  seed = config$cv$fold_seed
)


# seteup parallel ---------------------------------------------------------

parallel <- TRUE
workers <- min(config$cv$n_folds, parallel::detectCores(logical = FALSE)) 

# run inner CV to select lambda -------------------------------------------

cat("Running inner CV (", config$cv$n_folds, "folds,",
    config$cv$n_lambda, "lambdas)...\n")

plan(multisession, workers = workers)

inner_cv <- run_inner_cv(
  data = dat_train,
  folds = folds,
  en_alpha = config$cv$en_alpha,
  lambda_max_fun = cwexp_lambda_max_l1_tmb,
  lambda_max_args = list(
    formula = config$model$formula,
    cover_cols = config$model$cover_cols,
    dll_en = dll_en
  ),
  lambda_path_args = list(
    n_lambda = config$cv$n_lambda,
    include_zero = config$cv$include_zero
  ),
  fit_path_fun = cwexp_fit_lambda_path_tmb,
  fit_path_args = list(
    formula = config$model$formula,
    cover_cols = config$model$cover_cols,
    dll = dll_en,
    include_report = FALSE
  ),
  select_args = list(
    metric = config$cv$select_metric,
    rule = config$cv$select_rule
  ),
  keep_fold_results = TRUE,
  parallel = TRUE,
  dll_path = dll_path
)

plan(sequential)
selected_lambda <- inner_cv$selected$lambda[[1]]

cat("Selected lambda:", selected_lambda, "\n")
cat("Selection rule:", config$cv$select_rule, "\n")
cat("CV score summary:\n")
print(inner_cv$score_summary)
cat("\n")

# fit global model at selected lambda -------------------------------------

cat("Fitting global model at lambda =", selected_lambda, "...\n")

# use warm start from inner CV
start <- select_warm_start_par(inner_cv)

global_fit <- cwexp_fit_tmb(
  data = dat_train,
  formula = config$model$formula,
  cover_cols = config$model$cover_cols,
  dll = dll_en,
  penalty = "elastic_net",
  en_alpha = config$cv$en_alpha,
  lambda = selected_lambda,
  start = start,
  include_report = FALSE
)

cat("Global model convergence:", global_fit$tmb$opt$convergence, "\n")
cat("Estimated sigma:", round(global_fit$par$sigma, 3), "\n\n")

# save outputs ------------------------------------------------------------

v <- config$data$version

out <- list(
  fit = global_fit,
  cv = inner_cv,
  config = config,
  purer = list(
    thresholds = purer$thresholds,
    selection_summary = purer$selection_summary,
    n_selected = nrow(dat_train)
  ),
  clustering = list(
    k = config$clustering$k,
    vars = config$clustering$vars
  )
)

p_out <- file.path(
  paths$large,
  'Data_processed/BiomassQuantityData/Fit',
  paste0('fitted_model_', config$data$version, "_", config$mv, '.rds')
)

if(!test_run) {
  dir.create(dirname(p_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, p_out)
  cat("Saved fitted model to:", p_out, "\n")
}


