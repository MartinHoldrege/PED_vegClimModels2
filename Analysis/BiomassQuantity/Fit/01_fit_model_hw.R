# Purpose:
# Fit cwexp biomass models via single-level CV to select lambda,
# then refit global model at selected lambda. Supports fitting:
#   - Herbaceous model (RAP herb biomass ~ climate × herb cover, 1 group)
#   - Woody model (GEDI biomass ~ climate × tree/shrub cover, 2 groups)
# Each model is fit independently with its own data, formula, and cover_cols.
# Both are defined within a single model_spec entry (vm) which has $herb
# and $woody sub-specs.
#
# Input: training data .rds files created by 01_create_training_data.R
# Output: fitted model .rds files (one per model type)
#
# Author: Martin Holdrege

# NEXT: continue purer implementation
# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
library(future)
library(future.apply)

# params ------------------------------------------------------------------

test_run <- TRUE

# which models to fit (toggle on/off)
fit_herbaceous <- TRUE
fit_woody <- TRUE

# data version (from 08_create_training_data.R output)
vd <- "d05" # uses RAP cover, for code development/testing

# single model version — contains both herb and woody sub-specs
vm <- "m09"
vp <- 'p04'

stopifnot(vm %in% names(model_specs))
m_spec <- model_specs[[vm]]
p_spec <- purer_specs[[vp]]
stopifnot(all(c("herb", "woody") %in% names(m_spec)))

# CV settings (shared across models; override per-model if needed)
cv_settings <- list(
  n_folds = 4,
  en_alpha = 0.5,
  n_lambda = 10,
  include_zero = TRUE,
  select_rule = "1se",
  select_metric = "mae_log1p",
  fold_seed = 1
)

# clustering settings (for blocked CV folds)
clustering_settings <- list(
  k = 20,
  seed = 1
)

# optimizer control
control <- list(iter.max = 1000, eval.max = 1000)

# parallelization
use_parallel <- TRUE


# helper: fit one model (CV + global fit) ----------------------------------

fit_one_model <- function(dat_train, config, dll_en,
                          control, use_parallel,
                          purer_spec,
                          model_label = "model") {
  
  formula <- config$model$formula
  cover_cols <- config$model$cover_cols
  pred_vars <- config$model$pred_vars
  dll_path <- config$model$dll_path
  fix_alpha_pfts <- config$model$fix_alpha_pfts
  fix_alpha_filter <- config$model$fix_alpha_filter
  
  cat("\n============================================================\n")
  cat("Fitting:", model_label, "\n")
  cat("  Formula:", deparse(formula), "\n")
  cat("  Cover cols:", paste(cover_cols, collapse = ", "), "\n")
  cat("  N training:", nrow(dat_train), "\n")
  cat("============================================================\n\n")
  
  # collinearity check
  cat("Collinearity check:\n")
  print(check_collinearity(data = dat_train, formula = formula))
  cat("\n")
  
  # --- pre-estimate fixed alphas if configured ---
  fixed_alpha <- NULL
  alpha_prefit_info <- NULL
  

  # sample for fitting --------------------------------------

  
  selection <- select_training_pixels(
    dat = dat_train,
    cover_vars = cover_cols,
    purer_spec = purer_spec,
    region_col = "region"
  )
  
  dat_train <- selection$data
  
  # --- environmental clusters and CV folds ---
  cluster_vars <- pred_vars
  clust <- make_env_clusters(
    data = dat_train,
    vars = cluster_vars,
    k = config$clustering$k,
    seed = config$clustering$seed
  )
  dat_train$env_cluster <- clust$env_cluster
  
  folds <- make_cluster_folds(
    env_cluster = dat_train$env_cluster,
    n_folds = config$cv$n_folds,
    seed = config$cv$fold_seed
  )
  
  # --- parallel setup ---
  if (use_parallel) {
    workers <- min(config$cv$n_folds,
                   parallel::detectCores(logical = FALSE))
    plan(multisession, workers = workers)
  }
  
  # --- inner CV ---
  cat("Running inner CV (", config$cv$n_folds, "folds,",
      config$cv$n_lambda, "lambdas)...\n")
  
  inner_cv <- run_inner_cv(
    data = dat_train,
    folds = folds,
    en_alpha = config$cv$en_alpha,
    lambda_max_fun = cwexp_lambda_max_l1_tmb,
    lambda_max_args = list(
      formula = formula,
      cover_cols = cover_cols,
      dll_en = dll_en
    ),
    lambda_path_args = list(
      n_lambda = config$cv$n_lambda,
      include_zero = config$cv$include_zero
    ),
    fit_path_fun = cwexp_fit_lambda_path_tmb,
    fit_path_args = list(
      formula = formula,
      cover_cols = cover_cols,
      fixed_alpha = fixed_alpha,
      dll = dll_en,
      include_report = FALSE,
      control = list(iter.max = 500, eval.max = 500)
    ),
    select_args = list(
      metric = config$cv$select_metric,
      rule = config$cv$select_rule
    ),
    keep_fold_results = TRUE,
    parallel = use_parallel,
    dll_path = dll_path
  )
  
  if (use_parallel) plan(sequential)
  
  selected_lambda <- inner_cv$selected$lambda[[1]]
  cat("Selected lambda:", selected_lambda, "\n")
  cat("Selection rule:", config$cv$select_rule, "\n")
  cat("CV score summary:\n")
  print(inner_cv$score_summary)
  cat("\n")
  
  # --- fit global model at selected lambda ---
  cat("Fitting global model at lambda =", selected_lambda, "...\n")
  
  start <- select_warm_start_par(inner_cv)
  
  global_fit <- cwexp_fit_tmb(
    data = dat_train,
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en,
    penalty = "elastic_net",
    en_alpha = config$cv$en_alpha,
    lambda = selected_lambda,
    start = start,
    fixed_alpha = fixed_alpha,
    include_report = FALSE,
    control = control
  )
  
  cat("Global model convergence:", global_fit$tmb$opt$convergence, "\n")
  cat("Estimated sigma:", round(global_fit$par$sigma, 3), "\n")
  cat("Estimated alpha:\n")
  print(round(global_fit$par$alpha, 3))
  cat("\n")
  
  list(
    fit = global_fit,
    cv = inner_cv,
    config = config,
    alpha_prefit = alpha_prefit_info,
    clustering = list(
      k = config$clustering$k,
      vars = cluster_vars
    ),
    n_train = nrow(dat_train)
  )
}

# compile TMB -------------------------------------------------------------

dll_paths_needed <- unique(c(
  if (fit_herbaceous) m_spec$herb$dll_path,
  if (fit_woody) m_spec$woody$dll_path
))

dlls <- purrr::set_names(
  purrr::map_chr(dll_paths_needed, ~ cwexp_tmb_compile(.x, quiet = TRUE)),
  dll_paths_needed
)

# === HERBACEOUS MODEL =====================================================

if (fit_herbaceous) {
  
  sub_spec <- m_spec$herb

  # load training data
  p_herb <- file.path(
    paths$large,
    "Data_processed/BiomassQuantityData/analysis_ready",
    paste0("biomass_herb_sample_", vd, ".rds")
  )
  stopifnot(file.exists(p_herb))
  herb_obj <- readRDS(p_herb)
  
  dat_herb <- herb_obj$data
  if (test_run) dat_herb <- dplyr::sample_n(dat_herb, 3000)
  
  # build formula and config
  formula_herb <- make_model_formula(sub_spec)
  config_herb <- build_config(
    vm = vm, sub_spec = sub_spec,
    vp = vp,
    p_spec = p_spec$herb,
    formula = formula_herb,
    cv_settings = cv_settings,
    clustering_settings = clustering_settings,
    vd = vd, model_type = "herbaceous"
  )
  
  dll_en_herb <- dlls[[sub_spec$dll_path]]
  
  # fit
  herb_result <- fit_one_model(
    dat_train = dat_herb,
    config = config_herb,
    dll_en = dll_en_herb,
    control = control,
    use_parallel = use_parallel,
    model_label = "herbaceous_model",
    purer_specs
  )
  
  # save
  out_herb <- c(herb_result, list(scale = herb_obj$scale_df))
  
  suffix_herb <- if (test_run) "test_run_herb" else paste("herb", vd, vm, sep = "-")
  
  p_out_herb <- file.path(
    paths$large,
    "Data_processed/BiomassQuantityData/Fit",
    paste0("fitted_model_", suffix_herb, ".rds")
  )
  dir.create(dirname(p_out_herb), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out_herb, p_out_herb)
  cat("Saved herbaceous model to:", p_out_herb, "\n\n")
}

# === WOODY MODEL ==========================================================

if (fit_woody) {
  
  sub_spec <- m_spec$woody
  p_sub_spec <- pur
  # load training data
  p_woody <- file.path(
    paths$large,
    "Data_processed/BiomassQuantityData/analysis_ready",
    paste0("biomass_woody_sample_", vd, ".rds")
  )
  stopifnot(file.exists(p_woody))
  woody_obj <- readRDS(p_woody)
  
  dat_woody <- woody_obj$data
  if (test_run) dat_woody <- dplyr::sample_n(dat_woody, 3000)
  
  # build formula and config
  formula_woody <- make_model_formula(sub_spec)
  config_woody <- build_config(
    vm = vm, sub_spec = sub_spec,
    formula = formula_woody,
    cv_settings = cv_settings,
    clustering_settings = clustering_settings,
    vd = vd, model_type = "woody"
  )
  
  dll_en_woody <- dlls[[sub_spec$dll_path]]
  
  # fit
  woody_result <- fit_one_model(
    dat_train = dat_woody,
    config = config_woody,
    dll_en = dll_en_woody,
    control = control,
    use_parallel = use_parallel,
    model_label = "Woody model"
  )
  
  # save
  out_woody <- c(woody_result, list(scale = woody_obj$scale))
  
  suffix_woody <- if (test_run) "test_run_woody" else paste("woody", vd, vm, sep = "-")
  
  p_out_woody <- file.path(
    paths$large,
    "Data_processed/BiomassQuantityData/Fit",
    paste0("fitted_model_", suffix_woody, ".rds")
  )
  dir.create(dirname(p_out_woody), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out_woody, p_out_woody)
  cat("Saved woody model to:", p_out_woody, "\n\n")
}

cat("Done.\n")