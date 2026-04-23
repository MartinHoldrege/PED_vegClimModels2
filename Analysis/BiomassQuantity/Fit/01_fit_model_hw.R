# Purpose:
# Fit cwexp biomass models via single-level CV to select lambda,
# then refit global model at selected lambda. Supports fitting:
#   - Herbaceous model (RAP herb biomass ~ climate × herb cover, 1 group)
#   - Woody model (GEDI biomass ~ climate × tree/shrub cover, 2 groups)
# Each model is fit independently with its own data, formula, and cover_cols.
# Both are defined within a single model_spec entry (vm) which has $herb
# and $woody sub-specs.
#
# Input: training data .rds files created by 
# '06_sample_training_data.R'
# Output: fitted model .rds files (one per model type)
#
# Author: Martin Holdrege

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
library(future)
library(future.apply)

# params ------------------------------------------------------------------

test_run <- FALSE

# which models to fit (toggle on/off)
fit_herbaceous <- TRUE
fit_woody <- TRUE

# data version (from 08_create_training_data.R output)
vd <- opt$vd # uses RAP cover, for code development/testing

# single model version — contains both herb and woody sub-specs
vm <- opt$vm
vp <- opt$vp

suffix0 <- paste(vd, vp, vm, sep = "-")
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
    purer_spec = p_spec$herb,
    formula = formula_herb,
    cv_settings = cv_settings,
    clustering_settings = clustering_settings,
    vd = vd, model_type = "herb"
  )
  
  dll_en_herb <- dlls[[sub_spec$dll_path]]
  
  # fit
  herb_result <- fit_one_model_cwexp(
    data = dat_herb,
    config = config_herb,
    dll_en = dll_en_herb,
    control = control,
    use_parallel = use_parallel,
    model_label = "Herbaceous model"
  )
  
  # save
  out_herb <- c(herb_result, list(scale = herb_obj$scale_df))
  
  suffix_herb <- if (test_run) "test_run_herb" else paste0("herb_", suffix0)
  
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
    purer_spec = p_spec$woody,
    vp = vp,
    vd = vd, model_type = "woody"
  )
  
  dll_en_woody <- dlls[[sub_spec$dll_path]]
  
  # fit
  woody_result <- fit_one_model_cwexp(
    data = dat_woody,
    config = config_woody,
    dll_en = dll_en_woody,
    control = control,
    use_parallel = use_parallel,
    model_label = "Woody model"
  )
  
  # save
  out_woody <- c(woody_result, list(scale = woody_obj$scale))
  
  suffix_woody <- if (test_run) "test_run_woody" else paste0("woody_", suffix0)
  
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