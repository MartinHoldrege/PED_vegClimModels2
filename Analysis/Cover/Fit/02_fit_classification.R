# 02_fit_classification.R
#
# Fits one forest / non-forest classification model for the cover pipeline
# and saves it.
# The model is chosen by --cover_type (family), --cover_response (which model
# in that family) and --vmc (cover model version), which index cover_specs;
# --vc is the cover data version. Defaults are in params.R; main_cover.R can
# run other combinations.
#
# The folds and fitting are shared with 02_fit_cover.R (fit_engine.R).
#
# spec$engine picks the fitting engine: "glmnet" (penalized logistic
# regression), "ranger" (random forest, as a benchmark for how well the same
# predictors can do without the constraint of an equation), "grpreg"
# (hierarchical group lasso: x^2 can only enter with x, and x:z with x and
# z; see hier_groups()) or "gam" (one smooth per predictor, smoothing tuned
# by CV; see gam_helpers.R).
#
# spec$family picks the response. "binomial": the class (cover above
# spec$cover_threshold). "gaussian" or "poisson" (grpreg only): cover itself,
# after spec$response_transform; the predicted cover is then cut where the
# predicted fraction of forest matches the observed one, as a probability
# would be. "quasibinomial" (glmnet or gam): cover / 100 as a proportion, fit
# with glmnet's binomial likelihood on a two-column response cbind(1 - p, p)
# (the same estimates as quasibinomial), or mgcv's quasibinomial; predictions are mean cover (%), cut
# the same way as the continuous families.
#
# Predictors: log1p (where named "log1p_") -> standardized with the fixed
# global parameters -> squares and interactions (prepare_predictors()), so
# coefficients mean the same thing wherever the model is applied. Ranger specs
# turn off log1p, squares and interactions, so the same function returns the
# main effects alone; so do gam specs, which smooth the products of pairs
# in the formula instead (gam_formula()).
#
# Inputs:
#   cover_clim_soils_<vc>.csv  - 08_combine_cover_and_covariates.R
#   scale_params_*.csv         - 09_compute_scale_params.R
#   cover_specs                - Functions/models/cover_specs.R
#
# Output:
#   Data_processed/CoverData/Fit/<cover_type>_<cover_model>_<vc>-<vmc>.rds
#   (path from cover_model_path())
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

library(glmnet)

# spec --------------------------------------------------------------------
test_run <- FALSE
spec <- cover_specs[[opt$cover_type]][[opt$cover_model]][[opt$vmc]]

if (is.null(spec)) {
  stop("no spec for cover_specs$", opt$cover_type, "$", opt$cover_model, "$",
       opt$vmc)
}
stopifnot(opt$cover_type == "classification")

out_file <- cover_model_path(opt$cover_type, opt$cover_model, opt$vc,
                             opt$vmc)

# data --------------------------------------------------------------------

# original units; prepare_predictors() does the transforming and scaling
dat <- read_cover_training(opt$vc)

# which pixel-years to fit to
dat <- select_rows(dat, spec$rows)

# rows with the response and every source column; log1p_MAP comes from MAP
source_vars <- unique(c(str_remove(spec$pred_vars, "^log1p_"),
                        spec$cv$cluster_vars))
dat <- dat[complete.cases(dat[c(spec$response, source_vars)]), ]

if(test_run) {
  dat <- sample_n(dat, 1000)
}
# observed class: 1 where cover is above the threshold (e.g. forest). Called
# obs, not y, because y is also the cell-centre coordinate column
cover <- dat[[spec$response]]
obs <- as.integer(cover > spec$cover_threshold)

# what the model is fit to: the class, or (continuous families) cover itself
if (spec$family == "binomial") {
  y_fit <- obs
} else if (spec$family == "quasibinomial") {
  stopifnot(all(cover >= 0 & cover <= 100))
  y_fit <- cover / 100
} else {
  stopifnot(all(cover >= 0))
  y_fit <- switch(spec$response_transform,
                  identity = cover,
                  log1p = log1p(cover))
}

scale_df <- read_scale_params()  # climate and soils; no anomalies here

x <- prepare_predictors(dat,
                        pred_vars = spec$pred_vars,
                        scale_df = scale_df,
                        squares = spec$squares,
                        interactions = spec$interactions,
                        interact_log1p = spec$interact_log1p)
stopifnot(nrow(x) == length(obs))

# folds -------------------------------------------------------------------

# environmental blocking (see make_cover_folds())
folds <- make_cover_folds(dat, spec)
foldid <- folds$foldid
clusters <- folds$clusters
rm(folds)

# fit ---------------------------------------------------------------------

# engine and family from the spec (see fit_engine()). pred_in: in-sample
# (out-of-bag for ranger); pred_oof: out-of-fold
res <- fit_engine(spec, x = x, y_fit = y_fit, foldid = foldid)
fit <- res$fit
lambda <- res$lambda
pred_in <- res$pred_in
pred_oof <- res$pred_oof
rm(res)

# predictions -------------------------------------------------------------

# cutoff from pred_in (in-sample for glmnet and grpreg, out-of-bag for
# ranger). Predicted cover isn't a probability, so for continuous families
# the cutoff is the prediction quantile that matches the observed prevalence
threshold <- if (spec$family == "binomial") {
  choose_threshold(obs = obs, pred = pred_in, method = spec$threshold_method)
} else {
  prevalence_cutoff(obs = obs, score = pred_in)
}

# save --------------------------------------------------------------------

out <- list(
  fit = fit,
  lambda = lambda,
  threshold = threshold,
  config = list(cover_type = opt$cover_type,
                cover_model = opt$cover_model,
                vc = opt$vc, vmc = opt$vmc, spec = spec,
                x_colnames = colnames(x)),
  scale_df = scale_df,
  clustering = clusters,  # centers and scaling, for assign_to_clusters()
  predictions = tibble(cell = dat$cell, year = dat$year,
                       pred_in = pred_in, pred_oof = pred_oof)
)

if(!test_run) {
  class(out) <- "cover_classification"
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, out_file)
}