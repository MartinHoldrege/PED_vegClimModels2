# 02_fit_classification.R
#
# Fits one binomial classification model for the cover pipeline and saves it.
# The model is chosen by --cover_type (family), --cover_response (which model
# in that family) and --vmc (cover model version), which index cover_specs;
# --vc is the cover data version. Defaults are in params.R; main_cover.R can
# run other combinations.
#
# spec$engine picks the fitting engine: "glmnet" (penalized logistic
# regression) or "ranger" (random forest, as a benchmark for how well the same
# predictors can do without the constraint of an equation).
#
# Predictors: log1p (where named "log1p_") -> standardized with the fixed
# global parameters -> squares and interactions (prepare_predictors()), so
# coefficients mean the same thing wherever the model is applied. Ranger specs
# turn off log1p, squares and interactions, so the same function returns the
# main effects alone.
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
dat <- switch(
  spec$rows,
  all = dat,
  stop("unknown spec$rows: ", spec$rows)
)

# rows with the response and every source column; log1p_MAP comes from MAP
source_vars <- unique(c(str_remove(spec$pred_vars, "^log1p_"),
                        spec$cv$cluster_vars))
dat <- dat[complete.cases(dat[c(spec$response, source_vars)]), ]

# response: 1 where cover is above the threshold (e.g. forest). Called obs,
# not y, because y is also the cell-centre coordinate column
obs <- as.integer(dat[[spec$response]] > spec$cover_threshold)

scale_df <- read_scale_params()  # climate and soils; no anomalies here

x <- prepare_predictors(dat,
                        pred_vars = spec$pred_vars,
                        scale_df = scale_df,
                        squares = spec$squares,
                        interactions = spec$interactions)
stopifnot(nrow(x) == length(obs))

# folds -------------------------------------------------------------------

# environmental blocking: k-means clusters in climate space (the variables
# are standardized inside make_env_clusters()), one fold per cluster
clusters <- make_env_clusters(dat,
                              vars = spec$cv$cluster_vars,
                              iter.max = 500,
                              k = spec$cv$k_clusters,
                              seed = 1,
                              # forcing cells with multiple years to all go to the same fold
                              group = dat$cell)
foldid <- clusters$env_cluster

stopifnot(length(foldid) == nrow(dat), !anyNA(foldid))

# fit ---------------------------------------------------------------------

if (spec$engine == "glmnet") {
  
# standardize = TRUE (the default) rescales every column internally before
# penalizing, so squares and interactions are penalized on the same footing;
# the coefficients returned are on the scale of x
fit <- cv.glmnet(x = x, y = obs,
                 family = "binomial",
                 alpha = spec$alpha,
                 foldid = foldid,
                 keep = TRUE)  # keeps the out-of-fold predictions

lambda <- switch(spec$cv$select_rule,
                 "1se" = fit$lambda.1se,
                 "min" = fit$lambda.min,
                 stop("unknown select_rule: ", spec$cv$select_rule))

pred_in <- as.numeric(predict(fit, newx = x, s = lambda, type = "response"))

# fit$fit.preval: out-of-fold predictions on the link scale, one column per
# lambda. From the same folds that chose lambda, so mildly optimistic.
pred_oof <- plogis(fit$fit.preval[, match(lambda, fit$lambda)])

  # fit.preval holds out-of-fold predictions for every lambda (large); dropped
  # now that the column at the selected lambda has been taken
  fit$fit.preval <- NULL
  
} else if (spec$engine == "ranger") {
  
  lambda <- NA_real_  # no penalty to select
  
  fit_rf <- function(x, y) {
    do.call(ranger::ranger,
            c(list(x = x, y = factor(y, levels = c(0, 1)),
                   probability = TRUE, seed = 1),
              spec$ranger))
  }
  
  fit <- fit_rf(x, obs)
  
  # out-of-bag, not in-sample: in-sample forest predictions are near-perfect,
  # so the threshold taken from them would be meaningless
  pred_in <- fit$predictions[, "1"]
  
  # do.call() stores the evaluated arguments -- the full training matrix -- in
  # fit$call, and fit$predictions duplicates pred_in; predict() needs neither
  fit$call <- NULL
  fit$predictions <- NULL
  
  # out-of-fold predictions from the same environmental folds the lasso uses,
  # so the two engines are scored the same way. Costs one forest per fold.
  pred_oof <- rep(NA_real_, length(obs))
  for (f in unique(foldid)) {
    i <- foldid == f
    pred_oof[i] <- predict(fit_rf(x[!i, , drop = FALSE], obs[!i]),
                           data = x[i, , drop = FALSE])$predictions[, "1"]
  }
  stopifnot(!anyNA(pred_oof))
  
} else {
  stop("unknown engine: ", spec$engine)
}

# predictions -------------------------------------------------------------

# cutoff from pred_in (in-sample for glmnet, out-of-bag for ranger)
threshold <- choose_threshold(obs = obs, pred = pred_in,
                              method = spec$threshold_method)

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

class(out) <- "cover_classification"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, out_file)