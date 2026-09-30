# 02_fit_classification.R
#
# Fits one forest / non-forest classification model for the cover pipeline
# and saves it.
# The model is chosen by --cover_type (family), --cover_response (which model
# in that family) and --vmc (cover model version), which index cover_specs;
# --vc is the cover data version. Defaults are in params.R; main_cover.R can
# run other combinations.
#
# spec$engine picks the fitting engine: "glmnet" (penalized logistic
# regression), "ranger" (random forest, as a benchmark for how well the same
# predictors can do without the constraint of an equation) or "grpreg"
# (hierarchical group lasso: x^2 can only enter with x, and x:z with x and
# z; see hier_groups()).
#
# spec$family picks the response. "binomial": the class (cover above
# spec$cover_threshold). "gaussian" or "poisson" (grpreg only): cover itself,
# after spec$response_transform; the predicted cover is then cut where the
# predicted fraction of forest matches the observed one, as a probability
# would be. "quasibinomial" (glmnet only): cover / 100 as a proportion, fit
# with glmnet's binomial likelihood on a two-column response cbind(1 - p, p)
# (the same estimates as quasibinomial); predictions are mean cover (%), cut
# the same way as the continuous families.
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
# quasibinomial: a two-column matrix of proportions (second column is the
# target), so the binomial deviance is computed on cover rather than the class.
# gaussian: cover (or log1p cover), with CV error as mean squared error on
# that scale
y_glmnet <- switch(spec$family,
                   binomial = obs,
                   quasibinomial = cbind(1 - y_fit, y_fit),
                   gaussian = y_fit,
                   stop("family not set up for glmnet: ", spec$family))
glmnet_family <- if (spec$family == "gaussian") "gaussian" else "binomial"

fit <- cv.glmnet(x = x, y = y_glmnet,
                 family = glmnet_family,
                 alpha = spec$alpha,
                 foldid = foldid,
                 keep = TRUE)  # keeps the out-of-fold predictions

lambda <- switch(spec$cv$select_rule,
                 "1se" = fit$lambda.1se,
                 "min" = fit$lambda.min,
                 stop("unknown select_rule: ", spec$cv$select_rule))

# probability of forest, or (quasibinomial, gaussian) % cover
pred_in <- predict_score(list(fit = fit, lambda = lambda,
                              config = list(spec = spec)), x)

# fit$fit.preval: out-of-fold predictions on the link scale, one column per
# lambda. From the same folds that chose lambda, so mildly optimistic.
eta_oof <- fit$fit.preval[, match(lambda, fit$lambda)]
pred_oof <- switch(spec$family,
                   binomial = plogis(eta_oof),
                   quasibinomial = plogis(eta_oof) * 100,
                   gaussian = .untransform_cover(eta_oof,
                                                 spec$response_transform))
rm(eta_oof)

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
  
} else if (spec$engine == "grpreg") {
  
  # each term gets a group holding it and the terms it requires (x^2 needs x;
  # x:z needs x and z; log1p_x needs nothing), with shared columns repeated
  # once per group. With a group lasso, a square or interaction can then only
  # enter along with the terms it needs (latent overlapping groups, as in
  # glinternet)
  groups <- hier_groups(colnames(x))
  
  # grpreg standardizes and orthonormalizes each group internally, and
  # returns coefficients on the scale of x. cv.grpreg() needs folds 1..k
  fold <- as.integer(factor(foldid))
  # poisson fits warn "non-integer x" (from an AIC grpreg computes and
  # doesn't use); harmless
  cv_fit <- do.call(
    grpreg::cv.grpreg,
    c(list(X = x[, groups$term, drop = FALSE],
           y = y_fit,
           group = factor(groups$group, levels = unique(groups$group)),
           family = spec$family,
           penalty = "grLasso",
           alpha = spec$alpha,
           fold = fold,
           returnY = TRUE),  # out-of-fold predictions for every lambda
      spec$grpreg)
  )
  
  lambda <- cv_fit$lambda.min
  
  if (cv_fit$min == length(cv_fit$lambda)) {
    warning("CV error is lowest at the smallest lambda: the path may have",
            " ended early (raise max.iter in spec$grpreg) or need a smaller",
            " lambda.min")
  }
  
  # the fitted model, minus the parts that grow with the data (fitted values
  # for every row and lambda, and the response)
  grpreg_fit <- cv_fit$fit
  grpreg_fit$linear.predictors <- NULL
  grpreg_fit$y <- NULL
  
  # one coefficient per column of x (the sum over its copies), for the
  # coefficient table and the hierarchy check
  b <- coef(grpreg_fit, lambda = lambda)
  tmp <- tibble(coef = b, var = names(b)) |> 
    summarise(coef = sum(coef), .by = var)
  coefs <- tmp$coef
  names(coefs) <- tmp$var

  # every term in the model has the terms it requires
  in_model <- names(coefs)[-1][coefs[-1] != 0]
  required <- groups$term[groups$group %in% in_model]
  stopifnot("a selected term is missing a term it requires" =
              all(coefs[required] != 0))
  
  # cv: mean out-of-fold deviance per lambda (and its SE), for the report
  fit <- list(grpreg = grpreg_fit,
              groups = groups,
              coef = coefs,
              cv = tibble(lambda = cv_fit$lambda, cvm = cv_fit$cve,
                          cvsd = cv_fit$cvse))
  
  pred_in <- predict_score(list(fit = fit, lambda = lambda,
                                config = list(spec = spec)), x)
  
  # out-of-fold predictions at the selected lambda, from the same folds that
  # chose it (so mildly optimistic), back on the % cover scale
  pred_oof <- .untransform_cover(cv_fit$Y[, cv_fit$min],
                                 spec$response_transform)
  rm(cv_fit)
  
} else {
  stop("unknown engine: ", spec$engine)
}

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