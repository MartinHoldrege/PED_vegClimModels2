# Fitting shared by the cover models: CV folds, and the fitting engines
# (glmnet, ranger, grpreg, gam). Used by 02_fit_classification.R and
# 03_fit_cover.R. Moved from 02_fit_classification.R unchanged, except that
# the response is always y_fit (for binomial models, the 0/1 class), and
# ranger can also fit a regression forest (family "gaussian").





#' Fit a cover model with the engine in its spec
#'
#' @param spec Model spec (from `cover_specs`): `spec$engine` picks the
#'   engine, `spec$family` the response.
#' @param x Design matrix from `prepare_predictors()`.
#' @param y_fit What the model is fit to: the 0/1 class (binomial), cover / 100
#'   (quasibinomial), or cover after `spec$response_transform` (gaussian,
#'   poisson).
#' @param foldid Fold of each row of `x`, from `make_cover_folds()`.
#' @return List: `fit` (what `predict_score()` uses), `lambda` (the selected
#'   penalty, or gam smoothing multiplier; NA for ranger), `pred_in`
#'   (in-sample; out-of-bag for ranger) and `pred_oof` (out-of-fold).
#'   Predictions are probabilities for binomial models, % cover otherwise.
#' @examples
#' dat <- read_cover_training("c02")
#' spec <- cover_specs$classification$forest$m01
#' x <- prepare_predictors(dat, pred_vars = spec$pred_vars,
#'                         scale_df = read_scale_params(),
#'                         squares = spec$squares,
#'                         interactions = spec$interactions,
#'                         interact_log1p = spec$interact_log1p)
#' y_fit <- as.integer(dat$cov_tree > spec$cover_threshold)
#' foldid <- make_cover_folds(dat = dat, spec = spec)$foldid
#' res <- fit_engine(spec = spec, x = x, y_fit = y_fit, foldid = foldid)
fit_engine <- function(spec, x, y_fit, foldid) {
  stopifnot(nrow(x) == length(y_fit), length(foldid) == length(y_fit))
  
  if (spec$engine == "glmnet") {
    
    # standardize = TRUE (the default) rescales every column internally before
    # penalizing, so squares and interactions are penalized on the same footing;
    # the coefficients returned are on the scale of x
    # quasibinomial: a two-column matrix of proportions (second column is the
    # target), so the binomial deviance is computed on cover rather than the class.
    # gaussian: cover (or log1p cover), with CV error as mean squared error on
    # that scale
    y_glmnet <- switch(spec$family,
                       binomial = y_fit,
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
    
    # binomial: probability forest on the class; gaussian: regression forest
    # on cover (after spec$response_transform)
    classify <- spec$family == "binomial"
    fit_rf <- function(x, y, args = spec$ranger) {
      if (classify) y <- factor(y, levels = c(0, 1))
      do.call(ranger::ranger,
              c(list(x = x, y = y, probability = classify, seed = 1),
                args))
    }
    # probability of the class, or % cover
    to_score <- function(p) {
      if (classify) p[, "1"] else .fit_scale_to_score(p, spec)
    }
    
    fit <- fit_rf(x, y_fit)
    
    # out-of-bag, not in-sample: in-sample forest predictions are near-perfect,
    # so the threshold taken from them would be meaningless
    pred_in <- to_score(fit$predictions)
    
    # do.call() stores the evaluated arguments -- the full training matrix -- in
    # fit$call, and fit$predictions duplicates pred_in; predict() needs neither
    fit$call <- NULL
    fit$predictions <- NULL
    
    # out-of-fold predictions from the same environmental folds the lasso uses,
    # so the two engines are scored the same way. Costs one forest per fold.
    pred_oof <- rep(NA_real_, length(y_fit))
    for (f in unique(foldid)) {
      i <- foldid == f
      # importance (if in the spec) only on the global forest: it isn't needed
      # for predictions, and permutation importance is slow
      fold_args <- utils::modifyList(as.list(spec$ranger), list(importance = "none"))
      pred_oof[i] <- to_score(
        predict(fit_rf(x[!i, , drop = FALSE], y_fit[!i], fold_args),
                data = x[i, , drop = FALSE])$predictions)
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
    
  } else if (spec$engine == "gam") {
    
    # one smooth per standardized predictor, plus (spec$gam$interactions) one
    # per product of a pair, which mgcv computes from the predictor columns.
    # REML sets each term's smoothing; CV on the environmental folds picks one
    # multiplier on all of them, stored as lambda (see gam_helpers.R)
    gam_dat <- as.data.frame(x)
    gam_dat$.y <- y_fit  # the class, cover / 100, or transformed cover
    form <- gam_formula(".y", colnames(x), k = spec$gam$k, bs = spec$gam$bs,
                        interactions = spec$gam$interactions)
    family <- get(spec$family, mode = "function")()  # e.g. binomial()
    
    cv <- cv_gam_mult(form, data = gam_dat, family = family, foldid = foldid,
                      mults = 10^spec$gam$log10_mult,
                      metric = spec$gam$metric,
                      rule = spec$cv$select_rule)
    lambda <- cv$mult
    
    # global model: REML on all rows, then the selected multiplier
    sp_reml <- fit_gam(form, gam_dat, family)$sp
    gam_fit <- fit_gam(form, gam_dat, family, sp_reml = sp_reml, mult = lambda)
    
    # cv: mean out-of-fold score per multiplier (and its SE), in the columns
    # the grpreg report also uses
    fit <- list(gam = strip_gam(gam_fit),
                sp_reml = sp_reml,
                edf = gam_edf(gam_fit),
                cv = tibble(lambda = cv$summary$lambda,
                            cvm = cv$summary$score,
                            cvsd = cv$summary$score_se))
    rm(gam_fit, gam_dat)
    
    pred_in <- predict_score(list(fit = fit, config = list(spec = spec)), x)
    
    # out-of-fold predictions at the selected multiplier, from the same folds
    # that chose it (so mildly optimistic), on the scale predict_score() uses
    pred_oof <- .fit_scale_to_score(cv$oof, spec)
    rm(cv)
    
  } else {
    stop("unknown engine: ", spec$engine)
  }
  
  list(fit = fit, lambda = lambda, pred_in = pred_in, pred_oof = pred_oof)
}
