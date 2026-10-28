# Fitting shared by the cover models: CV folds, and the fitting engines
# (glmnet, ranger, grpreg, grplasso, gam). Used by 02_fit_classification.R and
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
#'dat <- dat |> 
#'  filter(!is.na(cov_tree)) |> 
#'  slice_sample(n = 1000)
#' #spec <- cover_specs$classification$forest$m01
#' spec <- cover_specs$cover$tree_forest$m01.51
#' x <- prepare_predictors(dat, pred_vars = spec$pred_vars,
#'                         scale_df = read_scale_params(),
#'                         squares = spec$squares,
#'                         interactions = spec$interactions,
#'                         interact_log1p = spec$interact_log1p)
#' # y_fit <- as.integer(dat$cov_tree > spec$cover_threshold)
#' y_fit <- dat$cov_tree/100
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
    # spec$grpreg$hierarchy picks the rule ("strict" by default; see
    # hier_groups()); the rest of spec$grpreg goes to cv.grpreg()
    grpreg_args <- spec$grpreg
    rule <- if (is.null(grpreg_args$hierarchy)) "strict" else grpreg_args$hierarchy
    grpreg_args$hierarchy <- NULL
    groups <- hier_groups(colnames(x), rule = rule)
    
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
        grpreg_args)
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
  
    # a group is in the model if any of its copies is non-zero; every term in
    # it must then be in the model (group names aren't term names under the
    # "variable" rule, so this uses the copies, not the summed coefficients)
    stopifnot(length(b) == nrow(groups) + 1)
    active <- unique(groups$group[b[-1] != 0])
    required <- groups$term[groups$group %in% active]
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
    
  } else if (spec$engine == "grplasso") {

    # latent overlapping groups, as for grpreg, with the hierarchy rule in
    # spec$grplasso$hierarchy (see hier_groups()). Unlike grpreg, grplasso
    # takes a proportion response (quasibinomial; see .grpl_model()). Intercept
    # unpenalized (index NA). grplasso centers and blockwise orthonormalizes
    # the columns internally, and returns coefficients on the scale of x
    groups <- hier_groups(colnames(x), rule = spec$grplasso$hierarchy)
    x_grp <- cbind("(Intercept)" = 1, x[, groups$term, drop = FALSE])
    index <- c(NA, as.integer(factor(groups$group,
                                     levels = unique(groups$group))))
    model <- .grpl_model(spec$family)
    ctrl <- grplasso::grpl.control(trace = 0, save.y = FALSE)

    # log-spaced path down from the smallest lambda that keeps every group out.
    # grplasso sums the loss over rows, so its lambda scales with n; these
    # lambdas are per row (as glmnet's) and are multiplied by the number of
    # rows in each fit
    lambda_max <- grplasso::lambdamax(x_grp, y = y_fit, index = index,
                                      model = model) / nrow(x_grp)
    lambdas <- lambda_max * spec$grplasso$lambda_min_ratio^
      seq(0, 1, length.out = spec$grplasso$nlambda)

    # out-of-fold linear predictors for every lambda; folds fit in parallel
    # (see .grpl_folds())
    eta_oof <- .grpl_folds(x_grp = x_grp, y_fit = y_fit, foldid = foldid,
                           index = index, lambdas = lambdas, model = model,
                           ctrl = ctrl)

    # CV error per lambda: mean out-of-fold deviance (binomial deviance for
    # binomial and quasibinomial, squared error for gaussian), and its SE
    # across folds
    loss <- .grpl_loss(y_fit, eta_oof, spec$family)
    fold_means <- rowsum(loss, foldid) / as.vector(table(foldid))
    cvm <- colMeans(loss)
    cvsd <- apply(fold_means, 2, sd) / sqrt(nrow(fold_means))
    rm(loss)

    i_min <- which.min(cvm)
    if (i_min == length(lambdas)) {
      warning("CV error is lowest at the smallest lambda: the path may need",
              " a smaller lambda_min_ratio in spec$grplasso")
    }
    # lambdas decrease along the path, so the first within 1 SE is the
    # largest
    i_sel <- switch(spec$cv$select_rule,
                    min = i_min,
                    "1se" = min(which(cvm <= cvm[i_min] + cvsd[i_min])),
                    stop("unknown select_rule: ", spec$cv$select_rule))
    lambda <- lambdas[i_sel]

    fit_all <- grplasso::grplasso(x_grp, y = y_fit, index = index,
                                  lambda = lambda * nrow(x_grp),
                                  model = model, control = ctrl)
    b <- fit_all$coefficients[, 1]  # one per column of x_grp
    converged <- fit_all$converged
    rm(fit_all)
    if (!isTRUE(converged)) {
      warning("grplasso did not converge at the selected lambda")
    }

    # one coefficient per column of x (the sum over its copies), which is all
    # predict_score() needs
    coefs <- c(tapply(b, factor(colnames(x_grp),
                                levels = c("(Intercept)", colnames(x))), sum))

    # hierarchy check. A group is selected if any of its columns has a
    # non-zero coefficient. Every term in a selected group must then be in the
    # final equation, e.g. MAT and MAP when the MAT:MAP group is selected
    selected_groups <- unique(groups$group[b[-1] != 0])
    required_terms <- groups$term[groups$group %in% selected_groups]
    stopifnot("a selected term is missing a term it requires" =
                all(coefs[required_terms] != 0))

    # cv: mean out-of-fold deviance per lambda (and its SE), in the columns
    # the grpreg report also uses
    fit <- list(coef = coefs,
                groups = groups,
                cv = tibble(lambda = lambdas, cvm = cvm, cvsd = cvsd))

    pred_in <- predict_score(list(fit = fit, lambda = lambda,
                                  config = list(spec = spec)), x)

    # out-of-fold predictions at the selected lambda, from the same folds that
    # chose it (so mildly optimistic)
    pred_oof <- .grpl_score(eta_oof[, i_sel], spec)
    rm(eta_oof, x_grp)

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


#' grplasso model for a spec's family
#'
#' grplasso's `LogReg()` checks that the response is 0/1, but its
#' likelihood, gradient and Hessian are the binomial ones, which also fit a
#' proportion (the quasi-likelihood glmnet fits with a two-column response).
#' For quasibinomial, the same model with the check relaxed to 0 <= y <= 1.
#'
#' @param family "binomial", "quasibinomial" or "gaussian".
#' @return A `grpl.model` object.
#' @examples
#' family <- "quasibinomial"
#' .grpl_model(family = family)
.grpl_model <- function(family) {
  switch(
    family,
    binomial = grplasso::LogReg(),
    gaussian = grplasso::LinReg(),
    quasibinomial = {
      m <- grplasso::LogReg()
      grplasso::grpl.model(invlink = m@invlink, link = m@link,
                           nloglik = m@nloglik, ngradient = m@ngradient,
                           nhessian = m@nhessian,
                           check = function(y) all(y >= 0 & y <= 1),
                           name = "Quasibinomial (logit link)",
                           comment = "Response is a proportion, 0-1")
    },
    stop("family not set up for grplasso: ", family)
  )
}


#' Out-of-fold loss for each row and lambda
#'
#' Binomial deviance (binomial and quasibinomial; for a proportion y, 0 log 0
#' is taken as 0) or squared error (gaussian), from linear predictors.
#'
#' @param y Response, as fit (0/1, proportion, or transformed cover).
#' @param eta Matrix of linear predictors, one row per element of `y`, one
#'   column per lambda.
#' @param family "binomial", "quasibinomial" or "gaussian".
#' @return Matrix the shape of `eta`.
#' @examples
#' y <- c(0, 0.2, 1)
#' eta <- cbind(c(-2, -1, 2), c(0, 0, 0))
#' family <- "quasibinomial"
#' .grpl_loss(y = y, eta = eta, family = family)
.grpl_loss <- function(y, eta, family) {
  if (family == "gaussian") return((y - eta)^2)
  # y repeated in every column (one per lambda), so it lines up with eta
  y_mat <- matrix(y, nrow = nrow(eta), ncol = ncol(eta))
  # binomial deviance per row, as glm() computes it (also quasibinomial's)
  stats::binomial()$dev.resids(y_mat, plogis(eta), 1)
}


# internal: grplasso linear predictor to what predict_score() returns
# (probability, or % cover)
.grpl_score <- function(eta, spec) {
  mu <- if (spec$family == "gaussian") eta else plogis(eta)
  .fit_scale_to_score(mu, spec)
}


#' Out-of-fold grplasso linear predictors, folds fit in parallel
#'
#' One grplasso path per fold, run with future.apply. If the caller has set
#' a future plan it is used; otherwise (the sequential default) one is set
#' here with `workers` workers and reset on exit: forked processes
#' (multicore) where supported, as on Linux/HPC outside RStudio, else
#' separate R sessions (multisession), each getting its own copy of `x_grp`,
#' so memory grows with the number of workers.
#'
#' @param x_grp Expanded design matrix, with an intercept column.
#' @param y_fit Response, as fit.
#' @param foldid Fold of each row.
#' @param index grplasso group index (NA for the intercept).
#' @param lambdas Per-row lambdas, multiplied by each fold's training rows.
#' @param model,ctrl grplasso model and control objects.
#' @param workers Number of workers; by default one per fold, capped at
#'   `future::availableCores()` (which follows SLURM_CPUS_PER_TASK etc.).
#' @return Matrix of linear predictors, one row per row of `x_grp`, one
#'   column per lambda.
#' @examples
#' set.seed(4402)
#' x_grp <- cbind("(Intercept)" = 1, a = rnorm(200), b = rnorm(200))
#' y_fit <- rbinom(200, 1, plogis(x_grp[, "a"]))
#' foldid <- rep(1:4, 50)
#' index <- c(NA, 1, 2)
#' lambdas <- c(0.05, 0.01, 0.001)
#' model <- grplasso::LogReg()
#' ctrl <- grplasso::grpl.control(trace = 0)
#' .grpl_folds(x_grp, y_fit, foldid, index, lambdas, model, ctrl, workers = 2)
.grpl_folds <- function(x_grp, y_fit, foldid, index, lambdas, model, ctrl,
                        workers = min(length(unique(foldid)),
                                      future::availableCores())) {
  if (inherits(future::plan(), "sequential") && workers > 1) {
    if (future::supportsMulticore()) {
      old_plan <- future::plan(future::multicore, workers = workers)
    } else {
      # explicit cluster, closed with stopCluster() on exit rather than left
      # for the garbage collector (which warns "closing unused connection")
      cl <- parallelly::makeClusterPSOCK(workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      old_plan <- future::plan(future::cluster, workers = cl)
    }
    on.exit(future::plan(old_plan), add = TRUE)
  }
  # x_grp is shipped to the workers, and can be several GB
  old_opt <- options(future.globals.maxSize = Inf)
  on.exit(options(old_opt), add = TRUE)

  folds <- unique(foldid)
  res <- future.apply::future_lapply(folds, function(f) {
    i <- foldid == f
    fit_f <- grplasso::grplasso(x_grp[!i, , drop = FALSE], y = y_fit[!i],
                                index = index, lambda = lambdas * sum(!i),
                                model = model, control = ctrl)
    list(eta = x_grp[i, , drop = FALSE] %*% fit_f$coefficients,
         converged = all(fit_f$converged))
  })

  if (!all(vapply(res, `[[`, logical(1), "converged"))) {
    warning("grplasso did not converge for some lambdas in some CV folds")
  }
  eta_oof <- matrix(NA_real_, nrow = nrow(x_grp), ncol = length(lambdas))
  for (k in seq_along(folds)) {
    eta_oof[foldid == folds[k], ] <- res[[k]]$eta
  }
  stopifnot(!anyNA(eta_oof))
  eta_oof
}
