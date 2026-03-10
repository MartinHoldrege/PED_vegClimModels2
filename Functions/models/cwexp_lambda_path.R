# functions for selecting the best lambda (penalization) for elastic net



#' Create an elastic-net lambda path
#'
#' Builds a decreasing sequence of lambda values for a given elastic-net
#' mixing parameter. The path is scaled relative to the alpha-specific
#' lambda maximum and can optionally include 0 as the final value.
#'
#' @param lambda_max_l1 Numeric scalar. L1-scale lambda maximum.
#' @param en_alpha Numeric scalar in (0, 1]. Elastic-net mixing parameter.
#' @param n_lambda Integer. Number of positive lambda values.
#' @param min_ratio Numeric scalar in (0, 1). Smallest positive lambda is
#'   `min_ratio * lambda_max`.
#' @param lambda_max_mult Numeric scalar > 0. Multiplier applied to the
#'   alpha-specific lambda maximum. should probably be
#'   >=1 in most cases
#' @param include_zero Logical. If `TRUE`, append 0 to the end of the path.
#'
#' @return Numeric vector of lambda values in decreasing order, with 0
#'   appended last if requested.
cwexp_lambda_path <- function(lambda_max_l1,
                              en_alpha,
                              n_lambda = 10,
                              min_ratio = 1e-3,
                              lambda_max_mult = 1.2,
                              include_zero = TRUE) {
  
  if (!is.finite(lambda_max_l1) || lambda_max_l1 <= 0) {
    stop("lambda_max_l1 must be a positive finite number.")
  }
  if (!is.finite(en_alpha) || en_alpha <= 0 || en_alpha > 1) {
    stop("en_alpha must be in (0, 1].")
  }
  if (!is.numeric(n_lambda) || length(n_lambda) != 1 || n_lambda < 1) {
    stop("n_lambda must be a positive integer.")
  }
  if (!is.finite(min_ratio) || min_ratio <= 0 || min_ratio >= 1) {
    stop("min_ratio must be in (0, 1).")
  }
  if (!is.finite(lambda_max_mult) || lambda_max_mult <= 0) {
    stop("lambda_max_mult must be > 0.")
  }
  
  lambda_max <- lambda_max_mult * (lambda_max_l1 / en_alpha)
  lambda_min <- min_ratio * lambda_max
  
  if (n_lambda == 1) {
    lam <- lambda_max
  } else {
    lam <- exp(seq(log(lambda_max), log(lambda_min), length.out = n_lambda))
  }
  
  if (include_zero) {
    lam <- c(lam, 0)
  }
  
  lam
}


#' Extract warm-start values from a cwexp TMB fit
#'
#' Creates a start list that can be passed back into `cwexp_fit_tmb()`.
#'
#' @param fit A `cwexp_tmb_fit` object.
#'
#' @return A named list with elements `alpha`, `B`, and `log_sigma`.
cwexp_tmb_start_from_fit <- function(fit) {
  if (!inherits(fit, "cwexp_tmb_fit")) {
    stop("fit must inherit from 'cwexp_tmb_fit'.")
  }
  
  list(
    alpha = fit$par$alpha,
    B = fit$par$B,
    log_sigma = log(fit$par$sigma)
  )
}

#' Fit cwexp TMB model across a lambda path
#'
#' Fits the elastic-net cwexp TMB model over a sequence of lambda values,
#' using the previous fit as a warm start for the next lambda.
#'
#' @param data Data frame used for fitting.
#' @param formula Model formula.
#' @param cover_cols Character vector of cover column names.
#' @param dll Character. Compiled TMB DLL name.
#' @param en_alpha Numeric scalar in [0, 1]. Elastic-net mixing parameter.
#' @param lambda_seq Numeric vector of lambda value (from cwexp_lambda_path)
#' @param start Optional named list of starting values for the first fit.
#' @param l1_eps Numeric scalar > 0. Smoothing constant for L1 penalty.
#' @param control List passed to `nlminb()` through `cwexp_fit_tmb()`.
#' @param include_report Logical. Passed to `cwexp_fit_tmb()`.
#'
#' @return A list with path settings, a summary data frame, and the full fit
#'   objects for each lambda.
#' @examples
#' # Simulate a small dataset
#' set.seed(1)
#' dummy <- cwexp_make_dummy_data(n = 500)
#'
#' dat <- dummy$data
#' cover_cols <- dummy$spec$cover_cols
#' formula <- totalBio ~ tmean + ppt + vpd + sand
#'
#' # Compile elastic-net TMB model
#' dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#'
#' # Compute lambda_max on L1 scale
#' lam_max <- cwexp_lambda_max_l1_tmb(
#'   data = dat,
#'   formula = formula,
#'   cover_cols = cover_cols,
#'   dll_en = dll_en
#' )
#'
#' # Build a short lambda path
#' lam_seq <- cwexp_lambda_path(
#'   lambda_max_l1 = lam_max$lambda_max_l1,
#'   en_alpha = 0.5,
#'   n_lambda = 5,
#'   include_zero = TRUE
#' )
#'
#' # Fit models across the lambda path
#' path_fit <- cwexp_fit_lambda_path_tmb(
#'   data = dat,
#'   formula = formula,
#'   cover_cols = cover_cols,
#'   dll = dll_en,
#'   en_alpha = 0.5,
#'   lambda_seq = lam_seq
#' )
#'
#' # True mean biomass from the simulation
#' mu_true <- predict(dummy, type = "mu")
#' # Compute correlation between true and predicted mu
#' cor_mu <- sapply(path_fit$fits, function(fit) {
#'   mu_hat <- predict(fit, dat, type = "mu")
#'   cor(mu_true, mu_hat)
#' })
#'
#' # Plot correlation vs lambda
#' plot(
#'   path_fit$summary$lambda,
#'   cor_mu,
#'   xlab = "Lambda",
#'   ylab = "cor(true mu, predicted mu)",
#'   main = "true mean recovery across lambda path",
#'   pch = 19
#' )
#' lines(path_fit$summary$lambda, cor_mu)
#' df_B <- cwexp_extract_B_path(path_fit)
#' plot(B ~ lambda, data = df_B)
#' # correlation between fitted and true B's
#' rmse_B <- unlist(lapply(path_fit$fits, function(fit) {
#'  yardstick::rmse_vec(as.numeric(fit$par$B), as.numeric(dummy$par$B))
#' }))
#' plot(
#'   path_fit$summary$lambda,
#'   rmse_B,
#'   xlab = "Lambda",
#'   ylab = "RMSE(true B, predicted B)",
#'   main = "true B recovery across lambda path",
#'   pch = 19
#' )
cwexp_fit_lambda_path_tmb <- function(data,
                                      formula,
                                      cover_cols,
                                      dll,
                                      en_alpha = 0.5,
                                      lambda_seq,
                                      start = NULL,
                                      l1_eps = 1e-8,
                                      control = list(iter.max = 200, eval.max = 200),
                                      include_report = TRUE) {
  
  if (!is.finite(en_alpha) || en_alpha < 0 || en_alpha > 1) {
    stop("en_alpha must be in [0, 1].")
  }
  if (!is.numeric(lambda_seq) || length(lambda_seq) < 1) {
    stop("lambda_seq must be a numeric vector of length >= 1.")
  }
  if (any(!is.finite(lambda_seq)) || any(lambda_seq < 0)) {
    stop("lambda_seq must contain finite values >= 0.")
  }
  
  n_fit <- length(lambda_seq)
  fits <- vector("list", n_fit)
  
  summary_list <- vector("list", n_fit)
  current_start <- start
  
  for (i in seq_along(lambda_seq)) {

    lam <- lambda_seq[i]
    
    fit_i <- cwexp_fit_tmb(
      data = data,
      formula = formula,
      cover_cols = cover_cols,
      dll = dll,
      start = current_start,
      penalty = "elastic_net",
      en_alpha = en_alpha,
      lambda = lam,
      l1_eps = l1_eps,
      control = control,
      include_report = include_report
    )
    
    fits[[i]] <- fit_i
    
    summary_list[[i]] <- data.frame(
      lambda = lam,
      en_alpha = en_alpha,
      convergence = fit_i$tmb$opt$convergence,
      objective = unname(fit_i$tmb$opt$objective),
      sigma = fit_i$par$sigma,
      sum_abs_B = sum(abs(fit_i$par$B)),
      max_abs_B = max(abs(fit_i$par$B)),
      nll = if (!is.null(fit_i$report)) fit_i$report$nll else NA_real_,
      penalty = if (!is.null(fit_i$report)) fit_i$report$penalty else NA_real_,
      obj_report = if (!is.null(fit_i$report)) fit_i$report$obj else NA_real_
    )
    
    current_start <- cwexp_tmb_start_from_fit(fit_i)
  }
  
  summary_df <- do.call(rbind, summary_list)
  
  out <- list(
    spec = list(
      formula = formula,
      cover_cols = cover_cols,
      dll = dll,
      en_alpha = en_alpha,
      lambda_seq = lambda_seq,
      l1_eps = l1_eps
    ),
    summary = summary_df,
    fits = fits
  )
  
  class(out) <- "cwexp_lambda_path_tmb_fit"
  out
}


#' Extract B coefficient path from a lambda-path fit
#'
#' Returns the slope coefficient path in a wide data-frame format with one row
#' per lambda and predictor, and one column for each group-specific coefficient.
#'
#' @param path_fit A `cwexp_lambda_path_tmb_fit` object.
#'
#' @return A data frame with columns `lambda`, `predictor`, and one column per
#'   group-specific coefficient.
cwexp_extract_B_path <- function(path_fit) {
  if (!inherits(path_fit, "cwexp_lambda_path_tmb_fit")) {
    stop("path_fit must inherit from 'cwexp_lambda_path_tmb_fit'.")
  }
  
  n_fit <- length(path_fit$fits)
  out <- vector("list", n_fit)
  
  # Try to recover predictor names from the first fit
  predictor_names <- path_fit$fits[[1]]$prep$x_cols
  if (is.null(predictor_names)) {
    predictor_names <- paste0("x", seq_len(nrow(path_fit$fits[[1]]$par$B)))
  }
  
  for (i in seq_along(path_fit$fits)) {
    fit_i <- path_fit$fits[[i]]
    B_i <- fit_i$par$B
    
    n_group <- ncol(B_i)
    group_names <- path_fit$spec$cover_cols
    
    out_i <- data.frame(
      lambda = path_fit$summary$lambda[i],
      predictor = predictor_names,
      B_i,
      row.names = NULL,
      check.names = FALSE
    )
    
    names(out_i)[-(1:2)] <- group_names
    out[[i]] <- tidyr::pivot_longer(out_i, cols = dplyr::all_of(group_names),
                                    values_to = 'B',
                                    names_to = 'group')
  }
  
  do.call(rbind, out)
}


#' Extract lightweight parameter objects from a lambda-path fit
#'
#' Returns one small parameter object per lambda, suitable for inspection or
#' later prediction.
#'
#' @param path_fit A `cwexp_lambda_path_tmb_fit` object.
#'
#' @return A list with one element per lambda. Each element contains
#'   `lambda`, `en_alpha`, and `par`.
extract_path_par_objects <- function(path_fit) {
  if (!inherits(path_fit, "cwexp_lambda_path_tmb_fit")) {
    stop("path_fit must inherit from 'cwexp_lambda_path_tmb_fit'.")
  }
  
  out <- vector("list", length(path_fit$fits))
  
  for (i in seq_along(path_fit$fits)) {
    fit_i <- path_fit$fits[[i]]
    
    out[[i]] <- list(
      lambda = path_fit$summary$lambda[i],
      en_alpha = path_fit$summary$en_alpha[i],
      par = fit_i$par
    )
  }
  
  out
}

