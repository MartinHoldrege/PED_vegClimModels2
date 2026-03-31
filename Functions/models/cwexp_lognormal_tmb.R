# cwexp_tmb_wrappers.R
# R wrappers to fit the cwexp *softplus* lognormal model using TMB as the engine.
#
# Model:
#   eta_ig = alpha_g + X_i %*% beta_g
#   mu_i   = sum_g C_ig * log(1 + exp(eta_ig))      # softplus
#   log(y_i) ~ Normal(log(mu_i), sigma^2)


# helpers -------------------


# stable-ish inverse softplus for start values:
# softplus(a) = log(1 + exp(a)) = t  -> a = log(expm1(t))
inv_softplus <- function(t) log(expm1(t))

# (optional) stable softplus in R for prediction checks
softplus_R <- function(x) log1p(exp(x))  # OK if x not enormous; use ifelse() variant if needed



# TMB compilation / loading -----------------


#' Compile and load a cwexp TMB model (only recompile if needed)
#'
#' @param cpp_file Path to .cpp file (e.g., "src/cwexp_lognormal_tmb.cpp")
#' @param quiet Logical; suppress compile output where possible.
#' @param force Logical; force recompilation even if up-to-date.
#' @return DLL base name (character), to pass as DLL to MakeADFun().
cwexp_tmb_compile <- function(cpp_file = "src/cwexp_lognormal_tmb.cpp",
                              quiet = FALSE,
                              force = FALSE) {
  stopifnot(file.exists(cpp_file))
  
  cpp_base <- tools::file_path_sans_ext(cpp_file)
  dll_base <- tools::file_path_sans_ext(basename(cpp_file))
  
  dynlib_path <- TMB::dynlib(cpp_base)
  
  cpp_mtime <- file.info(cpp_file)$mtime
  dylib_exists <- file.exists(dynlib_path)
  dylib_mtime <- if (dylib_exists) file.info(dynlib_path)$mtime else as.POSIXct(NA)
  
  needs_compile <- force || !dylib_exists || is.na(dylib_mtime) || (dylib_mtime < cpp_mtime)
  
  if (needs_compile) {
    if (!quiet) message("Compiling: ", cpp_file)
    TMB::compile(cpp_file)
  } else {
    if (!quiet) message("Up-to-date (skipping compile): ", dynlib_path)
  }
  
  if (!quiet) message("Loading: ", dynlib_path)
  dyn.load(dynlib_path)
  
  dll_base
}


# main fit wrapper (TMB engine) -----


#' Fit cwexp softplus-lognormal model with TMB
#'
#' Uses TMB::MakeADFun() + optimizer (nlminb) to estimate
#' alpha (G), B (P x G), and sigma (via log_sigma).
#'
#' @param data data.frame with outcome, predictors, and cover columns
#' @param formula y ~ x1 + x2 + ...
#' @param cover_cols character vector of cover column names (length G)
#' @param dll name of loaded DLL (base name only, e.g. "cwexp_lognormal_tmb")
#' @param start optional list of starting values: list(alpha=..., B=..., log_sigma=...)
#' @param eps_mu small positive value to keep log(mu) defined (passed to C++)
#' @param control list of control parameters for optimizer
#' @param include_report logical; return report objects from cpp
#'
#' @return object of class "cwexp_tmb_fit"
cwexp_fit_tmb <- function(data,
                          formula,
                          cover_cols,
                          dll = 'cwexp_lognormal_tmb',
                          start = NULL,
                          fixed_alpha = NULL,
                          eps_mu = 1e-12,
                          control = list(iter.max = 200, eval.max = 200),
                          include_gradient = TRUE,
                          penalty = c("none", "elastic_net"),
                          en_alpha = 0.5,
                          lambda = 0,
                          l1_eps = 1e-8,
                          include_report = TRUE) {
  
  penalty <- match.arg(penalty)
  
  # prep using your existing helper (formula-based)
  prep <- cwexp_prepare(data = data, formula = formula, cover_cols = cover_cols)
  
  # ----- starting values (softplus version) -----
  inv_softplus <- function(t) log(expm1(t))  # softplus(a)=t -> a=log(expm1(t))
  
  if (is.null(start)) {
    parameters <- cwexp_start_params(prep)
  } else {
    stopifnot(is.list(start))
    stopifnot(all(c("alpha", "B", "log_sigma") %in% names(start)))
    parameters <- start
  }
  
  # ----- base data for C++ -----
  data_tmb <- list(
    y = prep$y,
    X = prep$X,
    C = prep$C,
    eps_mu  = as.numeric(eps_mu)
  )
  
  # add elastic net inputs ONLY if requested
  if (penalty == "elastic_net") {
    if (!is.finite(en_alpha) || en_alpha < 0 || en_alpha > 1) {
      stop("en_alpha must be in [0, 1].")
    }
    if (!is.finite(lambda) || lambda < 0) {
      stop("lambda must be >= 0.")
    }
    if (!is.finite(l1_eps) || l1_eps <= 0) {
      stop("l1_eps must be > 0.")
    }
    
    data_tmb$en_alpha <- as.numeric(en_alpha)
    data_tmb$lambda   <- as.numeric(lambda)
    data_tmb$l1_eps   <- as.numeric(l1_eps)
  } 
  
  # dealing with fixed alphas (when previously estimated on a different/subset 
  # of the data)
  # ----- fix selected alphas if requested -----
  map <- list()
  if (!is.null(fixed_alpha)) {
    stopifnot(is.numeric(fixed_alpha),
              !is.null(names(fixed_alpha)),
              all(names(fixed_alpha) %in% cover_cols))
    
    fix_idx <- match(names(fixed_alpha), cover_cols)
    parameters$alpha[fix_idx] <- fixed_alpha
    
    map_vec <- rep(NA_integer_, length(cover_cols))
    free_idx <- setdiff(seq_along(cover_cols), fix_idx)
    map_vec[free_idx] <- seq_along(free_idx)
    map$alpha <- factor(map_vec)
  }
  
  
  # ----- build AD object -----
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    map = map, # (empty list when nothing is fixed)
    DLL = dll,
    silent = TRUE
  )
  
  # ----- optimize with nlminb -----
  
  gradient <- if(include_gradient) obj$gr  else NULL
  
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = gradient,
    control = control
  )
  
  # ----- unpack estimates -----
  est <- obj$env$parList(opt$par)
  est$sigma <- exp(est$log_sigma)
  
  # in this model actually optimizing for
  # sqrt(sigma)--this is a bit of a hack
  # to keep estimates >=0
  if(dll == 'cwexp_lognormal_tmb3') {
    est$alpha <- est$alpha^2
  }
  
  
  names(est$alpha) <- cover_cols
  rownames(est$B) <- prep$x_cols
  colnames(est$B) <- cover_cols
  
  # (optional) capture REPORT() outputs at the optimum
  rep <- NULL
  if (isTRUE(include_report)) {
    # forces the TMB object to update its internal parameter state to the optimized parameters
    obj$fn(opt$par)
    
    rep <- obj$report()
  }
  

  out <- list(
    call = match.call(),
    spec = list(
      formula = formula,
      cover_cols = cover_cols,
      dll = dll,
      penalty = penalty,
      en_alpha = if (penalty == "elastic_net") en_alpha else NULL,
      lambda   = if (penalty == "elastic_net") lambda else NULL,
      l1_eps   = if (penalty == "elastic_net") l1_eps else NULL,
      fixed_alpha = fixed_alpha
    ),
    prep = list(x_cols = prep$x_cols, outcome = prep$outcome, terms = prep$terms),
    par  = list(alpha = est$alpha, B = est$B, sigma = est$sigma),
    tmb  = list(obj = obj, opt = opt),
    report = rep
    )
  class(out) <- "cwexp_tmb_fit"
  out
}


#' Compute lambda_max_l1 for the cwexp TMB model
#'
#' Calculates the L1 regularization scale (`lambda_max_l1`) for the cover-weighted
#' exponential log-normal TMB model.It provides the upper bound of a
#' lambda path for LASSO-style regularization.
#'
#' The function fits the model with all slopes fixed at zero. I can't 
#' can't confirm the math, but mostly just needs to provide roughly
#' the right answer. Ie this estimates that lambda value that
#' if made any smaller the B estimates would start being non-zero (note
#' here 0 is approximate--B's won't fully go to zero)
#'
#' @param data A data frame containing the response and predictor variables.
#' @param formula A model formula describing the predictors used in the
#'   exponential component.
#' @param cover_cols Character vector giving the column names of cover
#'   variables used in the cover-weighted mixture.
#' @param dll_en Character name of the compiled TMB dynamic library for the
#'   elastic-net model.
#' @param eps_mu Small numeric value added inside the softplus transformation
#'   to avoid numerical issues.
#' @param l1_eps Small smoothing constant used in the smoothed L1 penalty
#'   approximation.
#' @param control List of control parameters passed to `nlminb()` for the
#'   nuisance-parameter optimization.
#'
#' @return A list with elements:
#' \describe{
#'   \item{lambda_max_l1}{Numeric scalar giving the L1 lambda value at which
#'   slopes begin to move away from zero.}
#'   \item{grad_B}{Matrix of gradients of the mean negative log-likelihood
#'   with respect to each slope parameter evaluated at `B = 0`.}
#'   \item{opt_fixB}{`nlminb` optimization result from fitting the model with
#'   slopes fixed at zero.}
#' }
#'
#' @details
#' The returned `lambda_max_l1` is typically used as the starting point for a
#' decreasing lambda path when fitting LASSO or elastic-net regularized models.
#'
#' @export
cwexp_lambda_max_l1_tmb <- function(data,
                                    formula,
                                    cover_cols,
                                    dll_en,              # elastic-net DLL
                                    eps_mu = 1e-12,
                                    l1_eps = 1e-8,
                                    control = list(iter.max = 200, eval.max = 200)) {
  
  if (!is.finite(l1_eps) || l1_eps <= 0) {
    stop("l1_eps must be > 0.")
  }
  
  prep <- cwexp_prepare(data = data, formula = formula, cover_cols = cover_cols)
  par0 <- cwexp_start_params(prep)
  par0$B[,] <- 0 # make sure setting to 0
  
  # Dimensions
  P <- prep$P
  G <- prep$G
  
  # For lambda_max_l1 we set en_alpha = 1 (pure L1), and lambda = 0 (penalty OFF)
  # so we can compute the gradient of mean NLL at B = 0.
  data_tmb <- list(
    y = prep$y,
    X = prep$X,
    C = prep$C,
    eps_mu = as.numeric(eps_mu),
    en_alpha = 1.0,
    lambda   = 0.0,
    l1_eps   = as.numeric(l1_eps)
  )
  
  # Step 1: optimize nuisance params with B fixed at 0
  map_fixB <- list(B = factor(matrix(NA, nrow = P, ncol = G)))
  
  obj_fixB <- TMB::MakeADFun(
    data = data_tmb,
    parameters = par0,
    map = map_fixB,# this makes TMB fix (not optimize) the betas. 
    # they're fixed at they're starting (0) value
    DLL = dll_en,
    silent = TRUE
  )
  
  opt_fixB <- nlminb(
    start     = obj_fixB$par,
    objective = obj_fixB$fn,
    gradient  = obj_fixB$gr,
    control   = control
  )
  
  if (!is.null(opt_fixB$convergence) && opt_fixB$convergence != 0) {
    warning(
      "nlminb did not report convergence when fitting with B fixed at 0. ",
      "lambda_max_l1 may be unreliable. convergence=", opt_fixB$convergence
    )
  }
  
  est_fixB <- obj_fixB$env$parList(opt_fixB$par)
  
  # Step 2: rebuild objective with B free; evaluate gradient at (B=0, nuisance params optimized)
  par_at0 <- list(
    alpha = est_fixB$alpha,
    B = matrix(0, nrow = P, ncol = G),
    log_sigma = est_fixB$log_sigma
  )
  
  obj_grad <- TMB::MakeADFun(
    data = data_tmb,
    parameters = par_at0, 
    DLL = dll_en,
    silent = TRUE
  )
  
  # this returns the partial derivative of the objective function with 
  # respect to each parameter, evaluated at the current parameter values
  # (i.e. 0 for the B's). So rougly this tells you how much the thing
  # negative log likelihood would change if you increased the given Beta value
  g <- obj_grad$gr(obj_grad$par)
  
  grad_B <- cwexp_unpack(g, G, P)$B# the slopes at B = 0

  lambda_max_l1 <- max(abs(as.numeric(grad_B)))
  
  list(
    lambda_max_l1 = as.numeric(lambda_max_l1),
    grad_B = grad_B,
    opt_fixB = opt_fixB
  )
}


# example usage ---

if (FALSE) {
  
  # create fake data
  dummy <- cwexp_make_dummy_data(n = 1000)
  dat <- dummy$data
  cover_cols <- dummy$spec$cover_cols
  
  
  # 1) compile + load
  dll <- cwexp_tmb_compile("src/cwexp_lognormal_tmb.cpp")
  
  # 2) fit
  formula = totalBio ~ tmean + ppt + vpd + sand
  fit <- cwexp_fit_tmb(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    dll = dll,
    control = list(iter.max = 200, eval.max = 200)
  )
  
  fit_intercept_only <- cwexp_fit_tmb(
    data = dat,
    formula = totalBio ~ 1,
    cover_cols = cover_cols,
    dll = dll,
    control = list(iter.max = 200, eval.max = 200)
  )
  
  plot(fit$par$alpha, fit_intercept_only$par$alpha)
  fit$par$sigma
  
  # 3) predict
  mu_hat <- predict(fit, dat, type = "mu")
  mu_truth <- predict(dummy,type = 'mu')
  plot(mu_truth, mu_hat)
  abline(0, 1)
  
  # fit w/ penalization
  dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp")
  
  fit_en <- cwexp_fit_tmb(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en,
    control = list(iter.max = 200, eval.max = 200),
    penalty = "elastic_net",
    en_alpha = 0.5,
    lambda = 0.005,
    l1_eps = 1e-8
  )

  # note that the mu hat predictions can
  # be better for elastic net even  within the
  # training dataset
  mu_hat_en <- predict(fit_en, dat, type = "mu")
  plot(mu_truth, mu_hat_en) 
  abline(0, 1)
  plot(as.numeric(fit_en$par$B) ~ as.numeric(fit$par$B))
  abline(0, 1)
  

  # finding max lambda ------------------------------------------------------

  max_obj <- cwexp_lambda_max_l1_tmb(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    dll_en = dll_en
  )
  
  en_alpha = 0.5
  lambda_max_en <-  max_obj$lambda_max_l1/en_alpha
  
  # should result in > 0 B's
  fit_at_max <- cwexp_fit_tmb(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en,
    penalty = "elastic_net",
    en_alpha = 0.8*en_alpha,
    lambda = lambda_max_en,
  )
  
  max(abs(fit_at_max$par$B))
  
  # should result in near 0 B's
  fit_at_max <- cwexp_fit_tmb(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en,
    penalty = "elastic_net",
    en_alpha = 1.2*en_alpha,
    lambda = lambda_max_en,
  )
  
  max(abs(fit_at_max$par$B))
}

