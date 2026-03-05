# cwexp_tmb_wrappers.R
# R wrappers to fit the cwexp *softplus* lognormal model using TMB as the engine.
#
# Model:
#   eta_ig = alpha_g + X_i %*% beta_g
#   mu_i   = sum_g C_ig * log(1 + exp(eta_ig))      # softplus
#   log(y_i) ~ Normal(log(mu_i), sigma^2)

# ---------------------------
# helpers
# ---------------------------

# stable-ish inverse softplus for start values:
# softplus(a) = log(1 + exp(a)) = t  -> a = log(expm1(t))
inv_softplus <- function(t) log(expm1(t))

# (optional) stable softplus in R for prediction checks
softplus_R <- function(x) log1p(exp(x))  # OK if x not enormous; use ifelse() variant if needed


# ---------------------------
# TMB compilation / loading
# ---------------------------

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

# ---------------------------
# main fit wrapper (TMB engine)
# ---------------------------

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
#' @param do_sdreport logical; compute sdreport() (can be slow)
#'
#' @return object of class "cwexp_tmb_fit"
cwexp_fit_tmb <- function(data,
                          formula,
                          cover_cols,
                          dll = 'cwexp_lognormal_tmb',
                          start = NULL,
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
  
  # ----- build AD object -----
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    DLL = dll
  )
  
  # ----- optimize with nlminb -----
  
  gradient <- if(include_gradient) obj$gr  else NULL
  
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = control
  )
  
  # ----- unpack estimates -----
  est <- obj$env$parList(opt$par)
  est$sigma <- exp(est$log_sigma)
  
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
      l1_eps   = if (penalty == "elastic_net") l1_eps else NULL
    ),
    prep = list(x_cols = prep$x_cols, outcome = prep$outcome, terms = prep$terms),
    par  = list(alpha = est$alpha, B = est$B, sigma = est$sigma),
    tmb  = list(obj = obj, opt = opt),
    report = rep
    )
  class(out) <- "cwexp_tmb_fit"
  out
}




# ---------------------------
# example usage
# ---------------------------
if (FALSE) {
  
  # create fake data
  dummy <- cwexp_make_dummy_data(n = 1000)
  dat <- dummy$data
  cover_cols <- dummy$spec$cover_cols
  
  
  # 1) compile + load
  dll <- cwexp_tmb_compile("src/cwexp_lognormal_tmb.cpp")
  
  # 2) fit

  fit <- cwexp_fit_tmb(
    data = dat,
    formula = totalBio ~ tmean + ppt + vpd + sand,
    cover_cols = cover_cols,
    dll = dll,
    control = list(iter.max = 200, eval.max = 200)
  )
  
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
    formula = totalBio ~ tmean + ppt + vpd + sand,
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
}

