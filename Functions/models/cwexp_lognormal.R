# cwexp_lognormal.R
# Cover-weighted exponential (cwexp) lognormal model
# implemented directly in R (likely slow), but for testing/development
# log(y_i) ~ Normal(log(mu_i), sigma^2)
# mu_i = sum_g C_ig * log(1 + exp(alpha_g + X_i %*% beta_g))
# where y_i is total biomass, C_ig is cover at the ith pixel of the gth 
# plant functional group, and X is a matrix of climate predictors
# and the beta's are the learned coefficients



#' Negative log-likelihood for cwexp lognormal model
#'
#' Computes the negative log-likelihood under the model:
#' log(y_i) ~ Normal(log(mu_i), sigma^2)
#'
#' @param par Numeric parameter vector (packed via cwexp_pack()).
#' @param prep Output list from cwexp_prepare().
#' @param eps_mu Small positive lower bound for mu.
#'
#' @return Scalar negative log-likelihood value.
#' @export
cwexp_nll_lognormal <- function(par, prep, eps_mu = 1e-12) {
  u <- cwexp_unpack(par, G = prep$G, P = prep$P)
  
  # optimize log sigma b/ sigma can't be negative 
  sigma <- exp(u$log_sigma)
  
  mu <- cwexp_mu(
    alpha = u$alpha,
    B = u$B,
    X = prep$X,
    C = prep$C,
    eps_mu = eps_mu
  )
  
  ll <- stats::dnorm(log(prep$y), mean = log(mu), sd = sigma, log = TRUE)
  nll <- -sum(ll)
  if (!is.finite(nll)) nll <- .Machine$double.xmax
  nll
}

# ---- public API ----

#' Fit cwexp lognormal model
#'
#' Fits the cover-weighted exponential (cwexp) model under a lognormal
#' likelihood using maximum likelihood estimation via \code{optim()}.
#'
#' Model:
#' \deqn{log(y_i) ~ Normal(log(mu_i), sigma^2)}
#' where
#' \deqn{mu_i = sum_g C_ig *log(1 +  exp(alpha_g + X_i %*% beta_g))}
#'
#' @param data A data.frame containing outcome, predictors, and cover columns.
#' @param formula A model formula of the form y ~ x1 + x2 + ....
#' @param cover_cols Character vector of column names corresponding to cover variables.
#' @param start Optional numeric parameter vector for optimization.
#' @param method Optimization method passed to \code{optim()}.
#' @param control Control list passed to \code{optim()}.
#' @param eps_mu Small positive lower bound for mu.
#'
#' @return An object of class \code{"cwexp_fit"} containing estimated
#' parameters and optimization results.
cwexp_fit_lognormal <- function(data,
                                formula,
                                cover_cols,
                                start = NULL,
                                method = "L-BFGS-B",
                                control = list(maxit = 200),
                                eps_mu = 1e-12) {
  
  prep <- cwexp_prepare(data = data, formula = formula, cover_cols = cover_cols)
  
  if (is.null(start)) {
    p <- cwexp_start_params(prep)
    par0 <- cwexp_pack(alpha = p$alpha, B = p$B, log_sigma = p$log_sigma)
  } else {
    par0 <- start
  }
  
  obj <- function(par) cwexp_nll_lognormal(par, prep = prep, 
                                           eps_mu = eps_mu)
  
  opt <- stats::optim(par = par0, fn = obj, method = method, control = control)
  u <- cwexp_unpack(opt$par, G = prep$G, P = prep$P)
  
  out <- list(
    call = match.call(),
    spec = list(formula = formula, cover_cols = cover_cols),
    prep = list(x_cols = prep$x_cols, outcome = prep$outcome, terms = prep$terms),
    par = list(alpha = u$alpha, B = u$B, sigma = exp(u$log_sigma)),
    optim = opt
  )
  class(out) <- "cwexp_fit"
  out
}




# ------------------------------------------------------------
# Dummy data + quick tests for cwexp_* functions
# ------------------------------------------------------------


if(FALSE) {
set.seed(1)
  source('Functions/data/simulate_data.R')
  source('Functions/models/cwexp_helpers.R')
  # consider putting some of this in a testthat type text framework
  dummy <- cwexp_make_dummy_data(n = 1000)
  dat <- dummy$data
  data <- dat
  cover_cols <- dummy$spec$cover_cols
  formula <- totalBio ~ tmean + ppt + vpd + sand
  
  # ---- 2) test cwexp_prepare ----
  prep <- cwexp_prepare(dat, formula = formula, cover_cols = cover_cols)
  
  str(prep)
  stopifnot(length(prep$y) == nrow(dat))
  stopifnot(nrow(prep$X) == nrow(dat))
  stopifnot(nrow(prep$C) == nrow(dat))
  stopifnot(prep$G == length(cover_cols))
  
  # ---- 3) test pack/unpack roundtrip ----
  G <- prep$G
  P <- prep$P
  
  alpha0 <- rnorm(G)
  B0 <- matrix(rnorm(P * G), nrow = P, ncol = G)
  log_sigma0 <- log(0.7)
  
  par0 <- cwexp_pack(alpha0, B0, log_sigma0)
  u0 <- cwexp_unpack(par0, G = G, P = P)
  
  stopifnot(isTRUE(all.equal(alpha0, u0$alpha)))
  stopifnot(isTRUE(all.equal(B0, u0$B)))
  stopifnot(isTRUE(all.equal(log_sigma0, u0$log_sigma)))
  
  # ---- 4) test cwexp_mu + nll ----
  mu0 <- cwexp_mu(alpha = alpha0, B = B0, X = prep$X, C = prep$C)
  stopifnot(length(mu0) == nrow(dat))
  stopifnot(all(is.finite(mu0)))
  stopifnot(all(mu0 > 0))
  
  nll0 <- cwexp_nll_lognormal(par0, prep)
  stopifnot(is.finite(nll0))
  print(nll0)
  
  # ---- 5) test cwexp_fit_lognormal end-to-end ----
  fit <- cwexp_fit_lognormal(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    control = list(maxit = 200)
  )
  
  print(fit$optim$convergence)  # 0 is good
  print(fit$par$sigma)
  
  # ---- 6) test predict() method ----
  mu_hat <- predict(fit, dat, type = "mu")
  log_mu_hat <- predict(fit, dat, type = "log_mu")
  
  stopifnot(length(mu_hat) == nrow(dat))
  stopifnot(length(log_mu_hat) == nrow(dat))
  stopifnot(all(mu_hat > 0))
  stopifnot(isTRUE(all.equal(log_mu_hat, log(mu_hat))))
  
  # ---- 7) test cwexp_coef ----
  coefs <- cwexp_coef(fit)
  str(coefs)
  
  # ---- 8) OPTIONAL: quick sanity check vs truth ----
  # NOTE: you shouldn't expect close recovery with small n / weak signal,
  # but signs/magnitudes often look reasonable.
  dummy$par$sigma
  fit$par$sigma
  dummy$par
  
  plot(fit$par$alpha ~ dummy$par$alpha) 
  plot(as.vector(fit$par$B) ~ as.vector(dummy$par$B)) 
  abline(0, 1)
  # Compare fitted vs true mu (on training data)
  prep_fit <- cwexp_prepare(dat, formula, cover_cols)
  mu_true <- predict(dummy) 
  cor(mu_true, mu_hat)  # should be fairly high in decent conditions
  plot(mu_true ~ mu_hat)
  plot(dummy$data$totalBio ~ mu_true)
  abline(0, 1)
}