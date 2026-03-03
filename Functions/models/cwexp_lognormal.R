# cwexp_lognormal.R
# Cover-weighted exponential (cwexp) lognormal model
# implemented directly in R (likely slow), but for testing/development
# log(y_i) ~ Normal(log(mu_i), sigma^2)
# mu_i = sum_g C_ig * exp(alpha_g + X_i %*% beta_g)
# where y_i is total biomass, C_ig is cover at the ith pixel of the gth 
# plant functional group, and X is a matrix of climate predictors
# and the beta's are the learned coefficients


# ---- helpers ----

#' Prepare data for cwexp model
#'
#' Builds outcome vector, predictor matrix, and cover matrix from a formula
#' and data frame. Performs basic input checks for the lognormal model.
#'
#' @param data A data.frame containing outcome, predictors, and cover columns.
#' @param formula A model formula of the form y ~ x1 + x2 + ...
#' @param cover_cols Character vector of column names corresponding to cover variables.
#'
#' @return A list containing:
#' \itemize{
#'   \item y: numeric outcome vector
#'   \item X: predictor matrix (no intercept column)
#'   \item C: cover matrix
#'   \item n: number of observations
#'   \item G: number of cover groups
#'   \item P: number of predictors
#' }
#' @export
cwexp_prepare <- function(data, formula, cover_cols) {
  stopifnot(is.data.frame(data))
  stopifnot(inherits(formula, "formula"))
  stopifnot(all(cover_cols %in% names(data)))
  
  mf <- model.frame(formula, data = data, na.action = na.fail)
  y <- model.response(mf) # returns first column of mf
  if (!is.numeric(y)) stop("Outcome (LHS of formula) must be numeric.")
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("Outcome must be finite and strictly > 0 for the lognormal model.")
  }
  
  X <- model.matrix(formula, data = mf)
  if (ncol(X) >= 1 && colnames(X)[1] == "(Intercept)") {
    X <- X[, -1, drop = FALSE]
  }
  storage.mode(X) <- "double"
  
  C <- as.matrix(data[rownames(mf), cover_cols, drop = FALSE])
  storage.mode(C) <- "double"
  if (any(!is.finite(C))) stop("Cover matrix has non-finite values.")
  if (any(C < 0)) stop("Cover values must be >= 0.")
  
  list(
    y = as.numeric(y),
    X = X,
    C = C,
    n = length(y),
    G = ncol(C),
    P = ncol(X),
    cover_cols = cover_cols,
    x_cols = colnames(X),
    terms = terms(formula),
    outcome = all.vars(formula)[1]
  )
}

#' Pack cwexp parameters into a vector
#'
#' Combines alpha, B, and log_sigma into a single numeric vector
#' suitable for optimization.
#'
#' @param alpha Numeric vector of length G (group intercepts).
#' @param B Numeric matrix of dimension P x G (coefficients).
#' @param log_sigma Log of the error standard deviation.
#'
#' @return Numeric parameter vector.
cwexp_pack <- function(alpha, B, log_sigma) {
  c(alpha, as.vector(B), log_sigma)
}

#' Unpack cwexp parameter vector
#'
#' Converts a parameter vector back into alpha, B, and log_sigma.
#'
#' @param par Numeric parameter vector.
#' @param G Number of cover groups.
#' @param P Number of predictors.
#'
#' @return A list with components alpha, B, and log_sigma.
#' @export
cwexp_unpack <- function(par, G, P) {
  stopifnot(length(par) == G + G * P + 1)
  alpha <- par[1:G]
  B <- matrix(par[(G + 1):(G + G * P)], nrow = P, ncol = G)
  log_sigma <- par[G + G * P + 1]
  list(alpha = alpha, B = B, log_sigma = log_sigma)
}

#' Compute cwexp mean parameter
#'
#' Computes mu_i = sum_g C_ig * exp(alpha_g + X_i %*% beta_g).
#'
#' @param alpha Numeric vector of group intercepts (length G).
#' @param B Numeric matrix of coefficients (P x G).
#' @param X Predictor matrix (N x P).
#' @param C Cover matrix (N x G).
#' @param cap_eta Numeric upper bound for linear predictor (prevents overflow).
#' @param eps_mu Small positive lower bound for mu (need to be able to calculated
#' log of mu for likelihood)
#'
#' @return Numeric vector of length N containing mu values.
#' @export
cwexp_mu <- function(alpha, B, X, C, cap_eta = 50, eps_mu = 1e-12) {
  eta <- X %*% B # (N x G)
  # adding alpha[1] to column 1, alpha[2] to column 2 etc. 
  eta <- sweep(eta, 2, alpha, FUN = "+") 
  eta <- pmin(eta, cap_eta)
  mu <- rowSums(C * exp(eta))
  pmax(mu, eps_mu)
}

#' Negative log-likelihood for cwexp lognormal model
#'
#' Computes the negative log-likelihood under the model:
#' log(y_i) ~ Normal(log(mu_i), sigma^2)
#'
#' @param par Numeric parameter vector (packed via cwexp_pack()).
#' @param prep Output list from cwexp_prepare().
#' @param cap_eta Numeric upper bound for linear predictor.
#' @param eps_mu Small positive lower bound for mu.
#'
#' @return Scalar negative log-likelihood value.
#' @export
cwexp_nll_lognormal <- function(par, prep, cap_eta = 50, eps_mu = 1e-12) {
  u <- cwexp_unpack(par, G = prep$G, P = prep$P)
  
  # optimize log sigma b/ sigma can't be negative 
  sigma <- exp(u$log_sigma)
  
  mu <- cwexp_mu(
    alpha = u$alpha,
    B = u$B,
    X = prep$X,
    C = prep$C,
    cap_eta = cap_eta,
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
#' \deqn{mu_i = sum_g C_ig * exp(alpha_g + X_i %*% beta_g)}
#'
#' @param data A data.frame containing outcome, predictors, and cover columns.
#' @param formula A model formula of the form y ~ x1 + x2 + ....
#' @param cover_cols Character vector of column names corresponding to cover variables.
#' @param start Optional numeric parameter vector for optimization.
#' @param method Optimization method passed to \code{optim()}.
#' @param control Control list passed to \code{optim()}.
#' @param cap_eta Upper bound for the linear predictor to prevent overflow.
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
                                cap_eta = 50,
                                eps_mu = 1e-12) {
  
  prep <- cwexp_prepare(data = data, formula = formula, cover_cols = cover_cols)
  
  if (is.null(start)) {
    # stable-ish starts
    y_bar <- mean(prep$y) # mean response
    c_bar <- colMeans(prep$C) # mean cover
    # calculating an alpha (same for each group) gives
    # you y_bar given the mean c_bar values 
    alpha0 <- rep(log(y_bar/sum(c_bar)), times = prep$G)
    B0 <- matrix(0, nrow = prep$P, ncol = prep$G) # starting w/ intercept only model
    log_sigma0 <- log(stats::sd(log(prep$y)) + 1e-6)
    par0 <- cwexp_pack(alpha0, B0, log_sigma0)
  } else {
    par0 <- start
  }
  
  obj <- function(par) cwexp_nll_lognormal(par, prep = prep, cap_eta = cap_eta, 
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


#' Predict from cwexp model
#'
#' Generates predictions from a fitted cwexp model.
#'
#' @param object A fitted \code{"cwexp_fit"} object.
#' @param new_data A data.frame containing predictor and cover columns.
#' @param type Prediction scale: \code{"mu"} (default) returns mu,
#'   \code{"log_mu"} returns log(mu).
#' @param ... Unused.
#'
#' @return Numeric vector of predictions.
#'
#' @export
predict.cwexp_fit <- function(object, new_data, type = c("mu", "log_mu"), ...) {
  type <- match.arg(type)
  prep <- cwexp_prepare(data = new_data, formula = object$spec$formula, 
                        cover_cols = object$spec$cover_cols)
  
  mu <- cwexp_mu(
    alpha = object$par$alpha,
    B = object$par$B,
    X = prep$X,
    C = prep$C
  )
  
  if (type == "mu") return(mu)
  log(mu)
}

# optional convenience
cwexp_coef <- function(fit) {
  stopifnot(inherits(fit, "cwexp_fit"))
  list(alpha = fit$par$alpha, B = fit$par$B, sigma = fit$par$sigma)
}


# ------------------------------------------------------------
# Dummy data + quick tests for cwexp_* functions
# ------------------------------------------------------------



# ---- 1) make a synthetic dataset consistent with the model ----
cwexp_make_dummy_data <- function(n = 500,
                                  cover_cols = paste0("cov_g", 1:6),
                                  x_cols = c("tmean", "ppt", "vpd", "sand"),
                                  sigma = 0.6) {
  G <- length(cover_cols)
  P <- length(x_cols)
  
  # predictors
  Xdat <- as.data.frame(replicate(P, rnorm(n)))
  names(Xdat) <- x_cols
  
  # covers: nonnegative, with many small values
  C <- matrix(rexp(n * G, rate = 2), nrow = n, ncol = G)
  colnames(C) <- cover_cols
  
  # "true" parameters (so you can check recovery)
  alpha_true <- rnorm(G, mean = 1, sd = 0.4)
  B_true <- matrix(rnorm(P * G, mean = 0, sd = 0.25), nrow = P, ncol = G)
  rownames(B_true) <- x_cols
  colnames(B_true) <- cover_cols
  
  # model mean
  X <- as.matrix(Xdat)
  eta <- X %*% B_true
  eta <- sweep(eta, 2, alpha_true, "+")
  mu <- rowSums(C * exp(eta))
  
  # lognormal response
  y <- exp(rnorm(n, mean = log(mu), sd = sigma))
  
  dat <- cbind.data.frame(totalBio = y, Xdat, as.data.frame(C))
  
  list(
    data = dat,
    truth = list(alpha = alpha_true, B = B_true, sigma = sigma),
    cover_cols = cover_cols,
    x_cols = x_cols
  )
}

if(FALSE) {
set.seed(1)
  # consider putting some of this in a testthat type text framework
  dummy <- cwexp_make_dummy_data(n = 10000)
  dat <- dummy$data
  data <- dat
  cover_cols <- dummy$cover_cols
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
  dummy$truth$sigma
  fit$par$sigma
  dummy$truth
  
  plot(fit$par$alpha ~ dummy$truth$alpha)
  plot(as.vector(fit$par$B) ~ as.vector(dummy$truth$B)) 
  abline(0, 1)
  # Compare fitted vs true mu (on training data)
  prep_fit <- cwexp_prepare(dat, formula, cover_cols)
  mu_true <- cwexp_mu(
    alpha = dummy$truth$alpha,
    B = dummy$truth$B,
    X = prep_fit$X,
    C = prep_fit$C
  )
  cor(mu_true, mu_hat)  # should be fairly high in decent conditions
  plot(mu_true ~ mu_hat)
  plot(dummy$data$totalBio ~ mu_true)
  abline(0, 1)
}