# helper functions for fitting cover weighted exponential lognormal
# models


# prep input data ---------------------------------------------------------


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


# parameter related funs --------------------------------------------------


#' Pick starting parameter values for cwexp  model
#'
#' Chooses a common intercept a_bar so that, with B = 0, the mean model
#' matches mean(y) at mean cover: sum_g mean(C_g) * softplus(a_bar) = mean(y).
#'
#' @param prep Output from cwexp_prepare().
#' @return A list with elements alpha (length G), B (P x G), log_sigma (scalar).
cwexp_start_params <- function(prep) {
  y_bar <- mean(prep$y)
  c_bar <- colMeans(prep$C)
  c_sum <- sum(c_bar)
  
  if (!is.finite(c_sum) || c_sum <= 0) {
    stop("sum(colMeans(C)) must be > 0 to construct default start values.")
  }
  
  # calculating an alpha (same for each group) gives
  # you y_bar given the mean c_bar values 
  a_bar <- log(expm1(y_bar/c_sum))
  B0 <- matrix(0, nrow = prep$P, ncol = prep$G)
  list(
    alpha = rep(a_bar, prep$G),
    B = B0,
    log_sigma = log(stats::sd(log(prep$y)) + 1e-6)
  )
}

# other helpers defined in cwexp_helpers.R

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

# predict on fitted object ------------------------------------------------

#' Compute cwexp mean parameter 
#' 
#' (also used in 
#' cwexp_lognormal.R for fitting in base R instead of TMB)
#'
#' Computes mu_i = sum_g C_ig *log(1 + exp(alpha_g + X_i %*% beta_g)).
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
cwexp_mu <- function(alpha, B, X, C, eps_mu = 1e-12) {
  eta <- X %*% B # (N x G)
  # adding alpha[1] to column 1, alpha[2] to column 2 etc. 
  eta <- sweep(eta, 2, alpha, FUN = "+") 
  mu <- rowSums(C * log1p(exp(eta)))
  pmax(mu, eps_mu)
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

#' @export
predict.cwexp_tmb_fit <- predict.cwexp_fit


# predicts mu from object created in cwexp_make_dummy_data
predict.cwexp_dummy <-  function(object, new_data = NULL, 
                                 type = c("mu", "log_mu"), ...) {

  data = if(is.null(new_data)) object$data else new_data
  
  predict.cwexp_fit(object = object, new_data = data, type = type, ...)
}

# optional convenience
cwexp_coef <- function(fit) {
  stopifnot(inherits(fit, "cwexp_fit") | inherits(fit, "cwexp_tmb_fit"))
  list(alpha = fit$par$alpha, B = fit$par$B, sigma = fit$par$sigma)
}





