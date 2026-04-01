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
  mu <- rowSums(C * ifelse(eta > 20, eta, log1p(exp(eta)))) # so doesn't blow up if eta large
  pmax(mu, eps_mu)
}

#' Predict from cwexp model
#'
#' Generates predictions from a fitted cwexp model.
#'
#' @param object A fitted \code{"cwexp_fit"} object.
#' @param newdata A data.frame containing predictor and cover columns.
#' @param type Prediction scale: \code{"mu"} (default) returns mu,
#'   \code{"log_mu"} returns log(mu).
#' @param ... Unused.
#'
#' @return Numeric vector of predictions.
#'
#' @export
predict.cwexp_fit <- function(object, newdata, type = c("mu", "log_mu"), ...) {
  type <- match.arg(type)
  prep <- cwexp_prepare(data = newdata, formula = object$spec$formula, 
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
predict.cwexp_dummy <-  function(object, newdata = NULL, 
                                 type = c("mu", "log_mu"), ...) {

  data = if(is.null(newdata)) object$data else newdata
  
  predict.cwexp_fit(object = object, newdata = data, type = type, ...)
}

# optional convenience
cwexp_coef <- function(fit) {
  stopifnot(inherits(fit, "cwexp_fit") | inherits(fit, "cwexp_tmb_fit"))
  list(alpha = fit$par$alpha, B = fit$par$B, sigma = fit$par$sigma)
}




#' Compute per-group mean contributions
#'
#' Returns an N x G matrix where column g is the contribution of group g
#' to total mu.
#'
#' @param alpha Numeric vector of group intercepts (length G).
#' @param B Numeric matrix of coefficients (P x G).
#' @param X Predictor matrix (N x P).
#' @param C Cover matrix (N x G).
#' @param weighted Logical. If TRUE (default), returns C_ig * softplus(...).
#'   If FALSE, returns softplus(alpha_g + X_i * beta_g) without cover weighting.
#'
#' @return N x G matrix. Column names inherited from C.
#' @export
cwexp_mu_by_group <- function(alpha, B, X, C, weighted = TRUE) {
  eta <- X %*% B
  eta <- sweep(eta, 2, alpha, FUN = "+")
  sp <- safe_softplus(eta)  # softplus, N x G
  
  if (weighted) {
    out <- C * sp
  } else {
    out <- sp
  }
  
  if (!is.null(colnames(C))) colnames(out) <- colnames(C)
  out
}


#' Predict per-group mu from a fitted cwexp model
#'
#' Convenience wrapper that extracts parameters and data matrices from a fitted
#' model object (or `cwexp_dummy`) and calls `cwexp_mu_by_group`.
#'
#' @param object A `cwexp_fit`, `cwexp_tmb_fit`, or `cwexp_dummy` object.
#' @param newdata Optional data frame. If NULL, uses the object's data
#'   (only works for `cwexp_dummy`).
#' @param weighted Logical. Passed to `cwexp_mu_by_group`.
#'
#' @return N x G matrix of per-group mu values.
#' @export
predict_by_group <- function(object, newdata = NULL, weighted = TRUE) {
  if ((inherits(object, "cwexp_dummy") || inherits(object, "cwexp_sim"))
      && is.null(newdata)) {
    newdata <- object$data
  }
  if (is.null(newdata)) stop("newdata is required for fitted model objects.")
  
  prep <- cwexp_prepare(
    data = newdata,
    formula = object$spec$formula,
    cover_cols = object$spec$cover_cols
  )
  
  cwexp_mu_by_group(
    alpha = object$par$alpha,
    B = object$par$B,
    X = prep$X,
    C = prep$C,
    weighted = weighted
  )
}


#' Compute per-observation log-likelihood for a cwexp model
#'
#' Computes the log-likelihood contribution of each observation under the
#' lognormal model. Supports both the standard `log(y) ~ Normal(log(mu), sigma^2)`
#' and the shifted `log(y+1) ~ Normal(log(mu+1), sigma^2)` formulations.
#'
#' @param fit Fitted cwexp model object (class `cwexp_tmb_fit`).
#' @param data Data frame with predictor, cover, and response columns.
#' @param response Character; name of the observed response column.
#'   Default, determined from fit by default
#' @param shift Numeric; value added before taking log. Default 1 for the
#'   `log(y+1)` formulation. Use 0 for the standard `log(y)` formulation.
#'
#' @return A numeric vector of per-observation log-likelihood values (length N).
#' @examples
#' dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#' dummy <- cwexp_make_dummy_data(n = 500)
#' fit <- cwexp_fit_tmb(dummy$data, dummy$spec$formula,
#'                      dummy$spec$cover_cols, dll = dll_en,
#'                      penalty = "elastic_net", en_alpha = 0.5, lambda = 0)
#' ll <- cwexp_loglik(fit, dummy$data)
#' plot(predict(fit, dummy$data), ll)
#' @export
cwexp_loglik <- function(fit, data, response = NULL, shift = NULL) {
  
  if (is.null(response)) {
    stopifnot(length(fit$spec$formula) == 3)
    response <-  fit$spec$formula[2] |> as.character()
  }
  stopifnot(response %in% names(data))
  
  y <- data[[response]]
  mu <- predict(fit, data, type = "mu")
  sigma <- fit$par$sigma
  
  log_y  <- log(y + shift)
  log_mu <- log(mu + shift)
  
  dnorm(log_y, mean = log_mu, sd = sigma, log = TRUE)
}


# predict using raster input ----------------------------------------------

#' Predict total biomass from a cwexp model onto a raster
#'
#' Computes `mu = sum_g(C_g * softplus(alpha_g + X * beta_g))` using raster
#' math. The input raster must have named layers matching the predictor
#' variables and cover columns in the fitted model. Predictor layers should
#' be pre-standardized to match the training data scale.
#'
#' @param fit Fitted cwexp model object (`cwexp_tmb_fit`).
#' @param rast A `SpatRaster` (terra) with named layers including all
#'   predictor variables (from `fit$prep$x_cols`) and all cover columns
#'   (from `fit$spec$cover_cols`).
#' @param type Character; `"total"` (default) returns a single-layer
#'   SpatRaster of total predicted biomass. `"by_group"` returns a
#'   multi-layer SpatRaster with one layer per PFT (cover-weighted).
#'   `"potential"` returns a multi-layer SpatRaster of per-PFT biomass
#'   at 100% cover (softplus(eta) without cover weighting).
#'
#' @return A `SpatRaster`.
#' @examples
#' dummy <- cwexp_make_dummy_data(n = 500)
#' dll <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb2.cpp", quiet = TRUE)
#' fit <- cwexp_fit_tmb(dummy$data, dummy$spec$formula,
#'                      dummy$spec$cover_cols, dll = dll,
#'                      penalty = "elastic_net", en_alpha = 0.5, lambda = 0)
#' rast <- create_dummy_raster(fit)
#' type = "total"
#' pred <- predict_raster(fit, rast, type = type)
#' terra::plot(pred)
#' @export
predict_raster <- function(fit, rast,
                           type = c("total", "by_group", "potential")) {
  
  type <- match.arg(type)
  
  stopifnot(inherits(rast, "SpatRaster"),
            str_detect(class(fit), 'cwexp'))
  
  x_cols <- fit$prep$x_cols
  cover_cols <- fit$spec$cover_cols
  alpha <- fit$par$alpha
  B <- fit$par$B
  
  required_layers <- c(x_cols, cover_cols)
  missing <- setdiff(required_layers, names(rast))
  if (length(missing) > 0) {
    stop("Missing layers in raster: ", paste(missing, collapse = ", "))
  }
  
  G <- length(cover_cols)
  P <- length(x_cols)
  
  # compute eta_g = alpha_g + sum_p(X_p * B_pg) for each group
  # then softplus(eta_g), then optionally multiply by cover_g
  
  # helper: safe softplus for rasters
  r_softplus <- function(r) {
    terra::ifel(r > 20, r, log1p(exp(r)))
  }
  
  # compute per-group contributions
  group_layers <- vector("list", G)
  for (g in seq_len(G)) {
    # start with alpha_g
    eta_g <- terra::rast(rast[[1]])  # template
    terra::values(eta_g) <- alpha[g]
    
    # add predictor contributions
    for (p in seq_len(P)) {
      eta_g <- eta_g + rast[[x_cols[p]]] * B[p, g]
    }
    
    sp_g <- r_softplus(eta_g)
    
    if (type == "potential") {
      group_layers[[g]] <- sp_g
    } else {
      group_layers[[g]] <- rast[[cover_cols[g]]] * sp_g
    }
  }
  
  names(group_layers) <- cover_cols
  group_layers <- terra::rast(group_layers)
  if (type == "total") {
    # sum across groups
    out2 <- terra::app(group_layers, fun = 'sum')
    names(out) <- paste('predicted', fit$prep$outcome, sep = "_")
    return(out)
  }
  
  # by_group or potential: stack into multi-layer raster
  out <- group_layers
  pft_labels <- stringr::str_remove(cover_cols, "Cov$")
  names(out) <- pft_labels
  out
}
