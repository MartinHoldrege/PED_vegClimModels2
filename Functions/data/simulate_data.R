# functions for simulation data
# (simulated data is for model testing)
# Started Jan 28, 2026 by Martin Holdrege

# helper: build formula from main effects + interaction tuples
make_formula <- function(pred_vars, inter = list()) {
  terms <- pred_vars
  
  # inter is list of character vectors, e.g. list(c("tmean_CLIM","precip_CLIM"))
  if (length(inter) > 0) {
    inter_terms <- purrr::map_chr(inter, ~ paste0(.x, collapse = ":"))
    terms <- c(terms, inter_terms)
  }
  
  rhs <- paste(terms, collapse = " + ")
  
  stats::as.formula(paste("~", rhs))
}


# helper: turn your coefs list into a coefficient matrix aligned to X column names
coefs_to_matrix <- function(coefs, X_colnames, pfts) {
  B <- matrix(0, nrow = length(X_colnames), ncol = length(pfts),
              dimnames = list(X_colnames, pfts))
  
  for (item in coefs) {
    var <- item$var
    
    # var can be "tmean_CLIM" (character) or c("tmean_CLIM","precip_CLIM") (interaction tuple)
    term <- if (is.character(var) && length(var) == 1) {
      var
    } else {
      paste0(var, collapse = ":")
    }
    
    if (!term %in% X_colnames) {
      stop("Coefficient term not found in model matrix: ", term)
    }
    
    B[term, ] <- item$coef[pfts]
  }
  
  B
}
safe_softplus <- function(x) ifelse(x > 20, x, log1p(exp(x)))
# simulate biomass for each pft
#'
#' Generates simulated per-PFT and total biomass using the cwexp model
#' (cover-weighted softplus). Per-PFT noise is added on the linear predictor
#' scale; observation noise is added on the log scale of totalBio. Optionally
#' adds region-level correlated bias, where all observations in a region share
#' a common shift — this creates the spatial autocorrelation structure that
#' makes blocked CV necessary.
#'
#' @param data Data frame with cover and predictor columns.
#' @param coefs List of coefficient specifications (from `coefs_to_matrix`).
#' @param intercepts Named numeric vector of per-PFT intercepts.
#' @param pred_vars Character vector of predictor variable names.
#' @param response_var Character; name of the total biomass column.
#' @param inter List of interaction tuples (e.g., `list(c("tmean", "precip"))`).
#' @param sigma_pft Numeric >= 0; SD of noise added to the linear predictor
#'   per PFT. Makes the problem harder by perturbing per-PFT contributions.
#' @param sigma_obs Numeric >= 0; SD of log-scale observation noise on totalBio.
#' @param sigma_region Numeric >= 0; SD of region-level which is bias multiplied by linear
#'   predictor. Each region gets one random draw shared across all observations
#'   and all PFTs in that region.
#' @param region Character or numeric vector of region assignments, one per row
#'   of `data`. Required when `sigma_region > 0`.
#' @param normalize Logical; if TRUE, standardize predictors before building X.
#'
#' @return A list of class `cwexp_sim` with:
#' \describe{
#'   \item{data}{The input data with simulated columns appended: `{pft}Bio`
#'     (noiseless per-PFT mu), `totalMu` (noiseless total), and `totalBio`
#'     (observed total with all noise sources).}
#'   \item{par}{List with `alpha` (intercepts, length G), `B` (slopes only,
#'     P x G matrix without intercept row — matches cwexp fitted model
#'     structure), `sigma_pft`, `sigma_obs`, `sigma_region`.}
#'   \item{spec}{List with `formula` (response ~ predictors, no intercept
#'     term), `cover_cols`, `pred_vars`.}
#'   \item{prep}{List with `x_cols` (predictor column names) and `outcome`.}
#' }
sim_bio <- function(data,
                    coefs,
                    intercepts,
                    pred_vars,
                    response_var = 'totalBio',
                    inter = list(),
                    sigma_pft = 0,
                    sigma_obs = 0,
                    sigma_region = 0,
                    region = NULL,
                    normalize = TRUE) {
  
  pfts <- names(intercepts)
  nms <- names(data)
  stopifnot(length(pfts) == length(intercepts),
            pred_vars %in% nms)
  
  # optionally standardize predictors
  dat <- data
  if (normalize) {
    vars_to_scale <- unique(c(pred_vars, unlist(inter)))
    dat <- dplyr::mutate(dat, dplyr::across(dplyr::all_of(vars_to_scale),
                                            ~ as.numeric(scale(.x))))
  }
  
  f <- make_formula(pred_vars = pred_vars, inter = inter)
  
  X <- stats::model.matrix(f, data = dat)  # includes "(Intercept)"
  B <- coefs_to_matrix(coefs, X_colnames = colnames(X), pfts = pfts)
  
  # set intercepts
  B["(Intercept)", ] <- intercepts
  
  eta <- X %*% B  # nrow(data) x n_pfts
  
  # region-level correlated bias (same shift for all obs/PFTs in a region)
  if (sigma_region > 0) {
    if (is.null(region)) {
      stop("region must be provided when sigma_region > 0.")
    }
    
    if(is.character(region) && region %in% names(data)) {
      region <- data[[region]]
    }
    
    if (length(region) != nrow(data)) {
      stop("region must have one value per row of data.")
    }
    region_ids <- unique(region)
    region_bias <- truncnorm::rtruncnorm(length(region_ids), a = 0,
                                         mean = 1, sd = sigma_region)
    names(region_bias) <- as.character(region_ids)
    bias_vec <- region_bias[as.character(region)]
    eta <- eta * bias_vec  # broadcasts across all PFT columns
  }
  
  # per-PFT noise on the linear predictor
  # (not in the model we're trying to capture, so this 
  # 'extra' change in the data generating model, that 
  # we'll be misspecifying)
  # per-PFT noise on the linear predictor (multiplicative, proportional to signal)
  if (sigma_pft > 0) {
    pft_noise <- matrix(truncnorm::rtruncnorm(length(eta), a = 0,
                                              mean = 1, sd = sigma_pft),
                        nrow = nrow(eta), ncol = ncol(eta))
    eta_noisy <- eta * pft_noise
  } else {
    eta_noisy <- eta
  }
  # noiseless per-PFT mu: cover * softplus(eta) — the validation truth
  cols_cov <- paste0(pfts, 'Cov')
  pft_mu <- as.matrix(dat[cols_cov]) * safe_softplus(exp(eta))
  
  # noisy total mu (includes per-PFT noise, before observation noise)
  pft_noisy <- as.matrix(dat[cols_cov]) * safe_softplus(exp(eta_noisy))
  total_mu_noisy <- rowSums(pft_noisy)
  
  # observation noise on log scale of total
  if (sigma_obs > 0) {
    total_bio <- exp(log(total_mu_noisy) +
                       stats::rnorm(nrow(dat), mean = 0, sd = sigma_obs))
  } else {
    total_bio <- total_mu_noisy
  }
  
  # assemble simulated columns
  sim_cols <- tibble::as_tibble(pft_mu, .name_repair = "minimal")
  names(sim_cols) <- stringr::str_replace(names(sim_cols), 'Cov', 'Bio')
  sim_cols$totalMu <- rowSums(pft_mu)
  sim_cols[[response_var]] <- total_bio
  
  # merge simulated columns back onto input data
  # (drop original response if it exists, replace with simulated)
  out_data <- data
  if (response_var %in% names(out_data)) {
    out_data[[response_var]] <- NULL
  }
  out_data <- dplyr::bind_cols(sim_cols, out_data)
  
  # build par list matching cwexp fitted model structure:
  # alpha = intercepts, B = slopes only (no intercept row)
  B_slopes <- B[rownames(B) != "(Intercept)", , drop = FALSE]
  
  # build formula: response ~ pred1 + pred2 (+ interactions)
  f_rhs <- make_formula(pred_vars = pred_vars, inter = inter)
  f_full <- stats::as.formula(paste(response_var, "~",
                                    as.character(f_rhs)[2]))
  
  # x_cols from model matrix without intercept
  x_cols <- colnames(X)[colnames(X) != "(Intercept)"]
  
  out <- list(
    data = out_data,
    par = list(
      alpha = intercepts,
      B = B_slopes,
      sigma_pft = sigma_pft,
      sigma_obs = sigma_obs,
      sigma_region = sigma_region
    ),
    spec = list(
      formula = f_full,
      cover_cols = cols_cov,
      pred_vars = pred_vars
    ),
    prep = list(
      x_cols = x_cols,
      outcome = response_var
    )
  )
  class(out) <- "cwexp_sim"
  out
}

if(FALSE) {
  # objects for stepping through sim_bio()
  # run after sourcing Functions/data/simulate_data.R
  
  library(dplyr)
  library(stringr)
  
  pfts <- c("shrub", "tree", "grass")
  
  # minimal data with cover and predictors
  set.seed(1)
  n <- 20
  data <- tibble(
    shrubCov = runif(n, 0, 0.4),
    treeCov = runif(n, 0, 0.5),
    grassCov = runif(n, 0, 0.3),
    tmean = rnorm(n, 10, 5),
    precip = rnorm(n, 500, 200),
    region = sample(1:3, n, replace = TRUE),
    totalBio = rlnorm(n, 3, 0.5)  # placeholder, gets replaced
  )
  
  pred_vars <- c("tmean", "precip")
  inter <- list()  # no interactions for simplicity
  
  # coefficients: one entry per predictor
  coefs <- list(
    list(var = "tmean",  coef = c(shrub = 0.3, tree = -0.1, grass = 0.2)),
    list(var = "precip", coef = c(shrub = 0.1, tree = 0.4, grass = -0.2))
  )
  
  intercepts <- c(shrub = 2.0, tree = 3.0, grass = 1.5)
  
  response_var <- "totalBio"
  sigma_pft <- 0.3
  sigma_obs <- 0.2
  sigma_region <- 0.2
  region <- data$region
  normalize <- TRUE
}


#' Predict from a cwexp simulation object
#'
#' Uses the true (noiseless) parameters to compute mu on the simulation data
#' or on new data. Compatible with `compare_to_truth()` and `predict_by_group()`.
#'
#' @param object A `cwexp_sim` object from `sim_bio()`.
#' @param new_data Optional data frame. If NULL, uses `object$data`.
#' @param type `"mu"` (default) or `"log_mu"`.
#' @param ... Unused.
#'
#' @return Numeric vector of predictions.
#' @export
predict.cwexp_sim <- function(object, newdata = NULL,
                              type = c("mu", "log_mu"), ...) {
  type <- match.arg(type)
  data <- if (is.null(newdata)) object$data else newdata
  
  prep <- cwexp_prepare(
    data = data,
    formula = object$spec$formula,
    cover_cols = object$spec$cover_cols
  )
  
  mu <- cwexp_mu(
    alpha = object$par$alpha,
    B = object$par$B,
    X = prep$X,
    C = prep$C
  )
  
  if (type == "mu") return(mu)
  log(mu)
}


# make a synthetic dataset consistent with the model ----
# (quick fully artificial dataset for simple testing)

#' Generate synthetic cwexp data for testing
#'
#' Creates a dataset with known cwexp parameters for model fitting verification.
#' Optionally adds collinear junk predictors (true beta = 0) to test
#' regularization.
#'
#' @param n Integer; number of observations.
#' @param cover_cols Character vector of cover column names.
#' @param x_cols Character vector of real predictor names.
#' @param sigma Numeric; log-scale SD for the lognormal observation model.
#' @param n_junk_pred Integer; number of junk (zero-effect) predictors to add.
#'   Each junk predictor is generated as a noisy copy of one of the real
#'   predictors (cycling through them), so they are collinear with the true
#'   signal but have no real effect. Default 0 (no junk predictors).
#'
#' @return A list of class `cwexp_dummy` with elements `data`, `par`, `spec`,
#'   and `prep`. When `n_junk_pred > 0`, the B matrix has additional rows of
#'   zeros, and the formula/x_cols include the junk predictor names.
cwexp_make_dummy_data <- function(n = 500,
                                  cover_cols = paste0("cov_g", 1:6),
                                  x_cols = c("tmean", "ppt", "vpd", "sand"),
                                  sigma = 0.6,
                                  n_junk_pred = 0) {
  G <- length(cover_cols)
  P <- length(x_cols)
  
  stopifnot(n_junk_pred >= 0)
  
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
  
  # --- junk predictors (correlated with real, true effect = 0) -------------
  all_x_cols <- x_cols
  
  if (n_junk_pred > 0) {
    junk_names <- paste0("junk", seq_len(n_junk_pred))
    X_real <- as.matrix(Xdat)
    
    Xjunk <- matrix(NA_real_, nrow = n, ncol = n_junk_pred)
    for (j in seq_len(n_junk_pred)) {
      # cycle through real predictors as "parents"
      parent_col <- ((j - 1) %% P) + 1
      noise_sd <- sd(X_real[, parent_col]) * 0.5
      Xjunk[, j] <- X_real[, parent_col] + rnorm(n, sd = noise_sd)
    }
    colnames(Xjunk) <- junk_names
    Xdat <- cbind(Xdat, as.data.frame(Xjunk))
    
    # extend B with zero rows for junk predictors
    B_junk <- matrix(0, nrow = n_junk_pred, ncol = G)
    rownames(B_junk) <- junk_names
    colnames(B_junk) <- cover_cols
    B_true <- rbind(B_true, B_junk)
    
    all_x_cols <- c(x_cols, junk_names)
  }
  
  # model mean
  X <- as.matrix(Xdat)
  eta <- X %*% B_true
  eta <- sweep(eta, 2, alpha_true, "+")
  mu <- rowSums(C * log1p (exp(eta)))
  
  # lognormal response
  y <- exp(rnorm(n, mean = log(mu), sd = sigma))
  
  dat <- cbind.data.frame(totalBio = y, Xdat, as.data.frame(C))
  formula = as.formula(paste0('totalBio ~ ', paste0(all_x_cols, collapse = ' + ')))
  
  out <- list(
    data = dat,
    par = list(alpha = alpha_true, B = B_true, sigma = sigma),
    spec = list(cover_cols = cover_cols,
                formula = formula),
    prep = list(x_cols = all_x_cols,
                outcome = 'totalBio')
    
  )
  class(out) <- 'cwexp_dummy'
  out
}

# 
#

