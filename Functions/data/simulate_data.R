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

# simulate biomass for each pft 
sim_bio <- function(data,
                    coefs,
                    intercepts,
                    pred_vars,
                    response_var = 'totalBio', 
                    inter = list(),
                    sigma = 0,
                    normalize = TRUE) {
  
  pfts <- names(intercepts)
  nms <- names(data)
  stopifnot(length(pfts) == length(intercepts),
            pred_vars %in% nms)
  # optionally standardize predictors (and any interaction components will reflect standardized inputs)
  dat <- data
  if (normalize) {
    vars_to_scale <- unique(c(pred_vars, unlist(inter)))
    dat <- dplyr::mutate(dat, dplyr::across(dplyr::all_of(vars_to_scale), ~ as.numeric(scale(.x))))
  }
  
  f <- make_formula(pred_vars = pred_vars, inter = inter)
  
  X <- stats::model.matrix(f, data = dat)  # includes "(Intercept)"
  B <- coefs_to_matrix(coefs, X_colnames = colnames(X), pfts = pfts)
  
  # set intercepts
  B["(Intercept)", ] <- intercepts
  
  eta <- X %*% B  # nrow(data) x n_pfts # linear predictions
  
  if (sigma > 0) {
    eta <- eta + matrix(stats::rnorm(length(eta), mean = 0, sd = sigma),
                        nrow = nrow(eta), ncol = ncol(eta))
  }
  
  eta2 <- exp(eta)
  cols_cov <- paste0(pfts, 'Cov')
  eta3 <- dat[cols_cov]*eta2
 
  
  # return as tibble, with one column per PFT
  out <- tibble::as_tibble(eta3, .name_repair = "minimal")
  names(out) <- str_replace(names(out), 'Cov', 'Bio')
  out[[response_var]] = rowSums(eta3)
  
  out
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

