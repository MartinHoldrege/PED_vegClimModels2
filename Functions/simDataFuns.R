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
                    normalize = FALSE) {
  
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



