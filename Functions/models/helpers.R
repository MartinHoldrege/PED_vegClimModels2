#' Check for high collinearity among predictors
#'
#' Computes pairwise correlations among predictor columns and throws an
#' error if any pair exceeds the threshold. Accepts either a formula
#' (which gets expanded via `model.matrix`, capturing interactions) or a
#' character vector of column names.
#'
#' @param data Data frame containing predictor columns.
#' @param formula Optional model formula (e.g., `~ x1 + x2 + x1:x2`). If
#'   provided, `model.matrix` is used to expand interactions. The intercept
#'   is automatically removed.
#' @param pred_vars Optional character vector of predictor column names. Used
#'   when `formula` is NULL. Ignored if `formula` is provided.
#' @param threshold Numeric; absolute correlation above which to flag.
#'   Default 0.9.
#' @param action Character; `"error"` (default) stops execution, `"warn"`
#'   issues a warning instead.
#'
#' @return Invisibly returns the correlation matrix. Throws an error or
#'   warning if any pair exceeds the threshold.
#' @examples
#' df <- data.frame(a = rnorm(100), b = rnorm(100))
#' df$c <- df$a + rnorm(100, sd = 0.1)
#'
#' # with pred_vars
#' check_collinearity(df, pred_vars = c("a", "b", "c"), threshold = 0.9)
#'
#' # with formula (catches interaction collinearity)
#' check_collinearity(df, formula = ~ a + b + a:b, threshold = 0.9)
#' @export
check_collinearity <- function(data, formula = NULL, pred_vars = NULL,
                               threshold = 0.9,
                               action = c("error", "warn")) {
  action <- match.arg(action)
  
  if (!is.null(formula)) {
    # strip response if present (e.g., y ~ x1 + x2 becomes ~ x1 + x2)
    if (length(formula) == 3) {
      formula <- formula[-2]
    }
    X <- stats::model.matrix(formula, data = data)
    # drop intercept
    if (ncol(X) >= 1 && colnames(X)[1] == "(Intercept)") {
      X <- X[, -1, drop = FALSE]
    }
    col_names <- colnames(X)
  } else if (!is.null(pred_vars)) {
    stopifnot(all(pred_vars %in% names(data)))
    X <- as.matrix(data[, pred_vars])
    col_names <- pred_vars
  } else {
    stop("Either formula or pred_vars must be provided.")
  }
  
  cor_mat <- cor(X, use = "pairwise.complete.obs")
  
  # upper triangle only (avoid duplicates and diagonal)
  cor_mat[!upper.tri(cor_mat)] <- 0
  high <- which(abs(cor_mat) > threshold, arr.ind = TRUE)

  if (nrow(high) > 0) {
    pairs <- apply(high, 1, function(idx) {
      paste0(col_names[idx[1]], " & ", col_names[idx[2]],
             " (r = ", round(cor_mat[idx[1], idx[2]], 3), ")")
    })
    msg <- paste0("High collinearity detected (|r| > ", threshold, "):\n  ",
                  paste(pairs, collapse = "\n  "))
    if (action == "error") stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }
  
  invisible(cor_mat)
}

safe_softplus <- function(x) ifelse(x > 20, x, log1p(exp(x)))

# where sub_spec is a list of model specifications
make_model_formula <- function(sub_spec) {
  pred_vars <- sub_spec$pred_vars
  inter <- sub_spec$inter
  
  inter_terms <- if (!is.null(inter)) {
    purrr::map_chr(inter, ~ paste(.x, collapse = ":"))
  } else {
    character(0)
  }
  all_terms <- c(pred_vars, inter_terms)
  as.formula(paste("totalBio ~", paste(all_terms, collapse = " + ")))
}


# helper: build config from sub-spec and settings
# for 01_fit_model_hw.R
build_config <- function(vm, sub_spec, formula,
                         cv_settings, clustering_settings, vd, vp,
                         purer_spec,
                         model_type) {
  list(
    vm = vm,
    model = list(
      formula = formula,
      cover_cols = sub_spec$cover_cols,
      pred_vars = sub_spec$pred_vars,
      dll_path = sub_spec$dll_path,
      fix_alpha_pfts = sub_spec$fix_alpha_pfts,
      fix_alpha_filter = sub_spec$fix_alpha_filter,
      type = model_type
    ),
    # purer pixel selection
    purer = list(
      vp = vp, # version of purer selection
      region_col = "region",
      q = purer_spec$q,
      min_raw_cover = purer_spec$min_raw_cover,
      n_sample = purer_spec$n_sample,
      max_n_per_region = NULL,  # NULL = no cap
      seed = 42
    ),
    cv = cv_settings,
    clustering = clustering_settings,
    data = list(version = vd, type = model_type)
  )
}
