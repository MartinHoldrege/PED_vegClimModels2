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
