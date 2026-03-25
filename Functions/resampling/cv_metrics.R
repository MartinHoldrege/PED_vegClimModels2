# metrics for calculating cross validation performance


#  functions of metrics ---------------------------------------------------


#' Mean absolute error
#'
#' Computes mean absolute error between truth and estimate.
#'
#' @param truth Numeric vector of observed values.
#' @param estimate Numeric vector of predicted values.
#'
#' @return Numeric scalar.
metric_mae <- function(truth, estimate) {
  mean(abs(truth - estimate))
}

#' Root mean squared error
#'
#' Computes root mean squared error between truth and estimate.
#'
#' @param truth Numeric vector of observed values.
#' @param estimate Numeric vector of predicted values.
#'
#' @return Numeric scalar.
metric_rmse <- function(truth, estimate) {
  sqrt(mean((truth - estimate)^2))
}

#' Mean absolute error on log scale
#'
#' Computes MAE between log-transformed truth and estimate.
#'
#' @param truth Numeric vector of observed values. Must be > 0.
#' @param estimate Numeric vector of predicted values. Must be > 0.
#'
#' @return Numeric scalar.
metric_mae_log <- function(truth, estimate) {
  mean(abs(log(truth) - log(estimate)))
}

# log(x + 1) scale
metric_mae_log1p <- function(truth, estimate) {
  mean(abs(log1p(truth) - log1p(estimate)))
}

#' RMSE on log scale
#'
#' Computes RMSE between log-transformed truth and estimate.
#'
#' @param truth Numeric vector of observed values. Must be > 0.
#' @param estimate Numeric vector of predicted values. Must be > 0.
#'
#' @return Numeric scalar.
metric_rmse_log <- function(truth, estimate) {
  sqrt(mean((log(truth) - log(estimate))^2))
}

# log(x + 1) scale
metric_rmse_log1p <- function(truth, estimate) {
  sqrt(mean((log1p(truth) - log1p(estimate))^2))
}

#' Correlation between truth and estimate
#'
#' Computes Pearson correlation.
#'
#' @param truth Numeric vector of observed values.
#' @param estimate Numeric vector of predicted values.
#'
#' @return Numeric scalar.
metric_cor <- function(truth, estimate) {
  stats::cor(truth, estimate)
}

#' Get a metric function by name
#'
#' Returns a metric function for use in tuning and scoring.
#'
#' @param metric Character scalar naming the metric.
#'
#' @return A function with arguments `truth` and `estimate`.
get_metric_fun <- function(metric = c("mae_log1p", "rmse_log1p",
                                      "mae_log", "rmse_log", "mae", 
                                      "rmse", "cor")) {
  metric <- match.arg(metric)
  
  switch(
    metric,
    mae_log1p = metric_mae_log1p,
    rmse_log1p = metric_rmse_log1p,
    mae_log = metric_mae_log,
    rmse_log = metric_rmse_log,
    mae = metric_mae,
    rmse = metric_rmse,
    cor = metric_cor
  )
}


# calculating metrics on holdout data -------------------------------------

#' Score a fitted model on a holdout dataset
#'
#' Computes one or more prediction metrics for a fitted model on a holdout
#' dataset. By default, the outcome name is taken from `fit$prep$outcome`
#' if available.
#'
#' @param fit A fitted model object.
#' @param newdata Data frame used for holdout scoring.
#' @param outcome Optional character scalar naming the observed response column.
#'   If `NULL`, the function tries to use `fit$prep$outcome`.
#' @param metrics Character vector of metric names.
#' @param predict_fun Function used to generate predictions.
#' @param predict_type Character scalar passed to `predict_fun`.
#'
#' @return Named list of metric values.
score_fit <- function(fit,
                      newdata,
                      outcome = NULL,
                      metrics = c("mae_log1p", "rmse_log1p",
                                  "mae_log", "rmse_log"),
                      predict_fun = predict,
                      predict_type = "mu") {
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame.")
  }
  
  if (is.null(outcome)) {
    outcome <- fit$prep$outcome
  }
  
  if (is.null(outcome) || length(outcome) != 1 || !is.character(outcome)) {
    stop("outcome must be provided or available as fit$prep$outcome.")
  }
  
  if (!outcome %in% names(newdata)) {
    stop("outcome column not found in newdata.")
  }
  
  truth <- newdata[[outcome]]
  estimate <- predict_fun(fit, newdata, type = predict_type)
  
  out <- vector("list", length(metrics))
  names(out) <- metrics
  
  for (i in seq_along(metrics)) {
    metric_fun <- get_metric_fun(metrics[i])
    out[[i]] <- metric_fun(truth = truth, estimate = estimate)
  }
  
  out
}

#' Score a fitted model on new data and return a one-row data frame
#'
#' Thin wrapper around `score_fit()`.
#'
#' @param fit A fitted model object.
#' @param newdata Data frame used for holdout scoring.
#' @param ... Additional arguments passed to `score_fit()`.
#'
#' @return One-row data frame of score metrics.
score_fit_df <- function(fit, newdata, ...) {
  scores <- score_fit(fit = fit, newdata = newdata, ...)
  
  if (!is.list(scores)) {
    stop("score_fit() must return a named list.")
  }
  
  out <- as.data.frame(scores, stringsAsFactors = FALSE)
  
  if (nrow(out) != 1L) {
    stop("score_fit_df() expected score_fit() to produce one row.")
  }
  
  out
}



#' Score a fitted model path on a holdout dataset
#'
#' Computes holdout metrics for each fitted model in a path object.
#'
#' @param path_fit A path-fit object containing elements `fits` and `summary`.
#' such as produced by the cwexp_fit_lambda_path_tmb function
#' @param newdata Data frame used for holdout scoring.
#' @param outcome Optional character scalar naming the observed response column.
#'   If `NULL`, the function tries to use `fit$prep$outcome` for each fit.
#' @param metrics Character vector of metric names.
#' @param predict_fun Function used to generate predictions.
#' @param predict_type Character scalar passed to `predict_fun`.
#'
#' @return Data frame with one row per fitted model and one column per metric.
score_path <- function(path_fit,
                       newdata,
                       outcome = NULL,
                       metrics = c("mae_log1p", "rmse_log1p", "mae_log", 
                                   "rmse_log", "mae", "rmse", "cor"),
                       predict_fun = predict,
                       predict_type = "mu") {
  if (!is.list(path_fit) || is.null(path_fit$fits) || is.null(path_fit$summary)) {
    stop("path_fit must be a list with elements 'fits' and 'summary'.")
  }
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame.")
  }
  
  n_fit <- length(path_fit$fits)
  out <- vector("list", n_fit)
  
  for (i in seq_along(path_fit$fits)) {
    fit_i <- path_fit$fits[[i]]
    
    score_i <- score_fit(
      fit = fit_i,
      newdata = newdata,
      outcome = outcome,
      metrics = metrics,
      predict_fun = predict_fun,
      predict_type = predict_type
    )
    
    row_i <- as.list(path_fit$summary[i, , drop = FALSE])
    out[[i]] <- c(row_i, score_i)
  }
  
  out_df <- dplyr::bind_rows(out)
  rownames(out_df) <- NULL
  out_df
}


# selecting the best scoring ----------------------------------------------


#' Select lambda by minimum metric
#'
#' @param score_df Data frame containing lambda values and metric columns
#'  from score_path()
#' @param metric Character scalar naming the metric column to optimize.
#'
#' @return One-row data frame corresponding to the selected lambda.
select_lambda_min_metric <- function(score_df, metric = "mae_log1p") {
  if (!"lambda" %in% names(score_df)) {
    stop("score_df must contain a 'lambda' column.")
  }
  if (!metric %in% names(score_df)) {
    stop("metric column not found in score_df.")
  }
  
  vals <- score_df[[metric]]
  best <- min(vals)
  
  candidates <- score_df[vals == best, , drop = FALSE]
  
  # if ties choose the one with stronger penalization
  candidates[which.max(candidates$lambda), , drop = FALSE]
}

#' Select largest lambda within tolerance of the best metric
#'
#' Selects the largest lambda whose metric is within a specified tolerance of
#' the minimum metric value.
#'
#' @param score_df Data frame containing lambda values and metric columns.
#' @param metric Character scalar naming the metric column to optimize.
#' @param tol Numeric scalar >= 0. Allowed deviation from the best metric.
#'
#' @return One-row data frame corresponding to the selected lambda.
select_lambda_within_tol <- function(score_df,
                                     metric = "mae_log1p",
                                     tol = 0) {
  if (!"lambda" %in% names(score_df)) {
    stop("score_df must contain a 'lambda' column.")
  }
  if (!metric %in% names(score_df)) {
    stop("metric column not found in score_df.")
  }
  if (!is.finite(tol) || tol < 0) {
    stop("tol must be a finite number >= 0.")
  }
  
  vals <- score_df[[metric]]
  best <- min(vals)
  
  keep <- vals <= best + tol
  score_keep <- score_df[keep, ]
  
  # among acceptable fits, prefer the one with stronger penalization
  score_keep[which.max(score_keep$lambda), ]
}

#' Select lambda using the 1-SE rule
#'
#' Selects the largest (most regularized) lambda whose mean metric is within
#' one standard error of the best mean metric. This is the standard "1-SE rule"
#' which favors parsimony when models are statistically indistinguishable.
#'
#' @param score_df Data frame from `summarize_scores()`, containing columns
#'   for `lambda`, the metric, and `{metric}_se`.
#' @param metric Character scalar naming the metric column to optimize.
#'
#' @return One-row data frame corresponding to the selected lambda.
#' @export
select_lambda_1se <- function(score_df, metric = "mae_log1p") {
  if (!"lambda" %in% names(score_df)) {
    stop("score_df must contain a 'lambda' column.")
  }
  if (!metric %in% names(score_df)) {
    stop("metric column not found in score_df.")
  }
  
  se_col <- paste0(metric, "_se")
  if (!se_col %in% names(score_df)) {
    stop("SE column '", se_col, "' not found in score_df. ",
         "Make sure score_df comes from summarize_scores().")
  } 
  
  ses  <- score_df[[se_col]]
  vals <- score_df[[metric]]
  
  best_idx <- which.min(vals)
  threshold <- vals[best_idx] + ses[best_idx]
  
  # candidates: all lambdas whose mean metric is within 1 SE of best
  keep <- vals <= threshold
  score_keep <- score_df[keep, , drop = FALSE]
  
  # among those, pick the largest lambda (most regularized)
  score_keep[which.max(score_keep$lambda), , drop = FALSE]
}


#' Select lambda from a scored path
#'
#' Applies a named selection rule to a scored lambda path.
#'
#' @param score_df Data frame containing lambda values and metric columns.
#' @param metric Character scalar naming the metric column to optimize.
#' @param rule Character scalar. One of `"min"`, `"largest_within_tol"`,
#'   or `"1se"`.
#' @param tol Numeric scalar >= 0. Used only for `"largest_within_tol"`.
#'
#' @return One-row data frame corresponding to the selected lambda.
#' @examples
#' # simulate a score_df like summarize_scores() would produce
#' score_df <- data.frame(
#'   lambda     = c(1.0, 0.5, 0.1, 0.01, 0),
#'   mae_log    = c(0.60, 0.52, 0.48, 0.47, 0.50),
#'   mae_log_se = c(0.05, 0.04, 0.03, 0.03, 0.04)
#' )
#'
#' # minimum: picks lambda with lowest mean metric
#' select_lambda(score_df, metric = "mae_log", rule = "min")
#'
#' # 1se: largest lambda within 1 SE of best
#' select_lambda(score_df, metric = "mae_log", rule = "1se")
#'
#' # within tolerance: largest lambda within fixed tol of best
#' select_lambda(score_df, metric = "mae_log", rule = "largest_within_tol",
#'               tol = 0.02)
#' @export
select_lambda <- function(score_df,
                          metric = "mae_log1p",
                          rule = c("1se", "min", "largest_within_tol"),
                          tol = 0) {
  rule <- match.arg(rule)
  
  if (rule == "min") {
    return(select_lambda_min_metric(score_df = score_df, metric = metric))
  }
  
  if (rule == "1se") {
    return(select_lambda_1se(score_df = score_df, metric = metric))
  }
  
  select_lambda_within_tol(
    score_df = score_df,
    metric = metric,
    tol = tol
  )
}

# misc --------------------------------------------------------------------


summarize_scores <- function(scores, metric_cols) {
  se <- function(x) sd(x)/sqrt(length(x))
  score_summary <- scores |>
    dplyr::group_by(.data$lambda) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(metric_cols),
        se,
        .names = "{.col}_se"
      ),
      dplyr::across(
        dplyr::all_of(metric_cols),
        mean,
        .names = "{.col}"
      ),
      n = dplyr::n(),
      .groups = "drop"
    )
  score_summary$summary_stat <- "mean"
  score_summary
}

cv_metrics_df <- function(data, observed = 'observed',
                             predicted = 'predicted') {
  
  stopifnot(
    c(observed, predicted) %in% names(data)
  )
  obs <- data[[observed]]
  pred <- data[[predicted]]
  
  metrics_df <- tibble(
    metric = c("MAE (log1p)", "RMSE (log1p)", "MAE (log)", "RMSE (log)", 
               "MAE", "RMSE", "Correlation", "N"),
    value = c(
      metric_mae_log1p(obs, pred),
      metric_rmse_log1p(obs, pred),
      metric_mae_log(obs, pred),
      metric_rmse_log(obs, pred),
      metric_mae(obs, pred),
      metric_rmse(obs, pred),
      metric_cor(obs, pred),
      length(obs)
    )
  ) 
  metrics_df
}


# compare to truth --------------------------------------------------------
# functions relevant for determining recovery of simulated data

#' Align two B matrices by row and column names
#'
#' Aligns a fitted B matrix to a truth B matrix by matching row names
#' (predictors) and column names (groups). Predictors in the fitted model
#' that don't exist in truth are compared against 0. Predictors in truth
#' that aren't in the fitted model are excluded from the comparison.
#'
#' @param B_true Matrix with named rows and columns (from simulation).
#' @param B_fit Matrix with named rows and columns (from fitted model).
#'
#' @return A list with `true` and `fit` matrices of the same dimensions,
#'   aligned by name, plus `extra_in_fit` (predictor names in fit but not
#'   truth, compared against 0) and `missing_in_fit` (predictor names in
#'   truth but not fit, excluded).
#' @examples A tibble with columns: model, group, rmse, cor.
#' # truth has 2 predictors x 3 groups
#' B_true <- matrix(c(0.3, 0.1, -0.2, 0.4, 0.1, -0.1, 0.1, 0.1, 0.1),
#'                  nrow = 3, ncol = 3,
#'                  dimnames = list(c("tmean", "precip", "other_var"),
#'                                 c("shrubCov", "treeCov", "grassCov")))
#'
#' # fitted model has 3 predictors (extra "junk1") x same 3 groups
#' B_fit <- matrix(c(0.28, 0.09, 0.02, -0.18, 0.38, 0.01,
#'                   0.12, -0.09, -0.01),
#'                 nrow = 3, ncol = 3,
#'                 dimnames = list(c("tmean", "precip", "junk1"),
#'                                c("shrubCov", "treeCov", "grassCov")))
#'
#' aligned <- align_B_matrices(B_true, B_fit)
#' aligned$true   # 3 x 3: tmean, precip, junk1 — junk1 row is 0
#' aligned$fit    # 3 x 3: all fitted values
#' aligned$extra_in_fit   # "junk1"
#' aligned$missing_in_fit # character(0)
align_B_matrices <- function(B_true, B_fit) {
  stopifnot(is.matrix(B_true), is.matrix(B_fit))
  
  if (is.null(rownames(B_true)) || is.null(rownames(B_fit))) {
    stop("Both B matrices must have row names (predictor names).")
  }
  if (is.null(colnames(B_true)) || is.null(colnames(B_fit))) {
    stop("Both B matrices must have column names (group/cover names).")
  }
  
  # align columns (groups) — use fitted model's column order
  shared_cols <- intersect(colnames(B_fit), colnames(B_true))
  if (length(shared_cols) == 0) {
    stop("No shared column names between B_true and B_fit.")
  }
  
  # align rows (predictors) — union of both, filling missing with 0
  all_rows <- union(rownames(B_fit), rownames(B_true))
  extra_in_fit <- setdiff(rownames(B_fit), rownames(B_true))
  missing_in_fit <- setdiff(rownames(B_true), rownames(B_fit))
  
  # build aligned matrices
  B_true_aligned <- matrix(0, nrow = length(all_rows), ncol = length(shared_cols),
                           dimnames = list(all_rows, shared_cols))
  B_fit_aligned <- B_true_aligned
  
  # fill in values that exist
  B_true_aligned[rownames(B_true), shared_cols] <- B_true[rownames(B_true), shared_cols]
  B_fit_aligned[rownames(B_fit), shared_cols] <- B_fit[rownames(B_fit), shared_cols]
  
  list(
    true = B_true_aligned,
    fit = B_fit_aligned,
    extra_in_fit = extra_in_fit,
    missing_in_fit = missing_in_fit
  )
}


#' Compare fitted models to truth (simulated or observed) for total mu
#'
#' Computes parameter- and prediction-level metrics comparing fitted models to
#' a known truth object or observed data columns. When observed `totalMu` is
#' present in the data, it can be used in place of predicted truth mu.
#'
#' @param truth Object with true parameters (`cwexp_sim`, `cwexp_dummy`, or
#'   similar list with `par` and `data`).
#' @param fits Named list of fitted model objects.
#' @param newdata Optional data frame. If NULL, uses `truth$data`. Must contain
#'   `totalMu` column if `use_observed = "force"`.
#' @param use_observed One of `"auto"` (use observed `totalMu` if present,
#'   otherwise predict from truth), `"force"` (require observed column, error
#'   if absent), or `"never"` (always predict from truth parameters).
#'
#' @return A tibble with one row per model and columns for each metric.
#' @export
compare_to_truth <- function(truth, fits, newdata = NULL,
                             use_observed = "auto") {
  stopifnot(is.list(truth), !is.null(truth$par),
            is.list(fits), length(fits) >= 1,
            use_observed %in% c("auto", "force", "never")
            )
  
  if (is.null(names(fits))) {
    names(fits) <- paste0("model_", seq_along(fits))
  }
  
  data <- if (!is.null(newdata)) newdata else truth$data
  alpha_true <- truth$par$alpha
  B_true <- truth$par$B
  
  has_obs <- "totalMu" %in% names(data)                      
  if (use_observed == "force" && !has_obs) {                   
    stop("use_observed = 'force' but 'totalMu' not found in data")  
  }                                                             
  use_obs_mu <- has_obs && use_observed != "never"             
  mu_true <- if (use_obs_mu) {                                 
    data$totalMu                                               
  } else {                                                     
    predict(truth, newdata = data, type = "mu")
  }   
  
  rows <- lapply(names(fits), function(nm) {
    fit <- fits[[nm]]
    mu_hat <- predict(fit, data, type = "mu")
    
    # align alpha by name (use cover_cols order from fit)
    fit_alpha <- fit$par$alpha
    if (!is.null(names(alpha_true)) && !is.null(names(fit_alpha))) {
      shared_groups <- intersect(names(fit_alpha), names(alpha_true))
      alpha_t <- alpha_true[shared_groups]
      alpha_f <- fit_alpha[shared_groups]
    } else {
      # fall back to positional if no names
      alpha_t <- alpha_true
      alpha_f <- fit_alpha
    }
    
    # align B matrices by name
    aligned <- align_B_matrices(B_true, fit$par$B)
    
    dplyr::tibble(
      model = nm,
      mu_source    = if (use_obs_mu) "observed" else "predicted",  
      alpha_rmse = metric_rmse(alpha_t, alpha_f),
      alpha_cor  = cor(alpha_t, alpha_f),
      B_rmse     = metric_rmse(as.numeric(aligned$true), as.numeric(aligned$fit)),
      B_cor      = cor(as.numeric(aligned$true), as.numeric(aligned$fit)),
      mu_rmse    = metric_rmse(mu_true, mu_hat),
      mu_cor     = cor(mu_true, mu_hat),
      mu_mae_log = metric_mae_log(mu_true, mu_hat),
      B_extra    = paste(aligned$extra_in_fit, collapse = ", "),
      B_missing  = paste(aligned$missing_in_fit, collapse = ", ")
    )
  })
  
  dplyr::bind_rows(rows)
}



#' Compare fitted models to truth (simulated or observed) for per-group mu
#'
#' Computes per-group prediction metrics comparing fitted models to a truth
#' object or observed per-group biomass columns. When `{group}Bio` columns are
#' present in the data they can be used in place of predicted group-level truth.
#' predicted group-level truth can be different from 'observed' if there
#' are some structural complications added in the simulation of truth (e.g. regional
#' intercepts differ around the 'true' or mean intercept)
#'
#' @param truth Object with true parameters (`cwexp_sim`, `cwexp_dummy`, or
#'   similar).
#' @param fits Named list of fitted model objects.
#' @param newdata Optional data frame. If NULL, uses `truth$data`. Must contain
#'   `{group}Bio` columns if `use_observed != "never"`.
#' @param weighted Logical. If TRUE, uses cover-weighted group mu.
#' @param use_observed One of `"auto"` (use `{group}Bio` columns if present,
#'   otherwise predict from truth), `"force"` (require observed columns, error
#'   if any are absent), or `"never"` (always predict from truth parameters).
#'
#' @return A tibble with one row per model-group combination and metric columns.
#' @export
compare_group_to_truth <- function(truth, fits, newdata = NULL,
                                   weighted = TRUE,
                                   use_observed = "auto") {    
  stopifnot(is.list(truth), !is.null(truth$par))
  stopifnot(is.list(fits), length(fits) >= 1)
  stopifnot(use_observed %in% c("auto", "force", "never"))  
  
  if (is.null(names(fits))) {
    names(fits) <- paste0("model_", seq_along(fits))
  }
  
  data        <- if (!is.null(newdata)) newdata else truth$data
  group_names <- names(truth$par$alpha)
  
  # --- resolve per-group truth ---                             
  if (use_observed != "never") {                              
    bio_cols    <- paste0(stringr::str_replace(group_names, 'Cov$', ''), 
                          "Bio")                 
    present     <- bio_cols %in% names(data)                   
    if (use_observed == "force" && !all(present)) {            
      missing <- bio_cols[!present]                            
      stop("use_observed = 'force' but columns not found: ",    
           paste(missing, collapse = ", "))                     
    }                                                           
    use_obs_grp <- any(present)                                 
  } else {                                                      
    use_obs_grp <- FALSE                                        
  }                                                             
  
  mu_true_grp <- if (use_obs_grp) {                            
    # build matrix from observed columns; NA for absent groups  
    bio_cols <- paste0(stringr::str_replace(group_names, 'Cov$', ''), "Bio")                      
    mat <- sapply(bio_cols, function(col) {                    
      if (col %in% names(data)) data[[col]] else rep(NA_real_, nrow(data)) 
    })                                                        
    colnames(mat) <- group_names                              
    mat                                                      
  } else {                                                     
    predict_by_group(truth, newdata = data, weighted = weighted)
  }                                                           
  
  rows <- lapply(names(fits), function(nm) {
    fit        <- fits[[nm]]
    mu_hat_grp <- predict_by_group(fit, newdata = data, weighted = weighted)
    
    dplyr::bind_rows(lapply(group_names, function(g) {
      true_vec <- mu_true_grp[, g]
      pred_vec <- mu_hat_grp[, g]
      ok       <- is.finite(true_vec) & is.finite(pred_vec)
      
      dplyr::tibble(
        model      = nm,
        group      = stringr::str_replace(g, "Cov", ''),
        mu_source  = if (use_obs_grp && !is.na(mu_true_grp[1, g])) 
          "observed" else "predicted",                 
        n_obs      = sum(ok),
        rmse       = metric_rmse(true_vec[ok], pred_vec[ok]),
        cor        = cor(true_vec[ok], pred_vec[ok]),
        mae_log    = metric_mae_log(true_vec[ok], pred_vec[ok])
      )
    }))
  })
  
  dplyr::bind_rows(rows)
}



