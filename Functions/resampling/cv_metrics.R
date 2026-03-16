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
get_metric_fun <- function(metric = c("mae_log", "rmse_log", "mae", "rmse", "cor")) {
  metric <- match.arg(metric)
  
  switch(
    metric,
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
                      metrics = c("mae_log", "rmse_log"),
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
                       metrics = c("mae_log", "rmse_log"),
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
select_lambda_min_metric <- function(score_df, metric = "mae_log") {
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
                                     metric = "mae_log",
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
select_lambda_1se <- function(score_df, metric = "mae_log") {
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
  
  vals <- score_df[[metric]]
  ses  <- score_df[[se_col]]
  
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
                          metric = "mae_log",
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
  se <- function(x) sd(x)/length(x)
  score_summary <- scores |>
    dplyr::group_by(.data$lambda) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(metric_cols),
        mean,
        .names = "{.col}"
      ),
      dplyr::across(
        dplyr::all_of(metric_cols),
        se,
        .names = "{.col}_se"
      ),
      n = dplyr::n(),
      .groups = "drop"
    )
  score_summary$summary_stat <- "mean"
  score_summary
}



#' Compare fitted models to known truth
#'
#' Given the true parameters (from simulation) and one or more fitted model
#' objects, computes recovery metrics for alpha, B, and mu. Intended for
#' evaluating model performance on simulated data where truth is known.
#'
#' @param truth Object of class `cwexp_dummy` (from `cwexp_make_dummy_data()`),
#'   or a list with elements `par` (containing `alpha`, `B`) and `data`.
#' @param fits Named list of fitted model objects (e.g., `cwexp_tmb_fit`), each
#'   with a `par` sub-list containing `alpha` and `B`.
#' @param newdata Optional data frame for computing mu predictions. If NULL,
#'   uses `truth$data`.
#'
#' @return A tibble with one row per model and columns for each metric.
#' @examples
#' \dontrun{
#' dummy <- cwexp_make_dummy_data(n = 500)
#' # ... fit models ...
#' compare_to_truth(
#'   truth = dummy,
#'   fits = list(unpenalized = fit0, regularized = final_fit)
#' )
#' }
compare_to_truth <- function(truth, fits, newdata = NULL) {
  stopifnot(is.list(truth), !is.null(truth$par))
  stopifnot(is.list(fits), length(fits) >= 1)
  
  if (is.null(names(fits))) {
    names(fits) <- paste0("model_", seq_along(fits))
  }
  
  data <- if (!is.null(newdata)) newdata else truth$data
  mu_true <- predict(truth, newdata = data, type = "mu")
  alpha_true <- truth$par$alpha
  B_true <- as.numeric(truth$par$B)
  
  rows <- lapply(names(fits), function(nm) {
    fit <- fits[[nm]]
    mu_hat <- predict(fit, data, type = "mu")
    
    dplyr::tibble(
      model = nm,
      alpha_rmse = metric_rmse(alpha_true, fit$par$alpha),
      alpha_cor  = cor(alpha_true, fit$par$alpha),
      B_rmse     = metric_rmse(B_true, as.numeric(fit$par$B)),
      B_cor      = cor(B_true, as.numeric(fit$par$B)),
      mu_rmse    = metric_rmse(mu_true, mu_hat),
      mu_cor     = cor(mu_true, mu_hat),
      mu_mae_log = metric_mae_log(mu_true, mu_hat)
    )
  })
  
  dplyr::bind_rows(rows)
}



#' Compare per-group mu recovery across models
#'
#' Computes RMSE and correlation between true and predicted per-group mu
#' for one or more fitted models. Uses cover-weighted group contributions
#' by default.
#'
#' @param truth Object with true parameters (`cwexp_dummy` or similar).
#' @param fits Named list of fitted model objects.
#' @param newdata Optional data frame. If NULL, uses `truth$data`.
#' @param weighted Logical. Passed to `predict_by_group`.
#'
#' @return A tibble with columns: model, group, rmse, cor.
#' @export
compare_group_to_truth <- function(truth, fits, newdata = NULL,
                                   weighted = TRUE) {
  stopifnot(is.list(fits), length(fits) >= 1)
  if (is.null(names(fits))) names(fits) <- paste0("model_", seq_along(fits))
  
  data <- if (!is.null(newdata)) newdata else truth$data
  
  mu_true_grp <- predict_by_group(truth, newdata = data, weighted = weighted)
  group_names <- colnames(mu_true_grp)
  
  rows <- lapply(names(fits), function(nm) {
    mu_hat_grp <- predict_by_group(fits[[nm]], newdata = data, weighted = weighted)
    
    grp_rows <- lapply(seq_along(group_names), function(g) {
      dplyr::tibble(
        model = nm,
        group = group_names[g],
        rmse = metric_rmse(mu_true_grp[, g], mu_hat_grp[, g]),
        cor = cor(mu_true_grp[, g], mu_hat_grp[, g])
      )
    })
    dplyr::bind_rows(grp_rows)
  })
  
  dplyr::bind_rows(rows)
}
