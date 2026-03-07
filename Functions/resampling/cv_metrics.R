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

#' Select lambda from a scored path
#'
#' Applies a named selection rule to a scored lambda path.
#'
#' @param score_df Data frame containing lambda values and metric columns.
#' @param metric Character scalar naming the metric column to optimize.
#' @param rule Character scalar. Currently `"min"` or `"largest_within_tol"`.
#' @param tol Numeric scalar >= 0. Used only for `"largest_within_tol"`.
#'
#' @return One-row data frame corresponding to the selected lambda.
select_lambda <- function(score_df,
                          metric = "mae_log",
                          rule = c("min", "largest_within_tol"),
                          tol = 0) {
  rule <- match.arg(rule)
  
  if (rule == "min") {
    return(select_lambda_min_metric(score_df = score_df, metric = metric))
  }
  
  select_lambda_within_tol(
    score_df = score_df,
    metric = metric,
    tol = tol
  )
}
