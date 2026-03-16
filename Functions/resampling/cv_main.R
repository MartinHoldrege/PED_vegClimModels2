# main functions for cross validation

#' Make a lambda sequence for one training dataset
#'
#' This is a lightweight wrapper that:
#' 1. computes `lambda_max_l1` on `data`, then
#' 2. builds a lambda path using `lambda_path_fun`.
#'
#' By default, `lambda_path_fun = cwexp_lambda_path`.
#'
#' @param data Data frame used for modeling.
#' @param en_alpha Elastic-net mixing parameter.
#' @param lambda_max_fun Function that computes `lambda_max_l1` on `data`.
#' (see example)
#' @param lambda_max_args Named list of additional arguments passed to
#'   `lambda_max_fun`.
#' @param lambda_path_fun Function that creates a lambda path.
#' @param lambda_path_args Named list of additional arguments passed to
#'   `lambda_path_fun`.
#'
#' @return Numeric lambda sequence.
#' @examples
#' # Simulate example data
#' dummy <- cwexp_make_dummy_data(n = 100)
#' dat <- dummy$data
#'
#' formula <- totalBio ~ tmean + ppt + vpd + sand
#' cover_cols <- dummy$spec$cover_cols
#'
#' # Compile elastic-net TMB model
#' dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#'
#' # Generate lambda path
#' lam_seq <- make_lambda_seq(
#'   data = dat,
#'   en_alpha = 0.5,
#'   lambda_max_fun = cwexp_lambda_max_l1_tmb,
#'   lambda_max_args = list(
#'     formula = formula,
#'     cover_cols = cover_cols,
#'     dll_en = dll_en
#'   ),
#'   lambda_path_fun = cwexp_lambda_path,
#'   lambda_path_args = list(
#'     n_lambda = 5,
#'     include_zero = TRUE
#'   )
#' )
#'
#' lam_seq
make_lambda_seq <- function(data,
                            en_alpha,
                            lambda_max_fun,
                            lambda_max_args = list(),
                            lambda_path_fun = cwexp_lambda_path,
                            lambda_path_args = list()) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.function(lambda_max_fun)) {
    stop("lambda_max_fun must be a function.")
  }
  if (!is.function(lambda_path_fun)) {
    stop("lambda_path_fun must be a function.")
  }
  if (!is.finite(en_alpha) || en_alpha <= 0 || en_alpha > 1) {
    stop("en_alpha must be in (0, 1].")
  }
  
  lambda_max_obj <- do.call(
    lambda_max_fun,
    c(list(data = data), lambda_max_args)
  )
  
  if (is.null(lambda_max_obj$lambda_max_l1)) {
    stop("lambda_max_fun output must contain 'lambda_max_l1'.")
  }
  
  lambda_seq <- do.call(
    lambda_path_fun,
    c(
      list(
        lambda_max_l1 = lambda_max_obj$lambda_max_l1,
        en_alpha = en_alpha
      ),
      lambda_path_args
    )
  )
  
  if (!is.numeric(lambda_seq) || length(lambda_seq) == 0L || any(!is.finite(lambda_seq))) {
    stop("lambda_path_fun must return a non-empty numeric vector of finite values.")
  }
  
  lambda_seq
}


#' Run path fitting, scoring, and lambda selection for one fold
#'
#' This function runs the full train/test workflow for one fold.
#' If `lambda_seq` is `NULL`, it computes `lambda_max_l1` on training data,
#' builds a lambda path, fits the path on training data, scores the path on
#' test data, and selects a best lambda.
#' If `lambda_seq` is supplied, that sequence is used directly.
#'
#' @param data Data frame used for modeling.
#' @param fold Fold object containing row indices, such as from
#'   `make_cluster_folds()`.
#' @param en_alpha Elastic-net mixing parameter.
#' @param lambda_seq Optional numeric lambda sequence. If supplied, this is used
#'   directly and no fold-specific lambda path is computed.
#' @param lambda_max_fun Function that computes `lambda_max_l1` on training data.
#'   Only used when `lambda_seq` is `NULL`.
#' @param lambda_max_args Named list of additional arguments passed to
#'   `lambda_max_fun`.
#' @param lambda_path_fun Function that creates a lambda path. Only used when
#'   `lambda_seq` is `NULL`.
#' @param lambda_path_args Named list of additional arguments passed to
#'   `lambda_path_fun`.
#' @param fit_path_fun Function that fits models across a lambda path.
#' @param fit_path_args Named list of additional arguments passed to
#'   `fit_path_fun`.
#' @param score_fun Function that scores a fitted path on test data.
#' @param score_args Named list of additional arguments passed to `score_fun`.
#' @param select_fun Function that selects the best lambda from scored results.
#' @param select_args Named list of additional arguments passed to `select_fun`.
#' @param keep_best_fit Logical. If `TRUE`, return the full fit object for the
#'   selected lambda.
#' @param keep_path_fit Logical. If `TRUE`, return the full path-fit object.
#'
#' @return A list containing fold metadata, scored results, selected lambda,
#'   lightweight parameter objects for all lambdas, and optionally the best
#'   full fit and/or the full path-fit object.
#' @examples
#'
#' # simulate small dataset
#' set.seed(1)
#' dummy <- cwexp_make_dummy_data(n = 1000)
#' dat <- dummy$data
#' cover_cols <- dummy$spec$cover_cols
#'
#' # environmental clusters
#' env_obj <- make_env_clusters(
#'   data = dat,
#'   vars = c("tmean", "ppt", "vpd", "sand"),
#'   k = 6,
#'   seed = 1
#' )
#'
#' # create cluster-based folds
#' folds <- make_cluster_folds(
#'   env_cluster = env_obj$env_cluster,
#'   n_folds = 2,
#'   seed = 1
#' )
#'
#' # compile model
#' dll <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#'
#' # run workflow for one fold
#' res <- run_fold_path(
#'   data = dat,
#'   fold = folds[[1]],
#'   en_alpha = 0.5,
#'
#'   lambda_max_fun = cwexp_lambda_max_l1_tmb,
#'   lambda_max_args = list(
#'     formula = totalBio ~ tmean + ppt + vpd + sand,
#'     cover_cols = cover_cols,
#'     dll_en = dll
#'   ),
#'
#'   lambda_path_args = list(
#'     n_lambda = 5,
#'     include_zero = TRUE
#'   ),
#'
#'   fit_path_fun = cwexp_fit_lambda_path_tmb,
#'   fit_path_args = list(
#'     formula = totalBio ~ tmean + ppt + vpd + sand,
#'     cover_cols = cover_cols,
#'     dll = dll,
#'     include_report = FALSE
#'   )
#' )
#'
#' # inspect scored lambda path
#' res$scores
#'
#' # selected lambda
#' res$selected
#'
#' # parameters for all lambdas
#' res$path_par
#'
#' # best fitted model
#' res$best_fit
run_fold_path <- function(data,
                          fold,
                          en_alpha,
                          lambda_seq = NULL,
                          lambda_max_fun = NULL,
                          lambda_max_args = list(),
                          lambda_path_fun = cwexp_lambda_path,
                          lambda_path_args = list(),
                          fit_path_fun,
                          fit_path_args = list(),
                          score_fun = score_path,
                          score_args = list(),
                          select_fun = select_lambda,
                          select_args = list(
                            metric = "mae_log",
                            rule = "min",
                            tol = 0
                          ),
                          keep_best_fit = TRUE,
                          keep_path_fit = FALSE) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.list(fold) || is.null(fold$train_rows) || is.null(fold$test_rows)) {
    stop("fold must be a list with elements 'train_rows' and 'test_rows'.")
  }
  if (!is.function(lambda_max_fun)) {
    stop("lambda_max_fun must be a function.")
  }
  if (!is.function(lambda_path_fun)) {
    stop("lambda_path_fun must be a function.")
  }
  if (!is.function(fit_path_fun)) {
    stop("fit_path_fun must be a function.")
  }
  if (!is.function(score_fun)) {
    stop("score_fun must be a function.")
  }
  if (!is.function(select_fun)) {
    stop("select_fun must be a function.")
  }
  if (!is.finite(en_alpha) || en_alpha <= 0 || en_alpha > 1) {
    stop("en_alpha must be in (0, 1].")
  }
  
  split_dat <- subset_fold_data(data = data, fold = fold)
  train_dat <- split_dat$train
  test_dat <- split_dat$test
  
  # start new
  lambda_max_l1 <- NA_real_
  
  if (is.null(lambda_seq)) {
    if (!is.function(lambda_max_fun)) {
      stop("lambda_max_fun must be a function when lambda_seq is NULL.")
    }
    if (!is.function(lambda_path_fun)) {
      stop("lambda_path_fun must be a function when lambda_seq is NULL.")
    }
    
    lambda_max_obj <- do.call(
      lambda_max_fun,
      c(list(data = train_dat), lambda_max_args)
    )
    
    if (is.null(lambda_max_obj$lambda_max_l1)) {
      stop("lambda_max_fun output must contain 'lambda_max_l1'.")
    }
    
    lambda_max_l1 <- lambda_max_obj$lambda_max_l1
    
    lambda_seq <- do.call(
      lambda_path_fun,
      c(
        list(
          lambda_max_l1 = lambda_max_l1,
          en_alpha = en_alpha
        ),
        lambda_path_args
      )
    )
  }
  
  if (!is.numeric(lambda_seq) || length(lambda_seq) == 0L || any(!is.finite(lambda_seq))) {
    stop("lambda_seq must be a non-empty numeric vector of finite values.")
  }

  path_fit <- do.call(
    fit_path_fun,
    c(
      list(
        data = train_dat,
        en_alpha = en_alpha,
        lambda_seq = lambda_seq
      ),
      fit_path_args
    )
  )
  
  score_df <- do.call(
    score_fun,
    c(
      list(
        path_fit = path_fit,
        newdata = test_dat
      ),
      score_args
    )
  )
  
  selected <- do.call(
    select_fun,
    c(
      list(score_df = score_df),
      select_args
    )
  )
  
  # lightweight parameter objects for every lambda
  path_par <- extract_path_par_objects(path_fit)
  
  # identify selected model
  best_idx <- match(selected$lambda[[1]], path_fit$summary$lambda)
  
  spec = path_fit$fits[[1]]$spec
  spec$lambda <- NA # generic specification info
  out <- list(
    fold_id = fold$fold_id,
    # train_clusters = fold$train_clusters,
    # test_clusters = fold$test_clusters,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    en_alpha = en_alpha,
    lambda_max_l1 = lambda_max_l1,
    lambda_seq = lambda_seq,
    scores = score_df,
    selected = selected,
    path_par = path_par,
    spec = spec
  )
  
  if (keep_best_fit) {
    out$best_fit <- path_fit$fits[[best_idx]]
  }
  
  if (keep_path_fit) {
    out$path_fit <- path_fit
  }
  
  class(out) <- "fold_path_result"
  out
}


#' Run full inner cross-validation for one outer-training dataset
#'
#' Runs `run_fold_path()` across all inner folds using one shared lambda path.
#' If `lambda_seq` is `NULL`, the path is created once from `data` using
#' `make_lambda_seq()`.
#'
#' By default, this function returns a lightweight object that keeps only the
#' pieces typically needed for downstream plotting, inspection, and model
#' extraction:
#' - one shared `spec`
#' - per-fold score tables
#' - per-fold parameter objects for each lambda
#' - per-fold selected rows
#' - aggregated score summary across folds
#' - overall selected lambda
#'
#' @param data Data frame used for modeling. This should be the outer-training
#'   data for one outer CV split.
#' @param folds List of fold objects, such as from `make_cluster_folds()`.
#' @param en_alpha Elastic-net mixing parameter.
#' @param lambda_seq Optional numeric lambda sequence shared across all inner
#'   folds. If `NULL`, it is created once from `data`.
#' @param lambda_max_fun Function that computes `lambda_max_l1` on `data`.
#'   Only used when `lambda_seq` is `NULL`.
#' @param lambda_max_args Named list of additional arguments passed to
#'   `lambda_max_fun`.
#' @param lambda_path_fun Function that creates a lambda path. Only used when
#'   `lambda_seq` is `NULL`.
#' @param lambda_path_args Named list of additional arguments passed to
#'   `lambda_path_fun`.
#' @param fit_path_fun Function that fits models across a lambda path.
#' @param fit_path_args Named list of additional arguments passed to
#'   `fit_path_fun`.
#' @param score_fun Function that scores a fitted path on test data.
#' @param score_args Named list of additional arguments passed to `score_fun`.
#' @param select_fun Function that selects the best lambda from scored results.
#' @param select_args Named list of additional arguments passed to `select_fun`.
#' @param keep_fold_results Logical. If `TRUE`, keep the full list returned by
#'   `run_fold_path()` for each fold. Useful for debugging.
#' @param run_fold_path_args Named list of additional arguments passed to
#'   `run_fold_path()`.
#' @param parallel Logical. If `TRUE`, run inner folds in parallel via the
#'   `future` backend. Requires `future` and `future.apply`. Default `FALSE`.
#' @param dll_path Character; path to the `.cpp` file for the TMB model
#'   (e.g., `"src/cwexp_lognormal_en_tmb.cpp"`). Needed when `parallel = TRUE`
#'   with `plan(multisession)` (Windows), because each worker is a separate R
#'   session without the DLL loaded. Not needed with `plan(multicore)` (Linux),
#'   where forked workers inherit loaded DLLs from the main process.
#' @return A lightweight inner-CV result.
#' @examples
#' # Simulate example data
#' dummy <- cwexp_make_dummy_data(n = 300)
#' dat <- dummy$data
#'
#' # Model spec
#' formula <- totalBio ~ tmean + ppt + vpd + sand
#' cover_cols <- dummy$spec$cover_cols
#' clust <- make_env_clusters(
#'   data = dat,
#'   vars = c("tmean", "ppt", "vpd"),
#'   k = 5,
#'   seed = 1
#' )
#' env_cluster <- clust$env_cluster
#' # Compile elastic-net TMB model
#' dll_path <- "src/cwexp_lognormal_en_tmb.cpp"
#' dll_en <- cwexp_tmb_compile(dll_path, quiet = TRUE)
#'
#' # Create inner CV folds
#' folds <- make_cluster_folds(env_cluster, n_folds = 3)
#'
#' # Run inner cross-validation
#' res <- run_inner_cv(
#'   data = dat,
#'   folds = folds,
#'   en_alpha = 0.5,
#'   lambda_max_fun = cwexp_lambda_max_l1_tmb,
#'   lambda_max_args = list(
#'     formula = formula,
#'     cover_cols = cover_cols,
#'     dll_en = dll_en
#'   ),
#'   fit_path_fun = cwexp_fit_lambda_path_tmb,
#'   fit_path_args = list(
#'     formula = formula,
#'     cover_cols = cover_cols,
#'     dll = dll_en
#'   )
#' )
#'
#' # Selected lambda across folds
#' res$selected
#' # for running in parallel:
#'# library(future)
# # plan(multisession, workers = 4)  # number of folds or cores whichever is less
# # the call to run_inner_cv(..., parallel = TRUE, dll_path = dll_path)
# # plan(sequential)  # reset when done
run_inner_cv <- function(data,
                         folds,
                         en_alpha,
                         lambda_seq = NULL,
                         lambda_max_fun = NULL,
                         lambda_max_args = list(),
                         lambda_path_fun = cwexp_lambda_path,
                         lambda_path_args = list(),
                         fit_path_fun,
                         fit_path_args = list(),
                         score_fun = score_path,
                         score_args = list(),
                         select_fun = select_lambda,
                         select_args = list(
                           metric = "mae_log",
                           rule = "1se",
                           tol = 0
                         ),
                         keep_fold_results = FALSE,
                         run_fold_path_args = list(),
                         parallel = FALSE,
                         dll_path = NULL
) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.list(folds) || length(folds) == 0L) {
    stop("folds must be a non-empty list.")
  }
  if (!is.finite(en_alpha) || en_alpha <= 0 || en_alpha > 1) {
    stop("en_alpha must be in (0, 1].")
  }
  if (!is.list(lambda_max_args)) {
    stop("lambda_max_args must be a list.")
  }
  if (!is.list(lambda_path_args)) {
    stop("lambda_path_args must be a list.")
  }
  if (!is.function(fit_path_fun)) {
    stop("fit_path_fun must be a function.")
  }
  if (!is.list(fit_path_args)) {
    stop("fit_path_args must be a list.")
  }
  if (!is.function(score_fun)) {
    stop("score_fun must be a function.")
  }
  if (!is.list(score_args)) {
    stop("score_args must be a list.")
  }
  if (!is.function(select_fun)) {
    stop("select_fun must be a function.")
  }
  if (!is.list(select_args)) {
    stop("select_args must be a list.")
  }
  if (!is.list(run_fold_path_args)) {
    stop("run_fold_path_args must be a list.")
  }
  
  if (isTRUE(parallel)) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages 'future' and 'future.apply' are required for parallel = TRUE.")
    }
  }
  
  # lambda sequence (computed ONCE, before any parallelism) ---
  if (is.null(lambda_seq)) {
    if (!is.function(lambda_max_fun)) {
      stop("lambda_max_fun must be a function when lambda_seq is NULL.")
    }
    if (!is.function(lambda_path_fun)) {
      stop("lambda_path_fun must be a function when lambda_seq is NULL.")
    }
    
    lambda_seq <- make_lambda_seq(
      data = data,
      en_alpha = en_alpha,
      lambda_max_fun = lambda_max_fun,
      lambda_max_args = lambda_max_args,
      lambda_path_fun = lambda_path_fun,
      lambda_path_args = lambda_path_args
    )
  }
  
  if (!is.numeric(lambda_seq) || length(lambda_seq) == 0L ||
      any(!is.finite(lambda_seq))) {
    stop("lambda_seq must be a non-empty numeric vector of finite values.")
  }
  
  .run_one_fold <- function(fold) {
    
    # ensure DLL is loaded on parallel workers
    if (isTRUE(parallel) && !is.null(dll_path)) {
      cwexp_tmb_compile(dll_path, quiet = TRUE)
    }
    
    if (isTRUE(parallel)) {
      # loading the needed functions (b/ own environment)
      source_functions(path = file.path(
            './Functions',
            c('grouping', 'models', 'resampling')
          )
          )
    }
    
    do.call(
      run_fold_path,
      c(
        list(
          data = data,
          fold = fold,
          en_alpha = en_alpha,
          lambda_seq = lambda_seq,
          lambda_max_fun = lambda_max_fun,
          lambda_max_args = lambda_max_args,
          lambda_path_fun = lambda_path_fun,
          lambda_path_args = lambda_path_args,
          fit_path_fun = fit_path_fun,
          fit_path_args = fit_path_args,
          score_fun = score_fun,
          score_args = score_args,
          select_fun = select_fun,
          select_args = select_args
        ),
        run_fold_path_args
      )
    )
  }
  
  
  if (isTRUE(parallel)) {
    fold_results <- future.apply::future_lapply(
      folds, .run_one_fold, future.seed = TRUE
    )
  } else {
    fold_results <- lapply(folds, .run_one_fold)
  }
  names(fold_results) <- paste0("fold", seq_along(folds))
  
  spec_list <- lapply(fold_results, \(x) x[['spec']])
  
  spec_norm <- lapply(spec_list, function(x) {
    x$lambda <- NULL
    x
  })
  
  spec_ref <- spec_norm[[1]]
  
  same_spec <- all(purrr::map_lgl(spec_norm, \(x) identical(x, spec_ref)))
  
  if (!same_spec) {
    stop("spec differs across folds after removing lambda.")
  }
  
  scores <- do.call(
    rbind,
    lapply(fold_results, function(x) {
      out <- x$scores
      out$fold_id <- x$fold_id
      out
    })
  )
  
  selected_by_fold <- do.call(
    rbind,
    lapply(fold_results, function(x) {
      out <- x$selected
      out$fold_id <- x$fold_id
      out
    })
  )
  
  # data frame with a list column of parameters
  path_par <- dplyr::bind_rows(
    lapply(fold_results, function(x) {
      out <- dplyr::tibble(
        fold_id = x$fold_id,
        lambda = purrr::map_dbl(x$path_par, \(l) l[["lambda"]]),
        en_alpha = purrr::map_dbl(x$path_par, \(l) l[["en_alpha"]])
      )
      # this is a list column
      out$par <- lapply(x$path_par, \(y) y$par)
      out
    })
  )


  numeric_cols <- names(scores)[vapply(scores, is.numeric, logical(1))]
  metric_cols <- setdiff(numeric_cols, c("lambda", "en_alpha", "fold_id"))
  
  if (length(metric_cols) == 0L) {
    stop("No numeric score columns found to summarize across folds.")
  }
  
  score_summary <- summarize_scores(scores, metric_cols)

  score_summary$en_alpha <- en_alpha

  # select the lambda that on average does best across folds
  selected <- do.call(
    select_fun,
    c(list(score_df = score_summary), select_args)
  )
  
  out <- list(
    en_alpha = en_alpha,
    lambda_seq = lambda_seq,
    spec = spec_ref,
    scores = scores,
    selected_by_fold = selected_by_fold,
    path_par = path_par,
    score_summary = score_summary,
    selected = selected
  )
  
  if (keep_fold_results) {
    out$fold_results <- fold_results
  }
  
  class(out) <- "inner_cv_result"
  out
}

if(FALSE) {
  # creating objects for all args
  # of run_inner_cv so i can step through

  # first run example code above that creates data
  # and folds
  # --- arguments for run_inner_cv ---
  
  data <- dat
  en_alpha <- 0.5
  
  lambda_seq <- NULL
  
  lambda_max_fun <- cwexp_lambda_max_l1_tmb
  
  lambda_max_args <- list(
    formula = formula,
    cover_cols = cover_cols,
    dll_en = dll_en
  )
  
  lambda_path_fun <- cwexp_lambda_path
  
  lambda_path_args <- list()
  
  fit_path_fun <- cwexp_fit_lambda_path_tmb
  
  fit_path_args <- list(
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en
  )
  
  score_fun <- score_path
  score_args <- list()
  
  select_fun <- select_lambda
  select_args <- list(
    metric = "mae_log",
    rule = "min",
    tol = 0
  )
  
  keep_fold_results <- FALSE
  
  run_fold_path_args <- list()
}



#' Select warm-start parameters from an inner CV result
#'
#' Chooses starting parameters for a final single-lambda fit using the stored
#' parameter objects in an `inner_cv_result`.
#'
#' Selection rule:
#' 1. Find the first model whose `lambda` matches the selected overall lambda.
#' 2. If none match, use the first available parameter object.
#'
#' @param inner_cv_result Object returned by `run_inner_cv()`.
#'
#' @return A named list of starting values with elements `alpha`, `B`,
#'   and `log_sigma`.
select_warm_start_par <- function(inner_cv_result) {
  if (!is.list(inner_cv_result)) {
    stop("inner_cv_result must be a list.")
  }
  if (is.null(inner_cv_result$path_par)) {
    stop("inner_cv_result must contain 'path_par'.")
  }
  if (is.null(inner_cv_result$selected)) {
    stop("inner_cv_result must contain 'selected'.")
  }
  
  path_par <- inner_cv_result$path_par
  selected_lambda <- inner_cv_result$selected$lambda[[1]]
  
  if (!is.data.frame(path_par) || !"par" %in% names(path_par) || !"lambda" %in% names(path_par)) {
    stop("inner_cv_result$path_par must be a data frame with columns 'lambda' and 'par'.")
  }
  if (nrow(path_par) == 0L) {
    stop("inner_cv_result$path_par must contain at least one row.")
  }
  
  idx <- which(path_par$lambda == selected_lambda)
  
  if (length(idx) == 0L) {
    idx <- 1L
  } else {
    idx <- idx[[1]]
  }
  
  par0 <- path_par$par[[idx]]
  
  if (!is.list(par0) || !all(c("alpha", "B") %in% names(par0))) {
    stop("Selected parameter object must contain at least 'alpha' and 'B'.")
  }
  
  sigma0 <- if ("log_sigma" %in% names(par0)) {
    par0$log_sigma
  } else if ("sigma" %in% names(par0)) {
    log(par0$sigma)
  } else {
    stop("Selected parameter object must contain either 'sigma' or 'log_sigma'.")
  }
  
  list(
    alpha = par0$alpha,
    B = par0$B,
    log_sigma = sigma0
  )
}



#' Run full nested outer cross-validation
#'
#' For each outer fold, this function:
#' 1. splits data into outer-train and outer-test,
#' 2. creates inner folds on the outer-train data,
#' 3. runs `run_inner_cv()` to select lambda,
#' 4. refits a single final model on all outer-train data at the selected lambda,
#' 5. scores that final model on outer-test data.
#'
#' By default, each outer-fold result stores:
#' - outer fold id
#' - selected lambda from inner CV
#' - outer test scores
#' - lightweight final model object (`spec`, `par`, `diag`)
#' - full `inner_cv_result`
#' - outer/inner cluster assignments (cluster ids only, not row ids)
#'
#' The top-level result also includes stacked outer scores, a summary of outer
#' performance across folds, the selected lambda for each outer fold, and a
#' final production lambda defined as the median selected lambda across outer
#' folds.
#'
#' @param data Data frame used for modeling.
#' @param outer_folds List of outer fold objects.
#' @param inner_fold_fun Function used to create inner folds from the
#'   outer-training data.
#' @param inner_fold_args Named list of additional arguments passed to
#'   `inner_fold_fun`.
#' @param en_alpha Elastic-net mixing parameter.
#' @param lambda_seq Optional numeric lambda sequence to use for all outer
#'   folds. If `NULL`, each outer fold creates its own lambda sequence from the
#'   corresponding outer-training data.
#' @param lambda_max_fun Function that computes `lambda_max_l1`.
#'   Only used when `lambda_seq` is `NULL`.
#' @param lambda_max_args Named list of additional arguments passed to
#'   `lambda_max_fun`.
#' @param lambda_path_fun Function used to create a lambda path.
#' @param lambda_path_args Named list of additional arguments passed to
#'   `lambda_path_fun`.
#' @param fit_path_fun Function that fits models across a lambda path.
#' @param fit_path_args Named list of additional arguments passed to
#'   `fit_path_fun`.
#' @param score_fun Function that scores a fitted path on test data.
#' @param score_args Named list of additional arguments passed to `score_fun`.
#' @param select_fun Function that selects the best lambda from scored results.
#' @param select_args Named list of additional arguments passed to `select_fun`.
#' @param final_fit_fun Function used to fit the final single-lambda model on
#'   the full outer-training data.
#' @param final_fit_args Named list of additional arguments passed to
#'   `final_fit_fun`.
#' @param final_score_fun Function used to score the final fitted model on the
#'   outer-test data. Should return a one-row data frame.
#' @param final_score_args Named list of additional arguments passed to
#'   `final_score_fun`.
#' @param use_warm_start Logical. If `TRUE`, use `warm_start_fun(inner_cv)` to
#'   generate starting values for the final outer-train fit.
#' @param keep_inner_cv Logical. If `TRUE`, store the full `inner_cv_result`
#'   for each outer fold.
#'
#' @return A list of class `outer_cv_result`.
#' @examples
#' # Simulate example data
#' dummy <- cwexp_make_dummy_data(n = 1000)
#' dat <- dummy$data
#'
#' # Model inputs
#' formula <- totalBio ~ tmean + ppt + vpd + sand
#' cover_cols <- dummy$spec$cover_cols
#'
#' # Compile elastic-net TMB model
#' dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#'
#' # Create environmental clusters and outer folds
#' clust <- make_env_clusters(dat, c("tmean", "ppt", "vpd", "sand"), k = 30)
#' dat$env_cluster <- clust$env_cluster
#'
#' outer_folds <- make_cluster_folds(
#'   env_cluster = dat$env_cluster,
#'   n_folds = 3
#' )
#'
#' # Run nested outer cross-validation
#' res <- run_outer_cv(
#'   data = dat,
#'   outer_folds = outer_folds,
#'   lambda_max_args = list(
#'     formula = formula,
#'     cover_cols = cover_cols,
#'     dll_en = dll_en
#'   ),
#'   lambda_path_args = list(
#'     n_lambda = 6,
#'     include_zero = TRUE
#'   ),
#'   fit_path_args = list(
#'     formula = formula,
#'     cover_cols = cover_cols,
#'     dll = dll_en
#'   ),
#'   final_fit_args = list(
#'     formula = formula,
#'     cover_cols = cover_cols,
#'     dll = dll_en,
#'     include_report = FALSE
#'   ),
#'   final_score_args = list(
#'     metrics = c("mae_log", "rmse_log"),
#'     predict_type = "mu"
#'   )
#' )
#'
#' # Final production lambda
#' res$lambda_final
run_outer_cv <- function(data,
                         outer_folds,
                         inner_fold_fun = make_cluster_folds,
                         inner_fold_args = list(n_fold = 3),
                         en_alpha = 0.5,
                         lambda_seq = NULL,
                         lambda_max_fun = cwexp_lambda_max_l1_tmb,
                         lambda_max_args = list(),
                         lambda_path_fun = cwexp_lambda_path,
                         lambda_path_args = list(),
                         fit_path_fun = cwexp_fit_lambda_path_tmb,
                         fit_path_args = list(),
                         score_fun = score_path,
                         score_args = list(),
                         select_fun = select_lambda,
                         select_args = list(
                           metric = "mae_log",
                           rule = "min",
                           tol = 0
                         ),
                         final_fit_fun = cwexp_fit_tmb,
                         final_fit_args = list(),
                         final_score_fun = score_fit_df,
                         final_score_args = list(),
                         use_warm_start = TRUE,
                         keep_inner_cv = TRUE) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.list(outer_folds) || length(outer_folds) == 0L) {
    stop("outer_folds must be a non-empty list.")
  }
  if (!is.function(inner_fold_fun)) {
    stop("inner_fold_fun must be a function.")
  }
  if (!is.list(inner_fold_args)) {
    stop("inner_fold_args must be a list.")
  }
  if (!is.finite(en_alpha) || en_alpha <= 0 || en_alpha > 1) {
    stop("en_alpha must be in (0, 1].")
  }
  if (!is.null(lambda_seq)) {
    if (!is.numeric(lambda_seq) || length(lambda_seq) == 0L ||
        any(!is.finite(lambda_seq)) || any(lambda_seq < 0)) {
      stop("lambda_seq must be NULL or a non-empty numeric vector of finite values >= 0.")
    }
  }
  if (!is.null(lambda_max_fun) && !is.function(lambda_max_fun)) {
    stop("lambda_max_fun must be a function.")
  }
  if (!is.list(lambda_max_args)) {
    stop("lambda_max_args must be a list.")
  }
  if (!is.function(lambda_path_fun)) {
    stop("lambda_path_fun must be a function.")
  }
  if (!is.list(lambda_path_args)) {
    stop("lambda_path_args must be a list.")
  }
  if (!is.function(fit_path_fun)) {
    stop("fit_path_fun must be a function.")
  }
  if (!is.list(fit_path_args)) {
    stop("fit_path_args must be a list.")
  }
  if (!is.function(score_fun)) {
    stop("score_fun must be a function.")
  }
  if (!is.list(score_args)) {
    stop("score_args must be a list.")
  }
  if (!is.function(select_fun)) {
    stop("select_fun must be a function.")
  }
  if (!is.list(select_args)) {
    stop("select_args must be a list.")
  }
  if (!is.function(final_fit_fun)) {
    stop("final_fit_fun must be a function.")
  }
  if (!is.list(final_fit_args)) {
    stop("final_fit_args must be a list.")
  }
  if (!is.function(final_score_fun)) {
    stop("final_score_fun must be a function.")
  }
  if (!is.list(final_score_args)) {
    stop("final_score_args must be a list.")
  }

  
  outer_results <- vector("list", length(outer_folds))
  names(outer_results) <- paste0("outer_fold", seq_along(outer_folds))
  
  for (i in seq_along(outer_folds)) {
    outer_fold <- outer_folds[[i]]
    
    split_dat <- subset_fold_data(data = data, fold = outer_fold)
    outer_train <- split_dat$train
    outer_test <- split_dat$test
    
    inner_folds <- do.call(
      inner_fold_fun,
      c(list(data = outer_train), inner_fold_args)
    )

    inner_cv <- run_inner_cv(
      data = outer_train,
      folds = inner_folds,
      en_alpha = en_alpha,
      lambda_seq = lambda_seq,
      lambda_max_fun = lambda_max_fun,
      lambda_max_args = lambda_max_args,
      lambda_path_fun = lambda_path_fun,
      lambda_path_args = lambda_path_args,
      fit_path_fun = fit_path_fun,
      fit_path_args = fit_path_args,
      score_fun = score_fun,
      score_args = score_args,
      select_fun = select_fun,
      select_args = select_args
    )
    
    selected_lambda <- inner_cv$selected$lambda[[1]]
    
    start <- NULL
    if (isTRUE(use_warm_start)) {
      start <- select_warm_start_par(inner_cv)
    }
    
    final_fit <- do.call(
      final_fit_fun,
      c(
        list(
          data = outer_train,
          start = start,
          penalty = "elastic_net",
          en_alpha = en_alpha,
          lambda = selected_lambda
        ),
        final_fit_args
      )
    )
    
    outer_scores <- do.call(
      final_score_fun,
      c(
        list(
          fit = final_fit,
          newdata = outer_test
        ),
        final_score_args
      )
    )
    
    if (!is.data.frame(outer_scores) || nrow(outer_scores) != 1L) {
      stop("final_score_fun must return a one-row data frame.")
    }
    
    outer_scores$outer_fold_id <- outer_fold$fold_id
    outer_scores$selected_lambda <- selected_lambda
    outer_scores$en_alpha <- en_alpha
    
    inner_fold_cluster_map <- lapply(inner_folds, function(x) {
      list(
        fold_id = x$fold_id,
        train_clusters = x$train_clusters,
        test_clusters = x$test_clusters
      )
    })
    
    final_model <- list(
      spec = final_fit$spec,
      par = final_fit$par,
      opt = list(
        convergence = final_fit$tmb$opt$convergence,
        objective = unname(final_fit$tmb$opt$objective),
        iterations = final_fit$tmb$opt$iterations,
        evaluations = final_fit$tmb$opt$evaluations,
        message = final_fit$tmb$opt$message
      )
    )
    
    out_i <- list(
      outer_fold_id = outer_fold$fold_id,
      selected_lambda = selected_lambda,
      outer_scores = outer_scores,
      final_model = final_model,
      outer_train_clusters = outer_fold$train_clusters,
      outer_test_clusters = outer_fold$test_clusters,
      inner_fold_cluster_map = inner_fold_cluster_map
    )
    
    if (keep_inner_cv) {
      out_i$inner_cv <- inner_cv
    }
    
    outer_results[[i]] <- out_i
  }
  
  outer_scores <- purrr::map_dfr(outer_results, \(x) x$outer_scores)

  lambda_final <- stats::median(outer_scores$selected_lambda)
  
  out <- list(
    en_alpha = en_alpha,
    outer_results = outer_results,
    outer_scores = outer_scores,
    lambda_final = lambda_final
  )
  
  class(out) <- "outer_cv_result"
  out
}


if(FALSE) {
  # example setup to step through run_outer_cv() line by line

  # use objects created in @example above
  
  # --- arguments for run_outer_cv() ---
  
  data <- dat
  
  inner_fold_fun <- make_cluster_folds
  inner_fold_args <- list(
    n_fold = 3
  )
  
  en_alpha <- 0.5
  
  lambda_seq <- NULL
  
  lambda_max_fun <- cwexp_lambda_max_l1_tmb
  lambda_max_args <- list(
    formula = formula,
    cover_cols = cover_cols,
    dll_en = dll_en
  )
  
  lambda_path_fun <- cwexp_lambda_path
  lambda_path_args <- list(
    n_lambda = 6,
    include_zero = TRUE
  )
  
  fit_path_fun <- cwexp_fit_lambda_path_tmb
  fit_path_args <- list(
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en
  )
  
  score_fun <- score_path
  score_args <- list()
  
  select_fun <- select_lambda
  select_args <- list(
    metric = "mae_log",
    rule = "min",
    tol = 0
  )
  
  final_fit_fun <- cwexp_fit_tmb
  final_fit_args <- list(
    formula = formula,
    cover_cols = cover_cols,
    dll = dll_en,
    include_report = FALSE
  )
  
  final_score_fun <- score_fit_df
  final_score_args <- list(
    metrics = c("mae_log", "rmse_log"),
    predict_type = "mu"
  )
  
  use_warm_start <- TRUE
  keep_inner_cv <- TRUE
}
