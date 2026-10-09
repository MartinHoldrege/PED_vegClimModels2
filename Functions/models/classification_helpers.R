# Helpers for the classification cover models (forest / non-forest).
#
# Binomial models predict a probability of forest. Continuous models (family
# "gaussian", "poisson" or "quasibinomial") predict tree cover (%), which is
# cut into classes the same way; "score" below means whichever of the two a
# model predicts.
#
# Each @examples block assigns every argument, so you can run it and then step
# through the function body line by line.


#' Pixel-years to fit to
#'
#' Applies `spec$rows` to the training data. "all" (or NULL) keeps every
#' row; otherwise `rows` is a quoted condition on the columns of `dat`, in
#' original units, e.g. `quote(MAP < 700)`. Rows where the condition is NA
#' are dropped, as in `dplyr::filter()`.
#'
#' @param dat Cover training data (e.g. from `read_cover_training()`).
#' @param rows "all", NULL, or a quoted expression.
#' @return `dat`, filtered.
#' @examples
#' dat <- data.frame(MAP = c(300, 650, 900, NA), cov_tree = c(0, 5, 40, 2))
#' rows <- quote(MAP < 700)
#' select_rows(dat = dat, rows = rows)
#' select_rows(dat = dat, rows = "all")
select_rows <- function(dat, rows) {
  if (is.null(rows) || identical(rows, "all")) {
    return(dat)
  }
  if (!is.language(rows)) {
    stop("spec$rows must be \"all\", NULL or a quoted expression, e.g. ",
         "quote(MAP < 700)")
  }
  out <- dplyr::filter(dat, !!rows)
  if (nrow(out) == 0) {
    stop("no rows left after filtering on ", deparse1(rows))
  }
  out
}


#' Probability cutoff that turns predictions into classes
#'
#' Wraps `PresenceAbsence::optimal.thresholds()`. With method "PredPrev=Obs"
#' the cutoff is the one where the predicted fraction in the class equals the
#' observed fraction in the data used here. Other methods from that function
#' can be named instead.
#'
#' @param obs Observed class, 0 or 1.
#' @param pred Predicted probability, same length as `obs`.
#' @param method Method name, as in the `Method` column that
#'   `optimal.thresholds()` returns.
#' @return Numeric cutoff.
#' @examples
#' set.seed(1)
#' obs <- rbinom(500, size = 1, prob = 0.3)
#' pred <- plogis(rnorm(500, mean = -1 + obs))
#' method <- "PredPrev=Obs"
#' choose_threshold(obs = obs, pred = pred, method = method)
choose_threshold <- function(obs, pred, method) {
  stopifnot(length(obs) == length(pred),
            all(obs %in% c(0, 1)),
            all(pred >= 0 & pred <= 1))
  
  dat <- data.frame(ID = seq_along(obs), obs = obs, pred = pred)
  
  # opt.methods = method computes only the requested rule; without it every
  # rule is computed, and the ones we don't use (ReqSens, ReqSpec, Cost) warn
  # about their defaults
  thresholds <- PresenceAbsence::optimal.thresholds(
    DATA = dat,
    threshold = 200,       # number of candidate cutoffs tested
    obs.prev = mean(obs),
    opt.methods = method
  )
  
  stopifnot(method %in% thresholds$Method)
  thresholds[thresholds$Method == method, 2]
}


#' Cutoff that makes the predicted fraction of forest match the observed one
#'
#' The "PredPrev=Obs" rule for predictions that aren't probabilities (e.g.
#' predicted % cover): the quantile of the predictions at 1 - observed
#' prevalence, so that about as many rows are above the cutoff as are
#' observed forest.
#'
#' @param obs Observed class, 0 or 1.
#' @param score Predictions, same length as `obs`; any scale.
#' @return Numeric cutoff, on the scale of `score`.
#' @examples
#' set.seed(1)
#' obs <- rbinom(500, size = 1, prob = 0.3)
#' score <- exp(rnorm(500, mean = 1 + obs))
#' cutoff <- prevalence_cutoff(obs = obs, score = score)
#' c(observed = mean(obs), predicted = mean(score > cutoff))
prevalence_cutoff <- function(obs, score) {
  stopifnot(length(obs) == length(score), all(obs %in% c(0, 1)))
  unname(quantile(score, probs = 1 - mean(obs)))
}


#' Prediction from a fitted cover classification model
#'
#' Dispatches on the engine the model was fit with, so callers don't need to
#' know which one it was. Binomial models return the probability of forest;
#' continuous ones return predicted tree cover (%), back-transformed from the
#' scale the model was fit on (see `score_name()`).
#'
#' @param fit Object of class "cover_classification".
#' @param x Design matrix from `prepare_predictors()`, with the same columns,
#'   in the same order, as the matrix the model was fit to.
#' @return Numeric vector of predictions, one per row of `x`.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' dat <- read_cover_training("c01")
#' spec <- fit$config$spec
#' x <- prepare_predictors(head(dat, 10), pred_vars = spec$pred_vars,
#'                         scale_df = fit$scale_df, squares = spec$squares,
#'                         interactions = spec$interactions,
#'                         interact_log1p = isTRUE(spec$interact_log1p))
#' predict_score(fit = fit, x = x)
predict_score <- function(fit, x) {
  engine <- fit$config$spec$engine
  
  switch(
    engine,
    glmnet = {
      requireNamespace("glmnet")  # to get the predict() method
      p <- as.numeric(predict(fit$fit, newx = x, s = fit$lambda,
                              type = "response"))
      family <- fit$config$spec$family
      # models fit before the family field existed were all binomial
      if (is.null(family)) family <- "binomial"
      switch(
        family,
        binomial = p,
        # fit to cover / 100, so back to % cover
        quasibinomial = p * 100,
        # fit to cover or log1p(cover), so back to % cover
        gaussian = .untransform_cover(p, fit$config$spec$response_transform),
        stop("family not set up for glmnet: ", family)
      )
    },
    ranger = {
      requireNamespace("ranger")  # to get the predict() method
      p <- predict(fit$fit, data = x)$predictions
      # probability forest: probability of the class; regression forest
      # (family "gaussian"): back to % cover
      family <- fit$config$spec$family
      if (is.null(family) || family == "binomial") {
        p[, "1"]
      } else {
        .fit_scale_to_score(p, fit$config$spec)
      }
    },
    grpreg = {
      requireNamespace("grpreg")  # to get the predict() method
      # columns repeated once per group, as the model was fit
      x_groups <- x[, fit$fit$groups$term, drop = FALSE]
      mu <- predict(fit$fit$grpreg, x_groups, lambda = fit$lambda,
                    type = "response")
      .untransform_cover(as.vector(mu), fit$config$spec$response_transform)
    },
    gam = {
      requireNamespace("mgcv")  # to get the predict() method
      # discrete = FALSE: exact predictions (discretized ones are approximate)
      p <- predict(fit$fit$gam, newdata = as.data.frame(x),
                   type = "response", discrete = FALSE)
      .fit_scale_to_score(as.numeric(p), fit$config$spec)
    },
    stop("unknown engine: ", engine)
  )
}


#' What a fitted classification model predicts
#'
#' @param fit Object of class "cover_classification".
#' @return "prob" for binomial models (probability of forest), "cover" for
#'   continuous ones (predicted % cover). Also the name of the first layer
#'   `predict_raster()` returns.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' score_name(fit = fit)
score_name <- function(fit) {
  family <- fit$config$spec$family
  # models fit before the family field existed were all binomial
  if (is.null(family) || family == "binomial") "prob" else "cover"
}


#' Outcome labels for a binary classification
#'
#' @param class_label Name of the positive class (e.g. "forest", "zero tree").
#' @return Character vector: true positive, true negative, false positive,
#'   false negative, e.g. "forest, correct", "not forest, correct",
#'   "false forest", "missed forest".
#' @examples
#' outcome_labels("zero tree")
outcome_labels <- function(class_label = "forest") {
  c(paste0(class_label, ", correct"),
    paste0("not ", class_label, ", correct"),
    paste("false", class_label),
    paste("missed", class_label))
}


#' Label each row by observed vs predicted class
#'
#' @param df Data frame with an observed class column (0/1, 1 = the positive
#'   class) and a prediction column.
#' @param pred_col Name of the prediction column (e.g. "pred_oof").
#' @param threshold Cutoff on the prediction; above it is predicted positive.
#' @param class_label Name of the positive class, used in the labels.
#' @param obs_col Name of the observed class column.
#' @return `df` with `predicted` (0/1) and `outcome`, a factor with levels
#'   from `outcome_labels(class_label)`.
#' @examples
#' df <- tibble(observed = c(1, 0, 0, 1), pred_oof = c(0.8, 0.1, 0.7, 0.2))
#' pred_col <- "pred_oof"
#' threshold <- 0.5
#' class_label <- "forest"
#' obs_col <- "observed"
#' classify_outcome(df = df, pred_col = pred_col, threshold = threshold,
#'                  class_label = class_label, obs_col = obs_col)
classify_outcome <- function(df, pred_col, threshold, class_label = "forest",
                             obs_col = "observed") {
  stopifnot(all(c(obs_col, pred_col) %in% names(df)),
            all(df[[obs_col]] %in% c(0, 1)))
  labels <- outcome_labels(class_label)
  obs <- df[[obs_col]]
  df |>
    mutate(predicted = as.integer(.data[[pred_col]] > threshold),
           outcome = case_when(
             obs == 1 & predicted == 1 ~ labels[1],
             obs == 0 & predicted == 0 ~ labels[2],
             obs == 0 & predicted == 1 ~ labels[3],
             obs == 1 & predicted == 0 ~ labels[4]
           ),
           outcome = factor(outcome, levels = labels))
}


# internal: undo the transformation applied to cover before fitting
.untransform_cover <- function(z, how) {
  switch(how,
         identity = z,
         log1p = expm1(z),
         stop("unknown response_transform: ", how))
}


# internal: predictions on the scale the model was fit on (probability,
# cover / 100, or transformed cover) to what predict_score() returns
# (probability, or % cover)
.fit_scale_to_score <- function(p, spec) {
  switch(spec$family,
         binomial = p,
         quasibinomial = p * 100,
         gaussian = .untransform_cover(p, spec$response_transform),
         stop("family not set up: ", spec$family))
}




#' Predict a cover classification model onto a raster
#'
#' Cells with all predictors present are converted to a data frame, run
#' through `prepare_predictors()` with the model's stored settings (log1p,
#' fixed scaling, squares and interactions), and predicted in chunks of rows
#' to limit memory.
#'
#' The first layer is named by `score_name()`: `prob` for binomial models,
#' `cover` (predicted %) for continuous ones.
#'
#' @param fit Object of class "cover_classification" (from
#'   `read_cover_model()`).
#' @param rast `SpatRaster` of predictors in original units, with a layer for
#'   each source variable (e.g. `MAP` for `log1p_MAP`).
#' @param chunk_size Number of cells predicted at a time.
#' @param ... Not used.
#' @return Two-layer `SpatRaster`: `prob` or `cover` (the prediction) and
#'   `class` (1 where the prediction > the model's threshold, else 0).
#'   Continuous cover models (class "cover_continuous", from 03_fit_cover.R)
#'   have no threshold, so only the `cover` layer.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' rast <- read_climate_raster("current") |> align_raster(read_mask())
#' chunk_size <- 1e6
#' pred <- predict_raster(fit = fit, rast = rast, chunk_size = chunk_size)
predict_raster.cover_classification <- function(fit, rast, chunk_size = 1e6,
                                                 ...) {
  spec <- fit$config$spec
  source_vars <- unique(str_remove(spec$pred_vars, "^log1p_"))
  stopifnot(inherits(rast, "SpatRaster"), all(source_vars %in% names(rast)))
  
  # one row per cell with every predictor present
  df <- terra::as.data.frame(rast[[source_vars]], cells = TRUE, na.rm = TRUE)
  
  rows <- seq_len(nrow(df))
  chunks <- split(rows, ceiling(rows / chunk_size))
  
  score <- map(chunks, \(i) {
    x <- prepare_predictors(df[i, ],
                            pred_vars = spec$pred_vars,
                            scale_df = fit$scale_df,
                            squares = spec$squares,
                            interactions = spec$interactions,
                            interact_log1p = isTRUE(spec$interact_log1p))
    # same columns, in the same order, as the fitted design matrix
    stopifnot(identical(colnames(x), fit$config$x_colnames))
    predict_score(fit, x)
  }) |>
    unlist(use.names = FALSE)
  
  r_score <- terra::rast(rast, nlyrs = 1)
  r_score[df$cell] <- score
  names(r_score) <- score_name(fit)
  
  if (is.null(fit$threshold)) return(r_score)
  
  r_class <- r_score > fit$threshold
  
  out <- c(r_score, r_class)
  names(out) <- c(score_name(fit), "class")
  out
}

# continuous cover models predict the same way, minus the class layer
predict_raster.cover_continuous <- predict_raster.cover_classification


#' Partial dependence for a cover classification model
#'
#' For each predictor (in original units, e.g. MAP for log1p_MAP), sets it to
#' each value on a grid for every row of a background sample, predicts, and
#' averages the prediction (probability, or % cover for continuous models).
#' Other predictors keep their observed values.
#'
#' @param fit Object of class "cover_classification".
#' @param dat Training data in original units (the background).
#' @param n_grid Number of grid values per predictor (quantiles from quantile_range)
#' @param n_background Rows sampled from `dat` as the background.
#' @return Tibble with `variable`, `x_value` and `yhat`, for `plot_pdp()`.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' dat <- read_cover_training("c01")
#' n_grid <- 20
#' n_background <- 2000
#' quantile_range <- c(0, 1)
#' pdp <- pdp_classification(fit = fit, dat = dat, n_grid = n_grid,
#'                           n_background = n_background)
pdp_classification <- function(fit, dat, n_grid = 50, n_background = 2000,
                               quantile_range = c(0, 1)) {
  spec <- fit$config$spec
  source_vars <- unique(str_remove(spec$pred_vars, "^log1p_"))
  
  dat <- dat[complete.cases(dat[source_vars]), ]
  set.seed(1)
  bg <- dat[sample(nrow(dat), min(n_background, nrow(dat))), ]
  
  predict_mean <- function(newdat) {
    x <- prepare_predictors(newdat, pred_vars = spec$pred_vars,
                            scale_df = fit$scale_df,
                            squares = spec$squares,
                            interactions = spec$interactions,
                            interact_log1p = isTRUE(spec$interact_log1p))
    mean(predict_score(fit, x))
  }
  
  map(source_vars, \(v) {
    value_range <- quantile(dat[[v]], probs = quantile_range, names = FALSE)
    grid <- seq(value_range[1], value_range[2], length.out = n_grid)
    yhat <- map_dbl(grid, \(value) {
      bg_v <- bg
      bg_v[[v]] <- value
      predict_mean(bg_v)
    })
    tibble(variable = v, x_value = grid, yhat = yhat)
  }) |>
    bind_rows()
}


# internal: prediction (see predict_score()) for data in original units
.predict_newdata <- function(fit, newdat) {
  spec <- fit$config$spec
  x <- prepare_predictors(newdat, pred_vars = spec$pred_vars,
                          scale_df = fit$scale_df,
                          squares = spec$squares,
                          interactions = spec$interactions,
                          interact_log1p = isTRUE(spec$interact_log1p))
  predict_score(fit, x)
}


# internal: ALE for one predictor over the rows of dat (see
# ale_classification()). Returns NULL if the predictor has too few distinct
# values in dat to form two bins.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' dat <- read_cover_training("c01")
#' n_bins <- 20
#' n_per_bin <- 250
#' quantile_range <- c(0.01, 0.99)
#' v = 'MAT'
#' rows <- rowSums(is.na(dat[, names(dat) %in% c(fit$config$spec$response, fit$config$spec$pred_vars)]))
#' dat <- dat[rows == 0, ]
.ale_one <- function(fit, dat, v, n_bins, n_per_bin, quantile_range) {
  # edges from all rows; type = 1 makes them observed values, so every bin
  # holds at least one row
  probs <- seq(quantile_range[1], quantile_range[2], length.out = n_bins + 1)
  edges <- unique(quantile(dat[[v]], probs = probs, type = 1, names = FALSE))
  n_edges <- length(edges)
  if (n_edges < 3) return(NULL)
  
  in_range <- dat[dat[[v]] >= edges[1] & dat[[v]] <= edges[n_edges], ]
  in_range$bin <- findInterval(in_range[[v]], edges,
                               rightmost.closed = TRUE, all.inside = TRUE)
  

  set.seed(1)
  smp <- slice_sample(in_range, n = n_per_bin, by = bin)
  bin <- smp$bin
  smp$bin <- NULL
  n <- nrow(smp)
  
  # each row at its bin's lower edge, its upper edge, and as observed
  lower <- smp
  lower[[v]] <- edges[bin]
  upper <- smp
  upper[[v]] <- edges[bin + 1]

  pred_obs <- .predict_newdata(fit, smp)
  pred_upper <- .predict_newdata(fit, upper)
  pred_lower <- .predict_newdata(fit, lower)
  change <- pred_upper - pred_lower

  bins <- factor(bin, levels = seq_len(n_edges - 1))
  
  # mean change per bin, then running sum from the lowest edge
  delta <- tapply(change, bins, mean)
  stopifnot(!anyNA(delta))
  ale <- c(0, cumsum(delta))
  
  # center so the curve averages zero of the sampled data
  ale_mid <- (ale[-n_edges] + ale[-1]) / 2
  ale <- ale - mean(ale_mid)
  
  # mean predicted probability
  mean_pred <- mean(pred_obs)
  
  tibble(x_value = edges, yhat = ale, mean_pred = mean_pred)
}


#' Accumulated local effects (ALE) for a cover classification model
#'
#' For each predictor (in original units, e.g. MAP for log1p_MAP), splits the
#' data into quantile bins of that predictor, samples rows from each bin,
#' moves each sampled row to its bin's lower and upper edge (other predictors
#' unchanged), and averages the change in prediction (probability, or % cover
#' for continuous models) within each bin. The running sum of those averages,
#' centered to mean zero over the data, is the ALE curve (Apley & Zhu 2020).
#' Unlike a PDP, the model is only evaluated near combinations of predictors
#' that occur in the data.
#'
#' @param fit Object of class "cover_classification".
#' @param dat Training data in original units.
#' @param n_bins Number of quantile bins per predictor (fewer where values
#'   tie).
#' @param n_per_bin Rows sampled from each bin (all rows if the bin has fewer).
#'   `Inf` uses every row. The model is predicted three times per sampled row.
#' @param quantile_range Range of each predictor covered, as quantiles.
#'   Trimming the tails keeps a few extreme rows from setting the outermost
#'   bin edges; rows outside the range are left out for that predictor.
#' @return Tibble with `variable`, `x_value` (bin edges), `yhat` (centered
#'   ALE, on the prediction's scale) and `mean_pred` (mean prediction over the
#'   data), for `plot_ale()`.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' dat <- read_cover_training("c01")
#' n_bins <- 20
#' n_per_bin <- 250
#' quantile_range <- c(0.01, 0.99)
#' ale <- ale_classification(fit = fit, dat = dat, n_bins = n_bins,
#'                           n_per_bin = n_per_bin,
#'                           quantile_range = quantile_range)
ale_classification <- function(fit, dat, n_bins = 20, n_per_bin = 250,
                               quantile_range = c(0.01, 0.99)) {
  source_vars <- unique(str_remove(fit$config$spec$pred_vars, "^log1p_"))
  dat <- dat[complete.cases(dat[source_vars]), ]
  
  map(source_vars, \(v) {
    .ale_one(fit, dat, v, n_bins = n_bins, n_per_bin = n_per_bin,
             quantile_range = quantile_range) |>
      mutate(variable = v, .before = 1)
  }) |>
    bind_rows()
}


#' ALE within the low and high percentiles of each filter variable
#'
#' Analogue of the filtered quantile plots: for each filter variable, splits
#' the data into rows below its `low` and above its `high` percentile, and
#' computes the ALE of every predictor within each subset (as in
#' `ale_classification()`). Differences in shape between the two subsets show
#' interactions in the fitted model.
#'
#' @param fit Object of class "cover_classification".
#' @param dat Training data in original units.
#' @param filter_vars Variables to split by; defaults to the predictors.
#' @param low,high Percentile cut-offs (as proportions) for the two subsets.
#' @param n_bins,n_per_bin,quantile_range As in `ale_classification()`, applied
#'   within each subset.
#' @return Tibble with `filter_var`, `percentile_category`, `variable`,
#'   `x_value`, `yhat` and `mean_pred` (mean prediction within the subset),
#'   for `plot_ale_filtered()`.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' dat <- read_cover_training("c01")
#' filter_vars <- c("MAT", "MAP")
#' low <- 0.2
#' high <- 0.8
#' n_bins <- 20
#' n_per_bin <- 100
#' quantile_range <- c(0.01, 0.99)
#' ale_f <- ale_classification_filtered(fit = fit, dat = dat,
#'                                      filter_vars = filter_vars,
#'                                      low = low, high = high,
#'                                      n_bins = n_bins,
#'                                      n_per_bin = n_per_bin,
#'                                      quantile_range = quantile_range)
ale_classification_filtered <- function(fit, dat, filter_vars = NULL,
                                        low = 0.2, high = 0.8,
                                        n_bins = 20, n_per_bin = 100,
                                        quantile_range = c(0.01, 0.99)) {
  source_vars <- unique(str_remove(fit$config$spec$pred_vars, "^log1p_"))
  if (is.null(filter_vars)) filter_vars <- source_vars
  stopifnot(all(filter_vars %in% names(dat)), low < high)
  
  dat <- dat[complete.cases(dat[union(source_vars, filter_vars)]), ]
  categories <- c(paste0("<", low * 100, "th"), paste0(">", high * 100, "th"))
  
  map(filter_vars, \(f) {
    pct <- ecdf(dat[[f]])(dat[[f]]) # percentile
    subsets <- list(dat[pct < low, ], dat[pct > high, ]) |>
      set_names(categories)
    
    imap(subsets, \(d, category) {
      map(source_vars, \(v) {
        .ale_one(fit, d, v, n_bins = n_bins, n_per_bin = n_per_bin,
                 quantile_range = quantile_range) |>
          mutate(filter_var = f, percentile_category = category,
                 variable = v, .before = 1)
      })
    })
  }) |>
    list_flatten() |>
    list_flatten() |>
    bind_rows() |>
    mutate(filter_var = factor(filter_var, levels = filter_vars),
           percentile_category = factor(percentile_category,
                                        levels = categories),
           variable = factor(variable, levels = source_vars))
}


#' Plot accumulated local effects
#'
#' Points at the bin edges, joined by lines, so each segment is the mean
#' change in predicted probability across one bin.
#'
#' @param ale Tibble from `ale_classification()`.
#' @param ylab y-axis label.
#' @return A ggplot object, one panel per predictor.
#' @examples
#' ale <- tibble(variable = rep(c("MAT", "MAP"), each = 3),
#'               x_value = c(0, 5, 10, 200, 400, 800),
#'               yhat = c(-0.1, 0, 0.1, 0.05, 0, -0.05))
#' ylab <- "ALE (change in predicted probability)"
#' plot_ale(ale = ale, ylab = ylab)
plot_ale <- function(ale, ylab = "ALE (change in predicted probability)") {
  ggplot(ale, aes(x = x_value, y = yhat)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_line() +
    geom_point(size = 1) +
    facet_wrap(~ variable, scales = "free_x",
               strip.position = 'bottom') +
    labs(x = NULL, y = ylab)+
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_blank()
    )
}


#' Plot ALE within the low and high percentiles of each filter variable
#'
#' Grid like the filtered quantile plots: columns are predictors, rows are
#' filter variables. Each line is the subset's mean predicted probability
#' plus its ALE, so the two lines can be compared in level as well as shape.
#' Points are bin edges; each segment is the mean change across one bin.
#'
#' @param ale Tibble from `ale_classification_filtered()`.
#' @param ylab y-axis label.
#' @param title Plot title.
#' @return A ggplot object.
#' @examples
#' ale <- tibble(filter_var = "MAT", variable = "MAP",
#'               percentile_category = rep(c("<20th", ">80th"), each = 3),
#'               x_value = rep(c(200, 400, 800), 2),
#'               yhat = c(-0.1, 0, 0.1, -0.2, 0, 0.2),
#'               mean_pred = rep(c(0.3, 0.6), each = 3))
#' ylab <- "Predicted probability (subset mean + ALE)"
#' title <- NULL
#' plot_ale_filtered(ale = ale, ylab = ylab, title = title)
plot_ale_filtered <- function(ale,
                              ylab = "ALE (subset mean + ALE)",
                              title = NULL) {
  ggplot(ale, aes(x = x_value, y = mean_pred + yhat,
                  colour = percentile_category)) +
    geom_line() +
    geom_point(size = 0.8) +
    facet_grid(filter_var ~ variable, scales = "free_x") +
    scale_colour_manual(name = NULL, values = c("#f03b20", "#0570b0")) +
    labs(x = NULL, y = ylab, title = title,
         caption = paste0("Columns: predictor variable. Rows: filtering ",
                          "variable (ALE computed only from pixels in its ",
                          "lowest/highest percentiles).")) +
    theme(legend.position = "top")
}