# Helpers for the binomial (classification) cover models.
#
# Each @examples block assigns every argument, so you can run it and then step
# through the function body line by line.


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


#' Predicted probability from a fitted cover classification model
#'
#' Dispatches on the engine the model was fit with, so callers don't need to
#' know which one it was.
#'
#' @param fit Object of class "cover_classification".
#' @param x Design matrix from `prepare_predictors()`, with the same columns,
#'   in the same order, as the matrix the model was fit to.
#' @return Numeric vector of predicted probabilities, one per row of `x`.
#' @examples
#' fit <- read_cover_model("classification", "forest", "c01", "m01")
#' dat <- read_cover_training("c01")
#' spec <- fit$config$spec
#' x <- prepare_predictors(head(dat, 10), pred_vars = spec$pred_vars,
#'                         scale_df = fit$scale_df, squares = spec$squares,
#'                         interactions = spec$interactions)
#' predict_prob(fit = fit, x = x)
predict_prob <- function(fit, x) {
  engine <- fit$config$spec$engine
  
  switch(
    engine,
    glmnet = {
      requireNamespace("glmnet")  # to get the predict() method
      as.numeric(predict(fit$fit, newx = x, s = fit$lambda,
                         type = "response"))
    },
    ranger = {
      requireNamespace("ranger")  # to get the predict() method
      predict(fit$fit, data = x)$predictions[, "1"]
    },
    stop("unknown engine: ", engine)
  )
}




#' Predict a cover classification model onto a raster
#'
#' Cells with all predictors present are converted to a data frame, run
#' through `prepare_predictors()` with the model's stored settings (log1p,
#' fixed scaling, squares and interactions), and predicted in chunks of rows
#' to limit memory.
#'
#' @param fit Object of class "cover_classification" (from
#'   `read_cover_model()`).
#' @param rast `SpatRaster` of predictors in original units, with a layer for
#'   each source variable (e.g. `MAP` for `log1p_MAP`).
#' @param chunk_size Number of cells predicted at a time.
#' @param ... Not used.
#' @return Two-layer `SpatRaster`: `prob` (predicted probability) and
#'   `class` (1 where prob > the model's threshold, else 0).
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
  
  prob <- map(chunks, \(i) {
    x <- prepare_predictors(df[i, ],
                            pred_vars = spec$pred_vars,
                            scale_df = fit$scale_df,
                            squares = spec$squares,
                            interactions = spec$interactions)
    # same columns, in the same order, as the fitted design matrix
    stopifnot(identical(colnames(x), fit$config$x_colnames))
    predict_prob(fit, x)
  }) |>
    unlist(use.names = FALSE)
  
  r_prob <- terra::rast(rast, nlyrs = 1)
  r_prob[df$cell] <- prob
  
  r_class <- r_prob > fit$threshold
  
  out <- c(r_prob, r_class)
  names(out) <- c("prob", "class")
  out
}


#' Partial dependence for a cover classification model
#'
#' For each predictor (in original units, e.g. MAP for log1p_MAP), sets it to
#' each value on a grid for every row of a background sample, predicts, and
#' averages the predicted probability. Other predictors keep their observed
#' values.
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
                            interactions = spec$interactions)
    mean(predict_prob(fit, x))
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


# internal: predicted probability for data in original units
.predict_newdata <- function(fit, newdat) {
  spec <- fit$config$spec
  x <- prepare_predictors(newdat, pred_vars = spec$pred_vars,
                          scale_df = fit$scale_df,
                          squares = spec$squares,
                          interactions = spec$interactions)
  predict_prob(fit, x)
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
#' unchanged), and averages the change in predicted probability within each
#' bin. The running sum of those averages, centered to mean zero over the
#' data, is the ALE curve (Apley & Zhu 2020). Unlike a PDP, the model is only
#' evaluated near combinations of predictors that occur in the data.
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
#'   ALE, probability scale) and `mean_pred` (mean predicted probability over
#'   the data), for `plot_ale()`.
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
#'   `x_value`, `yhat` and `mean_pred` (mean predicted probability within the
#'   subset), for `plot_ale_filtered()`.
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
    facet_wrap(~ variable, scales = "free_x") +
    labs(x = NULL, y = ylab)
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
