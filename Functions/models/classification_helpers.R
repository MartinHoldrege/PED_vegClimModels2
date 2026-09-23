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
  requireNamespace('glmnet') # to get predict() method
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
    as.numeric(predict(fit$fit, newx = x, s = fit$lambda, type = "response"))
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
  requireNamespace('glmnet')
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
    mean(predict(fit$fit, newx = x, s = fit$lambda, type = "response"))
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

