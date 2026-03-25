# Partial dependence plot functions for cwexp model
# Goes in Functions/viz/pdp_plots.R


#' Compute partial dependence for a cwexp model
#'
#' For each focal variable, creates a grid of values (percentiles of the
#' observed data plus min/max), then for each grid value substitutes it into
#' a background subsample, predicts, and averages. This gives the marginal
#' effect of the focal variable averaging over the joint distribution of all
#' other variables.
#'
#' @param fit Fitted cwexp model object.
#' @param data Data frame with predictor and cover columns.
#' @param focal_vars Character vector of predictor variables to compute PDPs
#'   for. Default: all main-effect predictors in the formula.
#' @param n_background Integer; number of background rows to subsample.
#'   Default 1000.
#' @param n_grid Integer; number of percentile grid points (min and max are
#'   always added). Default 100.
#' @param type Character; `"total"` for total predicted biomass, `"by_group"`
#'   for per-PFT predictions.
#' @param weighted Logical; only used when `type = "by_group"`. If TRUE,
#'   returns cover-weighted per-PFT biomass. If FALSE, returns potential
#'   biomass (cover = 100%).
#' @param seed Integer; random seed for background subsample.
#'
#' @return A tibble. For `type = "total"`: columns `variable`, `x_value`,
#'   `yhat`. For `type = "by_group"`: columns `variable`, `x_value`, `PFT`,
#'   `yhat`.
#' @examples
#' # note this is setup up so can steup throught the function
#' # (i.e. creating all args as objects)
#' dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#' set.seed(1)
#' dummy <- cwexp_make_dummy_data(n = 500)
#' data <- dummy$data
#' focal_vars = c("tmean", "ppt")
#' n_background = 200
#' n_grid = 20
#' type = 'total'
#' seed = 42
#' weighted = TRUE
#' fit <- cwexp_fit_tmb(
#'   data = data,
#'   formula = formula,
#'   cover_cols = cover_cols,
#'   dll = dll_en,
#'   penalty = "elastic_net",
#'   en_alpha = 0.5,
#'   lambda = 0.01
#' )
#'
#' # total biomass PDP
#' pdp_total <- cwexp_pdp(
#'   fit = fit, data = data,
#'   focal_vars = focal_vars,
#'   n_background = n_background, n_grid = n_grid,
#'   seed = seed,
#'   weighted = weighted
#' )
#' plot_cwexp_pdp(pdp_total)
#'
#' # per-PFT PDP (potential biomass)
#' pdp_pft <- cwexp_pdp(
#'   fit = fit, data = data,
#'   focal_vars = c("tmean", "ppt"),
#'   n_background = 200, n_grid = 20,
#'   type = "by_group", weighted = FALSE
#' )
#' plot_cwexp_pdp(pdp_pft)  
#' @export
cwexp_pdp <- function(fit, data, focal_vars = NULL,
                      n_background = 1000,
                      n_grid = 50,
                      type = c("total", "by_group"),
                      weighted = TRUE,
                      seed = 42) {

  type <- match.arg(type)
  
  formula <- fit$spec$formula
  cover_cols <- fit$spec$cover_cols
  
  stopifnot(
    !is.null(formula),
    !is.null(cover_cols)
  )

  # default focal vars: main-effect predictors from formula

  if (is.null(focal_vars)) {
    focal_vars <- attr(stats::terms(formula), "term.labels")
    # drop interaction terms (contain ":")
    focal_vars <- focal_vars[!grepl(":", focal_vars)]
  }

  stopifnot(all(focal_vars %in% names(data)))

  # subsample background data
  set.seed(seed)
  n <- nrow(data)
  bg_idx <- if (n <= n_background) seq_len(n) else sample.int(n, n_background)
  bg_data <- data[bg_idx, , drop = FALSE]

  pft_labels <- stringr::str_remove(cover_cols, "Cov$")

  # compute PDP for each focal variable
  results <- lapply(focal_vars, function(var) {

    # grid: percentiles + min/max
    x_obs <- data[[var]]
    x_grid <- seq(min(x_obs, na.rm = TRUE), 
                 max(x_obs, na.rm = TRUE), length.out = n_grid)

    # for each grid value, substitute into background and predict
    grid_results <- lapply(x_grid, function(x_val) {
      bg_mod <- bg_data
      bg_mod[[var]] <- x_val

      if (type == "total") {
        preds <- predict(fit, bg_mod, type = "mu")
        dplyr::tibble(
          variable = var,
          x_value = x_val,
          yhat = mean(preds)
        )
      } else {
        prep <- cwexp_prepare(bg_mod, formula, cover_cols)
        mu_grp <- cwexp_mu_by_group(
          alpha = fit$par$alpha,
          B = fit$par$B,
          X = prep$X,
          C = prep$C,
          weighted = weighted
        )
        means <- colMeans(mu_grp)
        dplyr::tibble(
          variable = var,
          x_value = x_val,
          PFT = pft_labels,
          yhat = unname(means)
        )
      }
    })

    dplyr::bind_rows(grid_results)
  })

  out <- dplyr::bind_rows(results)

  out
}


#' Plot partial dependence for a cwexp model
#'
#' @param pdp_df Output of `cwexp_pdp()`.
#' @param title Character; plot title.
#' @param ylab Character; y-axis label.
#' @param line_color Character; line color (used when no PFT grouping).
#' @param line_size Numeric; line width.
#'
#' @return A ggplot object. For total PDP: faceted by variable. For per-PFT
#'   PDP: `facet_grid(PFT ~ variable)` with variable labels on the bottom.
#' @export
plot_cwexp_pdp <- function(pdp_df,
                           title = NULL,
                           ylab = "Predicted biomass",
                           line_color = "steelblue",
                           line_size = 0.8) {

  has_pft <- "PFT" %in% names(pdp_df)

  if (has_pft) {
    g <- ggplot2::ggplot(pdp_df, ggplot2::aes(x = x_value, y = yhat)) +
      ggplot2::geom_line(color = line_color, linewidth = line_size) +
      ggplot2::facet_grid(PFT ~ variable, scales = "free", switch = 'x') +
      ggplot2::labs(
        x = NULL,
        y = ylab,
        title = title
      ) +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(angle = 0),
        strip.placement = "outside",
        strip.background.x = ggplot2::element_blank()
      )
  } else {
    g <- ggplot2::ggplot(pdp_df, ggplot2::aes(x = x_value, y = yhat)) +
      ggplot2::geom_line(color = line_color, linewidth = line_size) +
      ggplot2::facet_wrap(~ variable, scales = "free_x", 
                          strip.position = 'bottom') +
      ggplot2::labs(
        x = NULL,
        y = ylab,
        title = title
      ) +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background = ggplot2::element_blank()
      )
  }

  g 
}
