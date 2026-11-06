# trend_anomaly_check.R
#
# Do cover models with climate anomaly predictors (e.g. m05.5) track change in
# cover within a pixel better than models without them (e.g. m01.5)?
#
# Uses the RAP cells of the training data that have cover in every year. In
# each cell, for observed cover and for each model's in-sample (pred_in) and
# out-of-fold (pred_oof) predictions:
#   trend   = slope of the least-squares line of cover on year (% cover/year)
#   anomaly = residuals from that line (detrended year-to-year change)
#
# Summary per model and prediction type (trend_summary_<label>.csv):
#   trend_cor        correlation of observed and predicted trends (across cells)
#   anom_cor_median  median across cells of the per-cell correlation of
#                    observed and predicted anomalies (0 in a cell where the
#                    predictions don't vary from their trend line)
#   n_cells_no_cor   cells left out of anom_cor_median because the observed
#                    series doesn't vary from its trend line
#   anom_cor_pooled  correlation of observed and predicted anomalies (all
#                    cell-years)
#   anom_sd_ratio    sd(predicted anomalies) / sd(observed anomalies)
#
# Inputs: training data (read_cover_training()) and fitted models
#   (read_cover_model()) for each model version in vmcs
# Outputs (out_dir), <label> = <cover model>_<vc>_<vmcs>:
#   trend_<label>.pdf, trend_maps_<label>.pdf, trend_summary_<label>.csv
#
# November 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
theme_set(theme_custom1())

# params ------------------------------------------------------------------

vc <-  "c07" # "c08" #                          # data version the models were fit to
cover_models <- "bare_ground" # c("shrub", "herbaceous") # 
# model versions to compare
# e.g. 1.5 doesn't use climate anomalies as predictors, 5.5 does
vmcs <-    c("m04.2", "m04.3") #  c("m01.5", "m05.5") #  
years <- 2010:2024                        # cells need cover in every year
pred_cols <- c("pred_in", "pred_oof")     # columns of fit$predictions
n_example_cells <- 12
out_dir <- "Figures/Cover/trend_check"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# functions ---------------------------------------------------------------

#' Slope of the least-squares line of y on x
#' @param x,y Numeric vectors.
#' @return The slope (units of y per unit of x).
slope <- function(x, y) cov(x, y) / var(x)

#' Residuals from the least-squares line of y on x
#' (i.e. y minus prediction of y)
#' @param x,y Numeric vectors.
#' @return Vector of residuals, the length of y.
detrend <- function(x, y) y - mean(y) - slope(x, y) * (x - mean(x))

#' Correlation of observed and predicted anomalies, for a series that may not
#' vary from its trend line (e.g. cover that is 0 every year)
#'
#' NA if the observed anomalies don't vary (nothing to predict); 0 if only
#' the predicted ones don't (the model predicts no year-to-year change).
#' @param obs,pred Numeric vectors.
#' @param tol sd below which a vector is treated as constant. A tolerance
#'   rather than 0 because detrending a constant can leave rounding error.
#' @return The Pearson correlation, 0, or NA.
cor_obs_pred <- function(obs, pred, tol = 1e-8) {
  if (sd(obs) < tol) return(NA_real_)
  if (sd(pred) < tol) return(0)
  cor(obs, pred)
}

#' Binned scatter of predicted vs observed, one panel per model and
#' prediction type, with a dashed 1:1 line and the correlation (r)
#' @param d Data frame with columns model, pred_type, and the two below.
#' @param obs,pred Names of the observed and predicted columns.
#' @param title Plot title.
#' @return A ggplot.
plot_obs_pred <- function(d, obs, pred, title) {
  stopifnot(c("model", "pred_type", obs, pred) %in% names(d))
  lims <- quantile(c(d[[obs]], d[[pred]]), c(0.005, 0.995))  # drop extremes
  r <- d |>
    summarise(r = cor(.data[[obs]], .data[[pred]]),
              .by = c(model, pred_type)) |>
    mutate(label = paste("r =", round(r, 2)))

  ggplot(d, aes(.data[[obs]], .data[[pred]])) +
    geom_bin2d(bins = 80) +
    geom_abline(linetype = 2) +
    geom_text(data = r, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.2, vjust = 1.5, inherit.aes = FALSE) +
    scale_fill_viridis_c(trans = "log10") +
    facet_grid(pred_type ~ model) +
    coord_cartesian(xlim = lims, ylim = lims) +
    labs(title = title)
}

#' Map of a per-cell value, with a diverging colour scale from -lim to lim
#' @param d Data frame with x, y and the column to map.
#' @param col Column to map.
#' @param lim Colour scale limit.
#' @param title Plot title.
#' @return A ggplot.
plot_cell_map <- function(d, col, lim, title) {
  sf::st_as_sf(d, coords = c("x", "y"), crs = terra::crs(read_mask())) |>
    plot_points_conus(color_var = col, point_size = 0.3, title = title,
                      colorscale = colorscale_diverging(
                        name = "% cover\nper year", limits = c(-lim, lim),
                        aesthetics = "colour"))
}

# by cover model ----------------------------------------------------------

train <- read_cover_training(vc)
stopifnot(c("cell", "year", "x", "y", "sources") %in% names(train))

for (cover_model in cover_models) {

  label <- paste(cover_model, vc, paste(vmcs, collapse = "-"), sep = "_")
  fits <- map(vmcs, \(v) read_cover_model("cover", cover_model, vc, v)) |>
    set_names(vmcs)
  # same response, so the observed cover below applies to every model
  stopifnot(n_distinct(map_chr(fits, \(f) f$config$spec$response)) == 1)

  # * observed and predicted cover ----------------------------------------

  # observed cover (%) at RAP cell-years
  rap <- filter(train, sources == "rap", year %in% years)
  rap$obs <- response_pct(rap, fits[[1]]$config$spec)

  # one row per cell, year, model and prediction type
  d <- imap(fits, \(f, v) {
    select(f$predictions, cell, year, all_of(pred_cols)) |>
      mutate(model = v)
  }) |>
    bind_rows() |>
    pivot_longer(all_of(pred_cols), names_to = "pred_type",
                 values_to = "pred") |>
    inner_join(select(rap, cell, year, x, y, obs), by = c("cell", "year")) |>
    filter(!is.na(obs), !is.na(pred))
  stopifnot(!anyDuplicated(d[c("cell", "year", "model", "pred_type")]))

  # keep cells with every year for every model and prediction type
  d <- filter(d, n() == length(years) * length(vmcs) * length(pred_cols),
              .by = cell)
  message(cover_model, ": ", n_distinct(d$cell), " of ",
          n_distinct(rap$cell), " RAP cells have all years")
  stopifnot(n_distinct(d$cell) > n_example_cells)

  # * trends and anomalies ------------------------------------------------

  d <- mutate(d,
              obs_anom = detrend(year, obs),
              pred_anom = detrend(year, pred),
              .by = c(cell, model, pred_type))

  by_cell <- summarise(d,
                       obs_trend = slope(year, obs),
                       pred_trend = slope(year, pred),
                       anom_cor = cor_obs_pred(obs_anom, pred_anom),
                       .by = c(cell, x, y, model, pred_type))

  summary_tbl <- left_join(
    summarise(by_cell,
              n_cells = n(),
              trend_cor = cor(obs_trend, pred_trend),
              anom_cor_median = median(anom_cor, na.rm = TRUE),
              n_cells_no_cor = sum(is.na(anom_cor)),  # constant observed
              .by = c(model, pred_type)),
    summarise(d,
              anom_cor_pooled = cor(obs_anom, pred_anom),
              anom_sd_ratio = sd(pred_anom) / sd(obs_anom),
              .by = c(model, pred_type)),
    by = c("model", "pred_type"))

  print(summary_tbl)
  write_csv(summary_tbl,
            file.path(out_dir, paste0("trend_summary_", label, ".csv")))

  # * figures -------------------------------------------------------------

  pdf(file.path(out_dir, paste0("trend_", label, ".pdf")), width = 10,
      height = 7)

  ex_cells <- sample(unique(d$cell), n_example_cells)
  g_ex <- filter(d, cell %in% ex_cells, pred_type == "pred_oof") |>
    ggplot(aes(x = year)) +
    geom_line(aes(y = pred, colour = model)) +
    geom_line(aes(y = obs)) +
    facet_wrap(~cell, scales = "free_y") +
    labs(y = "Cover (%)",
         title = paste(cover_model, "- example cells: observed (black) and",
                       "out-of-fold predicted cover"))
  print(g_ex)

  print(plot_obs_pred(by_cell, "obs_trend", "pred_trend",
                      paste(cover_model, "- temporal trend (% cover per year)")))

  print(plot_obs_pred(d, "obs_anom", "pred_anom",
                      paste(cover_model, "- temporal anomaly (% cover)")))

  g_cor <- ggplot(by_cell, aes(anom_cor)) +
    geom_histogram(bins = 50) +
    geom_text(data = summary_tbl,
              aes(x = -Inf, y = Inf,
                  label = paste("median =", round(anom_cor_median, 2))),
              hjust = -0.1, vjust = 1.5, inherit.aes = FALSE) +
    facet_grid(pred_type ~ model) +
    labs(x = "Per-cell correlation of observed and predicted anomalies",
         y = "Cells", title = cover_model) +
    geom_vline(xintercept = 0, color = 'black', linetype = 2) +
    geom_vline(data = summary_tbl, aes(xintercept = anom_cor_median),
               color = 'red', linetype = 2)
  print(g_cor)

  dev.off()

  # * maps ----------------------------------------------------------------


  by_cell$diff <- by_cell$pred_trend - by_cell$obs_trend
  lim <- max(abs(quantile(by_cell$obs_trend, c(0.01, 0.99))))
  lim_diff <- max(abs(quantile(by_cell$diff, c(0.01, 0.99))))
  
  obs_map <- distinct(by_cell, cell, x, y, obs_trend)
  stopifnot(!anyDuplicated(obs_map$cell))
  
  # top row: observed trend, then each model's predicted trend (one colour
  # scale); bottom row: each model's trend residual, below its predicted trend
  by_series <- split(by_cell, paste(by_cell$model, by_cell$pred_type))
  top <- c(list(plot_cell_map(obs_map, "obs_trend", lim, "Observed trend")),
           unname(imap(by_series, \(d, s) {
             plot_cell_map(d, "pred_trend", lim, paste("Predicted trend,", s))
           })))
  
  bottom <- c(list(patchwork::plot_spacer()),
              unname(imap(by_series, \(d, s) {
                plot_cell_map(d, "diff", lim_diff,
                              paste("Trend residual (predicted - observed),",
                                    s))
              })))
  
  g_maps <- patchwork::wrap_plots(c(top, bottom), nrow = 2) +
    patchwork::plot_annotation(
      title = paste0(cover_model, " - per-cell trends, ", min(years), "-",
                     max(years)))
  
  ggsave(file.path(out_dir, paste0("trend_maps_", label, ".png")), g_maps,
         width = 6 * length(top), height = 9, dpi = 200)
}
