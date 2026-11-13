# soilwat_cover_compare.R
#
# Compares the herbaceous and shrub cover currently used as SOILWAT2 input
# with our cover models and with RAP. The SOILWAT2 inputs are fractions of
# herbaceous + shrub cover (they sum to 1), so every source is put on that
# scale:
#   herbaceous fraction = herbaceous / (herbaceous + shrub)
#   shrub fraction      = shrub / (herbaceous + shrub)
# Sources:
#   SOILWAT2: herbfrac and shrubfrac (gridMET 4 km grid, western US,
#     1979-2020 climatology)
#   model: current-climate predictions (1991-2020 normals, i.e. target year
#     2021) of the herbaceous and shrub cover models chosen below
#     (03_predict_raster.R); wall-to-wall (not masked)
#   RAP: 2021 cover, LCMAP- and 2021 fire-masked (read_rap_cover_raster())
# Compared on the SOILWAT2 4 km grid: the 1 km cover (model and RAP) is
# averaged onto the 4 km cells, then turned into fractions. Each comparison
# uses the cells where both of its sources have a value, so model vs SOILWAT2
# is wall-to-wall, and any comparison with RAP is limited to RAP's masked
# area.
#
# Outputs (out_dir), <label> = the model versions:
#   soilwat_compare_maps_<label>.png: herbaceous and shrub fraction of each
#     source, and the differences in herbaceous fraction (shrub fraction
#     differences are their negative), with summaries in the subtitles
#   soilwat_compare_scatter_<label>.png: herbaceous fraction of each source
#     against the others, with the 1:1 and least-squares lines, and summaries
#     in the subtitles
#
# November 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
theme_set(theme_custom1())

# params ------------------------------------------------------------------

# cover models (data version and model version), as in cover_pred_path()
herb_vc <- "c08"
herb_vmc <- "m05.51"
shrub_vc <- "c08"
shrub_vmc <- "m05.51"

rap_year <- 2021

out_dir <- "Figures/Cover/soilwat_compare"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
label <- paste0("herb-", herb_vc, "-", herb_vmc, "_shrub-", shrub_vc, "-",
                shrub_vmc)

# functions ---------------------------------------------------------------

#' Herbaceous and shrub fractions of their summed cover
#' @param herb,shrub SpatRasters of herbaceous and shrub cover (same units).
#' @return SpatRaster with layers herb and shrub (fractions that sum to 1;
#'   NA where both covers are 0).
#' @examples
#' herb <- terra::rast(nrows = 1, ncols = 3, vals = c(10, 0, 0))
#' shrub <- terra::rast(nrows = 1, ncols = 3, vals = c(30, 5, 0))
#' to_fractions(herb = herb, shrub = shrub)
to_fractions <- function(herb, shrub) {
  total <- herb + shrub
  out <- c(herb / total, shrub / total)
  names(out) <- c("herb", "shrub")
  out
}

#' Map of a fraction (0 to 1)
#' @param r Single-layer SpatRaster.
#' @param title Plot title.
#' @return A ggplot.
plot_frac_map <- function(r, title) {
  plot_map_conus(r, title = title,
                 colorscale = colorscale_cover(name = "Fraction",
                                               limits = c(0, 1)))
}

#' Map of a difference in fraction, with its mean and mean absolute value in
#' the subtitle
#' @param r Single-layer SpatRaster.
#' @param lim Colour scale runs from -lim to lim.
#' @param title Plot title.
#' @return A ggplot.
plot_diff_map <- function(r, lim, title) {
  v <- terra::values(r, mat = FALSE, na.rm = TRUE)
  plot_map_conus(r, title = title,
                 colorscale = colorscale_diverging(name = "Difference",
                                                   limits = c(-lim, lim))) +
    labs(subtitle = paste0("mean = ", round(mean(v), 3),
                           ", mean |difference| = ", round(mean(abs(v)), 3)))
}

#' Binned scatter of one source's herbaceous fraction against another's
#'
#' Uses the cells where both sources have a value. With the 1:1 line (dashed)
#' and the least-squares line (red). Subtitle: correlation, slope of the
#' least-squares line, mean difference (y - x), mean absolute difference, and
#' number of cells.
#' @param d Data frame with one column per source.
#' @param x,y Names of the sources on the x and y axes.
#' @return A ggplot.
plot_pair <- function(d, x, y) {
  d <- d[!is.na(d[[x]]) & !is.na(d[[y]]), ]
  difference <- d[[y]] - d[[x]]
  fit <- lm(d[[y]] ~ d[[x]])
  subtitle <- paste0("r = ", round(cor(d[[x]], d[[y]]), 2),
                     ", slope = ", round(coef(fit)[[2]], 2),
                     ", mean difference (y - x) = ", round(mean(difference), 3),
                     "\nmean |difference| = ", round(mean(abs(difference)), 3),
                     ", n = ", nrow(d), " cells")

  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_bin2d(bins = 60) +
    geom_abline(linetype = 2) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "red", se = FALSE) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = paste(x, "herbaceous fraction"),
         y = paste(y, "herbaceous fraction"),
         subtitle = subtitle)
}

# read --------------------------------------------------------------------

# SOILWAT2 inputs: herbfrac and shrubfrac, which sum to 1
sw_files <- list.files(file.path(paths$large, "data_raw/soilwat_input_cover"),
                       pattern = "\\.nc$", full.names = TRUE)
sw <- terra::rast(sw_files)[[c("herbfrac", "shrubfrac")]]
names(sw) <- c("herb", "shrub")
sw_sum <- terra::values(sw$herb + sw$shrub, mat = FALSE, na.rm = TRUE)
stopifnot(abs(sw_sum - 1) < 0.001)

# our models, current climate (% cover), wall-to-wall
herb_model <- cover_pred_path("cover", "herbaceous", herb_vc, herb_vmc,
                              "current") |>
  terra::rast()
shrub_model <- cover_pred_path("cover", "shrub", shrub_vc, shrub_vmc,
                               "current") |>
  terra::rast()
herb_model <- herb_model[["cover"]]
shrub_model <- shrub_model[["cover"]]

# RAP (% cover), LCMAP- and fire-masked
herb_rap <- read_rap_cover_raster("herbaceous", years = rap_year)
shrub_rap <- read_rap_cover_raster("shrub", years = rap_year)

# fractions on the SOILWAT2 grid -------------------------------------------

# 1 km cover averaged onto the 4 km SOILWAT2 cells (ignoring NA 1 km cells),
# then fractions
frac <- list(
  SOILWAT2 = sw,
  model = to_fractions(terra::project(herb_model, sw, method = "average"),
                       terra::project(shrub_model, sw, method = "average")),
  RAP = to_fractions(terra::project(herb_rap, sw, method = "average"),
                     terra::project(shrub_rap, sw, method = "average"))
)

# differences in herbaceous fraction (NA where either source is NA, so those
# with RAP are limited to RAP's masked area)
diffs <- list(
  "model - SOILWAT2" = frac$model$herb - frac$SOILWAT2$herb,
  "model - RAP"      = frac$model$herb - frac$RAP$herb,
  "SOILWAT2 - RAP"   = frac$SOILWAT2$herb - frac$RAP$herb
)

# maps ----------------------------------------------------------------------

# one colour scale for all three difference maps
diff_values <- unlist(map(diffs, \(r) terra::values(r, mat = FALSE,
                                                    na.rm = TRUE)))
lim <- max(abs(quantile(diff_values, c(0.01, 0.99))))

herb_maps <- imap(frac, \(r, source) {
  plot_frac_map(r$herb, paste(source, "- herbaceous fraction"))
})
shrub_maps <- imap(frac, \(r, source) {
  plot_frac_map(r$shrub, paste(source, "- shrub fraction"))
})
diff_maps <- imap(diffs, \(r, name) {
  plot_diff_map(r, lim, paste("Herbaceous fraction,", name))
})

g_maps <- patchwork::wrap_plots(unname(c(herb_maps, shrub_maps, diff_maps)),
                                ncol = 3) +
  patchwork::plot_annotation(
    title = paste0("Fractions of herbaceous + shrub cover. Model: herbaceous ",
                   herb_vc, "-", herb_vmc, ", shrub ", shrub_vc, "-",
                   shrub_vmc, "; RAP ", rap_year),
    caption = "Shrub fraction differences are the negative of the herbaceous ones")

ggsave(file.path(out_dir, paste0("soilwat_compare_maps_", label, ".png")),
       g_maps, width = 18, height = 16, dpi = 200)

# scatter plots -------------------------------------------------------------

# one row per 4 km cell; NA where a source has no value (plot_pair() uses the
# cells where both of its sources have one)
herb_frac <- tibble(SOILWAT2 = terra::values(frac$SOILWAT2$herb, mat = FALSE),
                    model = terra::values(frac$model$herb, mat = FALSE),
                    RAP = terra::values(frac$RAP$herb, mat = FALSE))

g_scatter <- patchwork::wrap_plots(
  plot_pair(herb_frac, x = "SOILWAT2", y = "model"),
  plot_pair(herb_frac, x = "RAP", y = "model"),
  plot_pair(herb_frac, x = "RAP", y = "SOILWAT2"),
  nrow = 1
) +
  patchwork::plot_annotation(
    title = paste0("Herbaceous fraction of herbaceous + shrub cover, 4 km ",
                   "cells. Model: herbaceous ", herb_vc, "-", herb_vmc,
                   ", shrub ", shrub_vc, "-", shrub_vmc, "; RAP ", rap_year),
    caption = paste("Dashed: 1:1 line; red: least-squares line. The shrub",
                    "fraction gives the same r and slope, and the negative",
                    "of the differences"))

ggsave(file.path(out_dir, paste0("soilwat_compare_scatter_", label, ".png")),
       g_scatter, width = 18, height = 7, dpi = 200)
