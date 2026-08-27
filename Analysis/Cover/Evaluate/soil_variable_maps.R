## Purpose:
# One multi-panel figure showing each soil covariate across CONUS.
# Panels have independent color scales (variables differ in units), with
# limits truncated at the 1st/99th percentiles and out-of-range values
# squished
#
# August, 2026

# dependencies ------------------------------------------------------------

source("Functions/init.R")
library(patchwork)
source_functions()

# params ------------------------------------------------------------------

# test_run = TRUE produces the same figure from a ~10 x 10 cell version of the
# raster, at low dpi, written to a separate file. For checking layout,
# labels and panel order without the full-resolution render.
test_run <- FALSE

path_soil <- file.path(paths$large, "Data_processed/soils",
                       "soil_covariates_solus100_1000m.tif")

out_dir <- file.path("Figures/Soils")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trunc_probs <- c(0.01, 0.99)   # color-limit truncation (squished by colorscale)
n_col       <- 3
fig_width   <- 16
fig_height  <- 9
dpi         <- if (test_run) 100 else 900
# 4381 x 2733 cells total; 5e6 is near-native for a 5-inch panel at 600 dpi
maxcell     <- if (test_run) 1e3 else 5e6

# panel order, titles, and legend units
var_info <- tibble::tribble(
  ~var,           ~label,                             ~units,
  "sand",         "Sand (profile mean)",              "%",
  "clay",         "Clay (profile mean)",              "%",
  "clay_surface", "Clay (0-3 cm)",                    "%",
  "coarse",       "Coarse fragments (profile mean)",  "% vol.",
  "carbon",       "Organic carbon (0-3 cm)",          "%",
  "soilDepth",    "Depth to restriction",             "cm",
  "AWHC",         "Available water holding capacity", "cm"
)

# load --------------------------------------------------------------------

r <- terra::rast(path_soil)

stopifnot(all(var_info$var %in% names(r)))

# coarse subsample so quantiles and rendering are fast; keeps extent and CRS
if (test_run) {
  r <- terra::spatSample(r, size = 100, method = "regular", as.raster = TRUE)
  message("test_run: using ", ncol(r), " x ", nrow(r), " cells")
}

# panels ------------------------------------------------------------------

panels <- purrr::pmap(var_info, function(var, label, units) {
  lyr <- r[[var]]

  limits <- c(raster_quantile(lyr, trunc_probs[1]),
              raster_quantile(lyr, trunc_probs[2]))

  plot_map_conus(
    lyr,
    colorscale = colorscale_biomass(name = units, limits = limits,
                                    option = "viridis"),
    title = label,
    maxcell = maxcell
  )
})

# assemble and write ------------------------------------------------------

fig <- patchwork::wrap_plots(panels, ncol = n_col)

out_file <- file.path(out_dir, "soil_covariate_maps.png")
ggsave(out_file, fig, width = fig_width, height = fig_height,
       dpi = dpi, bg = "white")
