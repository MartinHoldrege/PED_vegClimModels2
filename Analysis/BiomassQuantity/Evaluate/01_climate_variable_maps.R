## Purpose:
# For each climate variable, create a 6-cell figure comparing current climate
# to two end-of-century projections (BNU-ESM, IPSL-CM5A-MR; RCP 8.5):
#   top row:    current, future BNU-ESM, future IPSL-CM5A-MR  (absolute, magma)
#   bottom row: future - current for each model, then a density cell
# The 6th cell stacks two density plots: absolute values (3 scenarios) over
# differences (2 models). Color scales are pooled and shared within each row.
# Figures saved as  PNGs.
#
# June, 2026

# dependencies ------------------------------------------------------------

source("Functions/init.R")
library(patchwork)
source_functions()

# params ------------------------------------------------------------------
test_run <- FALSE
run_vpd_only <- FALSE
trunc_probs <- c(0.05, 0.95)   # color-limit truncation (squished by colorscale)
dens_probs  <- c(0.01, 0.99)   # density-panel x-limit truncation
rcp         <- "rcp85"
models      <- c("BNU-ESM", "IPSL-CM5A-MR")

# scenario colors, reused across both density panels
scen_colors <- c("Current"      = "grey30",
                 "BNU-ESM"      = "#2c7bb6",   # cool/wet-ish, blue
                 "IPSL-CM5A-MR" = "#d7191c")   # hot/dry, red

out_dir <- file.path("Figures/Climate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_width  <- 14
fig_height <- 8
dpi        <- 600

# load climate rasters ----------------------------------------------------
# read_climate_raster() applies climate_name_lookup() (short names). Drop soil
# layers (path_soil = NULL) so we map only climate variables.

r_cur <- read_climate_raster("current",  path_soil = NULL)
r_m1  <- read_climate_raster(models[1],   path_soil = NULL)
r_m2  <- read_climate_raster(models[2],   path_soil = NULL)

vars <- Reduce(intersect, list(names(r_cur), names(r_m1), names(r_m2)))

if(test_run) vars <- vars[1]

if(run_vpd_only) vars <- str_subset(vars, 'VPD')

# helpers -----------------------------------------------------------------

# Pooled limits across layers: min of lower-prob quantiles, max of upper.
# (Per-layer quantiles combined via min/max, not a true pooled quantile;
# adequate for setting visual limits.)
pooled_limits <- function(rast_list, probs) {
  lo <- min(vapply(rast_list, raster_quantile, numeric(1), prob = probs[1]))
  hi <- max(vapply(rast_list, raster_quantile, numeric(1), prob = probs[2]))
  c(lo, hi)
}

# loop over variables -----------------------------------------------------

for (v in vars) {
  rc <- r_cur[[v]]
  r1 <- r_m1[[v]]
  r2 <- r_m2[[v]]

  d1 <- r1 - rc   # future - current, model 1
  d2 <- r2 - rc   # future - current, model 2

  # --- absolute maps: shared magma scale ---
  abs_col_lim <- pooled_limits(list(rc, r1, r2), trunc_probs)
  abs_scale <- colorscale_biomass(name = v, limits = abs_col_lim,
                                  option = "magma")

  g_cur <- plot_map_conus(rc, colorscale = abs_scale, title = "Current")
  g_m1  <- plot_map_conus(r1, colorscale = abs_scale,
                          title = paste("Future:", models[1]))
  g_m2  <- plot_map_conus(r2, colorscale = abs_scale,
                          title = paste("Future:", models[2]))

  # --- difference maps: shared symmetric diverging scale ---
  diff_q <- max(abs(pooled_limits(list(d1, d2), trunc_probs)))
  diff_scale <- colorscale_diverging(name = paste0("\u0394 ", v),
                                     limits = c(-diff_q, diff_q))

  g_d1 <- plot_map_conus(d1, colorscale = diff_scale,
                         title = paste(models[1], "\u2212 Current"))
  g_d2 <- plot_map_conus(d2, colorscale = diff_scale,
                         title = paste(models[2], "\u2212 Current"))

  # --- density panels for the 6th cell ---
  abs_dens_lim  <- pooled_limits(list(rc, r1, r2), dens_probs)
  diff_dens_lim <- pooled_limits(list(d1, d2), dens_probs)

  g_dens_abs <- climate_density_panel(
    rast_list = setNames(list(rc, r1, r2), c("Current", models)),
    colors = scen_colors,
    xlim = abs_dens_lim,
    xlab = v
  )
  g_dens_diff <- climate_density_panel(
    rast_list = setNames(list(d1, d2), models),
    colors = scen_colors,
    xlim = diff_dens_lim,
    vline = 0,
    xlab = paste0("\u0394 ", v)
  )
  g_dens <- g_dens_abs / g_dens_diff   # stacked, absolute on top

  # --- assemble: top row 3 maps, bottom row 2 maps + density cell ---
  design <- "ABC\nDEF"
  panel <- g_cur + g_m1 + g_m2 + g_d1 + g_d2 + g_dens +
    plot_layout(design = design, guides = "collect") +
    plot_annotation(title = v)

  out_file <- file.path(out_dir, paste0("climate_comparison_", v, ".png"))
  ggsave(out_file, panel, width = fig_width, height = fig_height,
         dpi = dpi, bg = "white")
  message("Wrote ", out_file)
}
