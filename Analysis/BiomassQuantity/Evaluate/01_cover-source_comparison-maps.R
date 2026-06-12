# Purpose:
# Create comparison maps of vegetation cover from RAP vs modelled (GLM)
# sources. Outputs one PNG per functional group showing RAP cover,
# modelled cover, and the difference. Difference maps are masked to
# the same valid areas used for model training.
#
# script started: June, 2026

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
library(patchwork)

# params ------------------------------------------------------------------

out_dir <- file.path("Figures/CoverDatFigures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_width <- 15
fig_height <- 10

# load cover rasters -------------------------------------------------------

r_rap <- load_cover("rap")
r_mod <- load_cover("model")

# load masks (same as used in 06_sample_training_data.R)
rasters <- load_conus_rasters()
mask_lcmap <- rasters$mask_lcmap
mask_fire <- rasters$mask_fire

# align extents
aligned <- align_raster_extents(rast_list = list(
  rap = r_rap, model = r_mod, lcmap = mask_lcmap, fire = mask_fire
))
r_rap <- aligned$rap
r_mod <- aligned$model
mask_lcmap <- aligned$lcmap
mask_fire <- aligned$fire

cover_cols <- names(r_rap)

# which covers are "woody" (get fire mask too)
woody_cols <- c("totalTreeCov", "totalShrubCov")


# create maps per functional group -----------------------------------------

for (col in cover_cols) {
  cat("Mapping:", col, "\n")
  
  r_rap_i <- r_rap[[col]]
  r_mod_i <- r_mod[[col]]
  r_diff <- r_mod_i - r_rap_i
  names(r_diff) <- "difference"
  
  # mask difference map: LCMAP for all, + fire mask for woody
  r_diff <- terra::mask(r_diff, mask_lcmap)
  if (col %in% woody_cols) {
    r_diff <- terra::mask(r_diff, mask_fire)
  }
  
  # shared cover limits (99th percentile across both)
  cov_max <- max(
    terra::global(r_rap_i, fun = function(x) quantile(x, 0.99, na.rm = TRUE))[[1]],
    terra::global(r_mod_i, fun = function(x) quantile(x, 0.99, na.rm = TRUE))[[1]]
  )
  cov_limits <- c(0, cov_max)
  
  # symmetric difference limits (from masked diff)
  diff_q99 <- terra::global(
    abs(r_diff),
    fun = function(x) quantile(x, 0.99, na.rm = TRUE)
  )[[1]]
  
  # clean label
  pft_label <- stringr::str_remove(col, "Cov$") |>
    stringr::str_replace("total", "Total ")
  
  mask_note <- if (col %in% woody_cols) "LCMAP + fire mask" else "LCMAP mask"
  
  # maps
  g_rap <- plot_map_conus(
    r_rap_i,
    colorscale = colorscale_cover(name = "Cover", limits = cov_limits),
    title = paste(pft_label, "— RAP")
  )
  
  g_mod <- plot_map_conus(
    r_mod_i,
    colorscale = colorscale_cover(name = "Cover", limits = cov_limits),
    title = paste(pft_label, "— Modelled")
  )
  
  g_diff <- plot_map_conus(
    r_diff,
    colorscale = colorscale_diverging(
      name = "Mod \u2212 RAP",
      limits = c(-diff_q99, diff_q99)
    ),
    title = paste(pft_label, "— Difference (", mask_note, ")")
  )
  
  # layout: RAP and modelled side by side, difference below
  g_combined <- (g_rap + g_mod) / g_diff +
    plot_layout(heights = c(1, 1))+
    plot_annotation(title = paste("Cover comparison:", pft_label)) &
    theme(plot.margin = margin(0, 0, 0, 0)) 
  
  # save
  p_out <- file.path(out_dir, paste0("cover_comparison_", col, ".png"))
  ggsave(p_out, g_combined, width = fig_width, height = fig_height,
         dpi = 600, bg = "white")
  
}

