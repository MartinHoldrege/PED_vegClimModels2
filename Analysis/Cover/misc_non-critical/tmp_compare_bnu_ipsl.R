# Temporary: compare BNU-ESM and IPSL-CM5A-MR (MACAv2-METDATA, RCP 8.5,
# 2070-2099) for MAT and MAP. For each variable, a 3-row figure:
#   row 1: future - historical, for each GCM (shared colour scale)
#   row 2: IPSL - BNU for the historical period, and for the future period
#          (shared colour scale)
#   row 3: difference of the changes, IPSL change - BNU change (row 1 right
#          minus row 1 left; also equals row 2 right minus row 2 left)
# MACA is bias-corrected to gridMET over the historical period, so the
# historical IPSL - BNU difference should be near zero.

# dependencies ------------------------------------------------------------

source("Functions/init.R")
library(patchwork)
source_functions()

# params ------------------------------------------------------------------

gcms <- c("BNU-ESM", "IPSL-CM5A-MR")
vars <- c("MAT", "MAP")
fut <- "rcp85_2070-2099"
hist <- "historical_1971-2000"
trunc_probs <- c(0.02, 0.98) # colour limits (values beyond are squished)

maca_dir <- file.path(paths$large, "Data_processed/WallToWallClimateData/maca")
out_dir <- "Figures/Climate/tmp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# read --------------------------------------------------------------------

#' Read one MACA normals raster, keeping `vars`
#'
#' Longitudes are rotated from 0-360 to -180-180 so the state outlines in
#' plot_map_conus() line up.
#'
#' @param gcm GCM name.
#' @param scen_period e.g. "rcp85_2070-2099".
#' @return SpatRaster with layers `vars`.
read_maca <- function(gcm, scen_period) {
  p <- file.path(maca_dir, paste0("maca_", gcm, "_", scen_period, "_CLIM.tif"))
  stopifnot(file.exists(p))
  r <- terra::rast(p)
  names(r) <- climate_name_lookup(names(r))
  terra::rotate(r[[vars]])
}

r_hist <- purrr::map(setNames(gcms, gcms), read_maca, scen_period = hist)
r_fut <- purrr::map(setNames(gcms, gcms), read_maca, scen_period = fut)

# helpers -----------------------------------------------------------------

#' Symmetric colour limits pooled across rasters
#'
#' @param rast_list List of single-layer SpatRasters.
#' @return c(-q, q), where q is the largest absolute truncated quantile.
sym_limits <- function(rast_list) {
  q <- purrr::map_dbl(rast_list, \(r) max(abs(c(
    raster_quantile(r, trunc_probs[1]),
    raster_quantile(r, trunc_probs[2])))))
  c(-max(q), max(q))
}

# figures -----------------------------------------------------------------

for (v in vars) {
  d_bnu <- r_fut[[gcms[1]]][[v]] - r_hist[[gcms[1]]][[v]]
  d_ipsl <- r_fut[[gcms[2]]][[v]] - r_hist[[gcms[2]]][[v]]
  hist_diff <- r_hist[[gcms[2]]][[v]] - r_hist[[gcms[1]]][[v]]
  fut_diff <- r_fut[[gcms[2]]][[v]] - r_fut[[gcms[1]]][[v]]
  # difference of the changes; equals fut_diff - hist_diff
  dd <- d_ipsl - d_bnu
  
  scale_delta <- colorscale_diverging(name = paste("Change in", v),
                                      limits = sym_limits(list(d_bnu, d_ipsl)))
  scale_gcm <- colorscale_diverging(name = paste(v, "\nIPSL - BNU"),
                                    limits = sym_limits(list(hist_diff,
                                                             fut_diff)))
  scale_dd <- colorscale_diverging(name = paste("Change in", v,
                                                "\nIPSL - BNU"),
                                   limits = sym_limits(list(dd)))
  
  g <- plot_map_conus(d_bnu, colorscale = scale_delta,
                      title = paste(gcms[1], ": future - historical")) +
    plot_map_conus(d_ipsl, colorscale = scale_delta,
                   title = paste(gcms[2], ": future - historical")) +
    plot_map_conus(hist_diff, colorscale = scale_gcm,
                   title = "Historical (1971-2000): IPSL - BNU") +
    plot_map_conus(fut_diff, colorscale = scale_gcm,
                   title = "Future (RCP 8.5, 2070-2099): IPSL - BNU") +
    plot_map_conus(dd, colorscale = scale_dd,
                   title = "Change: IPSL - BNU (top right - top left)") +
    plot_layout(design = "AB\nCD\nE#", guides = "collect") +
    plot_annotation(title = v)
  
  out_file <- file.path(out_dir, paste0("bnu_vs_ipsl_", v, ".png"))
  ggsave(out_file, g, width = 12, height = 12, dpi = 300, bg = "white")
  message("Wrote ", out_file)
}
