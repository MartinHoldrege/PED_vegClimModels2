# rap_vs_field_tree_cover.R
#
# Compares RAP tree cover with the field tree cover in the c01 training data,
# at the same pixels. Motivation: needleleaf tree cover from the combined
# model (RAP-based tree cover x modelled needleleaf share) is higher than
# observed needleleaf cover in c01, which suggests RAP tree cover is higher
# than field tree cover where the needleleaf data are.
#
# Two sets of pixels:
#   needle  pixel-years with a field needleleaf fraction (frac_needle not NA)
#   all     every pixel-year with field tree cover
#
# RAP is the 2021 tree cover raster (the one c02 was sampled from), so field
# years (2000-2024) don't match it. Pixels where the RAP zero-tree layer says
# zero tree (fracZeroTree > .thresh_zt) are handled two ways, each with its
# own set of figures:
#   dropped  those pixels left out: only where RAP isn't forced to 0
#   zeroed   kept, with RAP set to 0, as in the combined tree cover
#
# Field tree cover is averaged over years at each pixel (RAP is one value per
# pixel). Difference = RAP - field, in percentage points.
#
# Outputs (PNG, fig_dir; <set> is needle or all, <zt> dropped or zeroed):
#   rap_minus_field_tree_<set>_<zt>.png  histograms, all sources and by
#                                             source, with the mean difference
#                                             (bias) and map of difference
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")  # .thresh_zt
library(patchwork)

# params ------------------------------------------------------------------

fig_dir <- "Figures/Cover/rap_vs_field"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# data --------------------------------------------------------------------

field_vc <- 'c01'
field <- read_cover_training(field_vc) |>
  filter(!is.na(cov_tree))

snap <- read_mask()
rap_tree <- terra::rast(file.path(paths$large, "Data_processed/CoverData/rap",
                                  "RAP_v3_cover_2021_1000m.tif"))[["tree"]] |>
  align_raster(snap)  # GEE exports are one cell larger than the snap grid
frac_zero <- read_zero_tree_raster("fracZeroTree") |>
  align_raster(snap)

# RAP at each field pixel (x, y are snap-grid cell centres)
xy <- as.matrix(field[c("x", "y")])
field$rap_tree <- terra::extract(rap_tree, xy)[[1]]
field$zero_tree <- 100 * terra::extract(frac_zero, xy)[[1]] > .thresh_zt

field <- field |>
  filter(!is.na(rap_tree), !is.na(zero_tree), sources != 'rap') |>
  mutate(rap_tree = if_else(zero_tree, 0, rap_tree))

message("pixel-years in RAP zero-tree pixels: ", sum(field$zero_tree), " of ",
        nrow(field))

#' One row per pixel: mean field tree cover over years, RAP, and RAP - field
#'
#' @param d Pixel-years with cov_tree, rap_tree, sources, cell, x, y.
#' @return Tibble; `source` is the pixel's field source, or "mixed" if more
#'   than one contributed.
by_pixel <- function(d) {
  d |>
    summarise(field_tree = mean(cov_tree),
              rap_tree = first(rap_tree),
              sources = paste(sort(unique(unlist(str_split(sources, ",")))),
                              collapse = ","),
              .by = c(cell, x, y)) |>
    mutate(diff = rap_tree - field_tree,
           source = if_else(str_detect(sources, ","), "mixed", sources))
}

# every combination of pixel set and zero-tree handling
combos <- expand_grid(set = c("needle", "all"), zt = c("dropped", "zeroed"))
px <- pmap(combos, \(set, zt) {
  d <- if (zt == "dropped") filter(field, !zero_tree) else field
  if (set == "needle") d <- filter(d, !is.na(frac_needle))
  by_pixel(d)
})
mutate(combos, n_pixels = map_int(px, nrow))

# figures -----------------------------------------------------------------

set_title <- c(needle = "Pixels with field needleleaf fraction",
               all    = "All pixels with field tree cover")
zt_title <- c(dropped = "RAP zero-tree pixels dropped",
              zeroed  = "RAP set to 0 in zero-tree pixels")
lim <- map_dbl(px, \(d) max(abs(quantile(d$diff, c(0.1, 0.9), names = FALSE)))) |> 
  max()
plots <- pmap(mutate(combos, d = px), \(set, zt, d) {
  title <- paste0(set_title[[set]], " (", zt_title[[zt]], ")")

  # histograms: all sources, then each source
  d_long <- bind_rows(mutate(d, source = "all sources"), d) |>
    mutate(source = fct_relevel(factor(source), "all sources"))
  bias <- d_long |>
    summarise(bias_mean = mean(diff), 
              bias_median = median(diff), 
              n = n(), .by = source) |>
    mutate(label = paste0("mean = ", round(bias_mean, 1),
                          "\nmedian = ", round(bias_median, 1),
                          "\nn = ", n))

  g_hist <- ggplot(d_long, aes(diff)) +
    geom_histogram(bins = 50, fill = "grey60") +
    geom_vline(xintercept = 0, colour = "grey30") +
    geom_vline(data = bias, aes(xintercept = bias_mean), colour = "darkred",
               linetype = 2) +
    geom_text(data = bias, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.3, size = 3) +
    facet_wrap(~ source, scales = "free_y") +
    labs(x = "RAP - field tree cover (% points)", y = "Pixels",
         title = title,
         subtitle = "RAP 2021 vs field (mean over years); dashed: bias (mean difference)")

  # map
  
  g_map <- d |>
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(snap)) |>
    plot_points_conus(color_var = "diff", point_size = 0.15,
                      colorscale = colorscale_diverging(
                        name = "RAP - field (all sources)\n(% points)",
                        limits = c(-lim, lim), aesthetics = "colour"),
                      title = NULL,
                      inset_hist = TRUE, inset_drop_zero = FALSE)
  
  # wrap_elements() so the border survives nesting in wrap_plots()
  wrap_plots(g_hist, g_map, ncol = 1, heights = c(0.7, 1))
})

comb <- patchwork::wrap_plots(plots, nrow = 2)
ggsave(file.path(fig_dir, paste0("rap_minus_field_tree_", field_vc, ".png")),
       comb, width = 25, height = 20, dpi = 300)
