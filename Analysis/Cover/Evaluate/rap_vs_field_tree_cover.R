# rap_vs_field_cover.R
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
# Shrub and herbaceous cover are compared too, zero-tree pixels kept, at two
# sets of pixels:
#   all        every pixel with field shrub (or herbaceous) cover
#   nonforest  only where RAP tree cover is < 10%: RAP can't see the
#              understorey, while field shrub and herbaceous cover include it
#
# Field cover is averaged over years at each pixel (RAP is one value per
# pixel). Difference = RAP - field, in percentage points.
#
# Outputs (PNG, fig_dir): histograms, all sources and by source, with the mean
# difference (bias), and a map of the difference
#   rap_minus_field_tree_<field_vc>.png        tree, each set x zero-tree
#                                              handling
#   rap_minus_field_shrub_<field_vc>.png       shrub, all and nonforest
#   rap_minus_field_herbaceous_<field_vc>.png  herbaceous, all and nonforest
#   rap_minus_field_tree_relationships_<field_vc>.png
#                                              RAP - field tree cover against
#                                              field needleleaf share, RAP
#                                              tree cover and field tree cover
#                                              (all sources pooled)
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
field <- read_cover_training(field_vc)

snap <- read_mask()
rap <- terra::rast(file.path(paths$large, "Data_processed/CoverData/rap",
                             "RAP_v3_cover_2021_1000m.tif"))
rap <- rap[[c("tree", "shrub", "herbaceous")]] |>
  align_raster(snap)  # GEE exports are one cell larger than the snap grid
frac_zero <- read_zero_tree_raster("fracZeroTree") |>
  align_raster(snap)

# RAP at each field pixel (x, y are snap-grid cell centres)
xy <- as.matrix(field[c("x", "y")])
field <- bind_cols(field, rename_with(terra::extract(rap, xy),
                                      \(n) paste0("rap_", n)))
field$zero_tree <- 100 * terra::extract(frac_zero, xy)[[1]] > .thresh_zt

field <- field |>
  filter(!is.na(rap_tree), !is.na(zero_tree), sources != 'rap') |>
  mutate(rap_tree = if_else(zero_tree, 0, rap_tree))

message("pixel-years in RAP zero-tree pixels: ", sum(field$zero_tree), " of ",
        nrow(field))

#' One row per pixel: mean field cover over years, RAP, and RAP - field
#'
#' @param d Pixel-years with cov_<var>, rap_<var>, sources, cell, x, y.
#' @param var "tree", "shrub" or "herbaceous".
#' @return Tibble; `source` is the pixel's field source, or "mixed" if more
#'   than one contributed.
by_pixel <- function(d, var) {
  d |>
    summarise(field = mean(.data[[paste0("cov_", var)]]),
              frac_needle = mean(frac_needle),  # NA unless every year has it
              rap = first(.data[[paste0("rap_", var)]]),
              sources = paste(sort(unique(unlist(str_split(sources, ",")))),
                              collapse = ","),
              .by = c(cell, x, y)) |>
    mutate(diff = rap - field,
           source = if_else(str_detect(sources, ","), "mixed", sources))
}

# tree: every combination of pixel set and zero-tree handling; shrub and
# herbaceous: all pixels with field cover, and those with RAP tree cover < 10%
# (zero-tree pixels kept)
combos <- bind_rows(
  expand_grid(var = "tree", set = c("needle", "all"),
              zt = c("dropped", "zeroed")),
  expand_grid(var = c("shrub", "herbaceous"), set = c("all", "nonforest"),
              zt = "none")
)
px <- pmap(combos, \(var, set, zt) {
  d <- filter(field, !is.na(.data[[paste0("cov_", var)]]))
  if (zt == "dropped") d <- filter(d, !zero_tree)
  if (set == "needle") d <- filter(d, !is.na(frac_needle))
  if (set == "nonforest") d <- filter(d, rap_tree < 10)
  by_pixel(d, var)
})
mutate(combos, n_pixels = map_int(px, nrow))

# figures -----------------------------------------------------------------

set_title <- c(needle    = "Pixels with field needleleaf fraction",
               all       = "All pixels with field",
               nonforest = "Pixels with RAP tree cover < 10% and field")
zt_title <- c(dropped = " (but RAP zero-tree pixels dropped)",
              zeroed  = " (RAP set to 0 in zero-tree pixels)",
              none    = "")
# one colour scale per cover type
lim <- map_dbl(px, \(d) max(abs(quantile(d$diff, c(0.1, 0.9), names = FALSE)))) |> 
  tapply(combos$var, max)
plots <- pmap(mutate(combos, d = px), \(var, set, zt, d) {
  title <- paste0(set_title[[set]],
                  if (set != "needle") paste0(" ", var, " cover"),
                  zt_title[[zt]])

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
    labs(x = paste0("RAP - field ", var, " cover (% points)"), y = "Pixels",
         title = title,
         subtitle = "RAP 2021 vs field (mean over years); dashed: bias (mean difference)")

  # map
  
  g_map <- d |>
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(snap)) |>
    plot_points_conus(color_var = "diff", point_size = 0.15,
                      colorscale = colorscale_diverging(
                        name = "RAP - field (all sources)\n(% points)",
                        limits = c(-lim[[var]], lim[[var]]),
                        aesthetics = "colour"),
                      title = NULL,
                      inset_hist = TRUE, inset_drop_zero = FALSE)
  
  # wrap_elements() so the border survives nesting in wrap_plots()
  wrap_plots(g_hist, g_map, ncol = 1, heights = c(0.7, 1))
})

is_tree <- combos$var == "tree"
comb <- patchwork::wrap_plots(plots[is_tree], nrow = 2)
ggsave(file.path(fig_dir, paste0("rap_minus_field_tree_", field_vc, ".png")),
       comb, width = 25, height = 20, dpi = 300)

# shrub and herbaceous: one file each, all pixels and nonforest side by side
walk(c("shrub", "herbaceous"), \(v) {
  comb <- patchwork::wrap_plots(plots[combos$var == v], nrow = 1)
  ggsave(file.path(fig_dir, paste0("rap_minus_field_", v, "_", field_vc,
                                   ".png")),
         comb, width = 25, height = 10, dpi = 300)
})

# relationships -------------------------------------------------------------

# does the RAP - field difference in tree cover depend on the needleleaf share,
# or on how much tree cover there is? All sources pooled. Note y contains x
# in the last two: if both datasets have independent error, RAP - field rises
# with RAP and falls with field cover even with no real bias pattern

pick <- \(set, zt_scatter) px[[which(combos$var == "tree" & combos$set == set &
                           combos$zt == zt_scatter)]]

scatter_diff <- function(d, x, xlab) {
  ggplot(d, aes(.data[[x]], diff)) +
    geom_point(alpha = 0.1, size = 0.2) +
    geom_hline(yintercept = 0, colour = "grey30") +
    geom_smooth(colour = "darkred", se = FALSE) +
    labs(x = xlab, y = "RAP - field tree cover",
         subtitle = paste0("n = ", nrow(d), " pixels"))
}

#' Field against RAP tree cover, with a 1:1 line
scatter_1to1 <- function(d) {
  ggplot(d, aes(rap, field)) +
    geom_point(alpha = 0.1, size = 0.2) +
    geom_abline(colour = "grey30") +
    geom_smooth(colour = "darkred", se = FALSE) +
    coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(x = "RAP tree cover (%)", y = "Field tree cover (%)",
         subtitle = paste0("n = ", nrow(d), " pixels; grey: 1:1"))
}


for (zt_scatter in c("dropped", "zeroed")) {
  d_needle <- mutate(pick("needle", zt_scatter), needle_pct = 100 * frac_needle)
  d_all <- pick("all", zt_scatter)
  
  comb <- scatter_diff(d_needle, "needle_pct",
                       "Field needleleaf share of tree cover (%)") +
    scatter_1to1(d_all) +
    scatter_diff(d_all, "rap", "RAP tree cover (%)") +
    scatter_diff(d_all, "field", "Field tree cover (%)") +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = paste0("RAP 2021 - field tree cover", zt_title[[zt_scatter]]),
      caption = paste("Field cover: mean over years at each pixel.",
                      "Top left: pixels with a field needleleaf fraction."))
  ggsave(file.path(fig_dir, paste0("rap_minus_field_tree_relationships_",
                                   zt_scatter, "_", field_vc, ".png")),
         comb, width = 11, height = 11, dpi = 300)
}
