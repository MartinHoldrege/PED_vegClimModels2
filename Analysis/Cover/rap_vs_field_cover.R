# rap_vs_field_cover.R
#
# Compares RAP cover with field cover (the c04 field data), for each cover type
# both have: tree, shrub, herbaceous and bare ground. RAP is the 1 km annual
# RAP cover of the same year as the field data, so the comparison is by
# pixel-year (2000-2024; field years outside RAP's are dropped). RAP is LCMAP-
# and fire-masked, as in the training data: in a pixel that is mostly e.g.
# farmland, the 1 km RAP mean can't be compared with a field plot in a patch
# of natural vegetation.
#
# Sets of pixel-years, per cover type:
#   all        every field pixel-year with that cover type
#   nonforest  shrub, herbaceous and bare ground: only where RAP tree cover is
#              < tree_cutoff in that year. RAP can't see the understorey,
#              while field shrub and herbaceous cover include it
#   needle     if used applies to tree: only pixel-years with a field needleleaf fraction (the
#              original motivation: combined needleleaf cover, RAP-based tree
#              cover x modelled needleleaf share, was higher than observed)
#
# Per cover type, one PNG with a column per set:
#   field against RAP, with the 1:1 line (r, bias and n in the subtitle)
#   histograms of RAP - field, all sources and by source, with the mean
#     (bias) and median
#   map of RAP - field, averaged over each pixel's years
# Difference = RAP - field, in percentage points.
#
# Output (fig_dir): rap_vs_field_<cover type>_<field_vc>.png
#
# October, 2026 (rewritten November 2026: c04 field data, same-year annual
# RAP, every cover type)


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
library(patchwork)

# params ------------------------------------------------------------------

field_vc <- "c04"
years <- 2000:2024   # years of the annual RAP cover
tree_cutoff <- 10    # nonforest set: RAP tree cover < tree_cutoff (%)

fig_dir <- "Figures/Cover/rap_vs_field"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# data --------------------------------------------------------------------

field <- read_cover_training(field_vc) |>
  filter(year %in% years)
stopifnot(!any(field$sources == "rap"))  # field data only

# cover types in both the field data and RAP
cover_types <- c("tree", "shrub", "herbaceous", "bare_ground")
cover_types <- cover_types[paste0("cov_", cover_types) %in% names(field)]
message("cover types compared: ", paste(cover_types, collapse = ", "))
stopifnot("tree" %in% cover_types)  # needed for the nonforest set

# same-year RAP cover at each field pixel-year (rap_<type> columns), LCMAP-
# and fire-masked (NA where masked)
for (type in cover_types) {
  r <- read_rap_cover_raster(type, years = years)
  field[[paste0("rap_", type)]] <- rap_at(r, field)
}

# field source of each pixel-year; "mixed" if more than one
field$source <- if_else(str_detect(field$sources, ","), "mixed",
                        field$sources)

#' Pixel-years of one set for one cover type
#'
#' @param type Cover type, e.g. "shrub".
#' @param set "all", "nonforest" or "needle" (see header).
#' @return Rows of `field` with field_cover, rap_cover and diff (RAP - field).
#' @examples
#' type <- "shrub"
#' set <- "nonforest"
#' pixel_set(type = type, set = set)
pixel_set <- function(type, set) {
  d <- field |>
    mutate(field_cover = .data[[paste0("cov_", type)]],
           rap_cover = .data[[paste0("rap_", type)]]) |>
    filter(!is.na(field_cover), !is.na(rap_cover))
  if (set == "nonforest") d <- filter(d, rap_tree < tree_cutoff)
  if (set == "needle") d <- filter(d, !is.na(frac_needle))
  mutate(d, diff = rap_cover - field_cover)
}

# tree: all and needle; the other cover types: all and nonforest
combos <- bind_rows(
  tibble(type = "tree", set = c("all", "nonforest")), # or c('all', 'needle')
  expand_grid(type = setdiff(cover_types, "tree"), set = c("all", "nonforest"))
)
px <- pmap(combos, pixel_set)
print(mutate(combos, n_pixel_years = map_int(px, nrow)))

# figures -----------------------------------------------------------------

set_label <- c(all = "all field pixel-years",
               nonforest = paste0("RAP tree cover < ", tree_cutoff, "%"),
               needle = "pixel-years with a field needleleaf fraction")

#' Scatter, histograms and map for one cover type and set
#'
#' @param d Pixel-years, from pixel_set().
#' @param type Cover type.
#' @param set Name of the set.
#' @param lim Colour scale of the map runs from -lim to lim.
#' @return A patchwork.
plot_set <- function(d, type, set, lim) {
  type_label <- str_replace(type, "_", " ")

  # field against RAP
  g_scatter <- ggplot(d, aes(rap_cover, field_cover)) +
    geom_bin2d(bins = 50) +
    geom_abline(colour = "grey30") +
    scale_fill_viridis_c(trans = "log10", name = "Pixel-years") +
    coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(x = paste("RAP", type_label, "cover (%)"),
         y = paste("Field", type_label, "cover (%)"),
         title = paste0(str_to_sentence(type_label), " cover, ",
                        set_label[[set]]),
         subtitle = paste0("r = ", round(cor(d$rap_cover, d$field_cover), 2),
                           ", bias (RAP - field) = ", round(mean(d$diff), 1),
                           ", n = ", nrow(d), " pixel-years; grey: 1:1"))

  # histograms of RAP - field: all sources, then each source
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
    labs(x = paste0("RAP - field ", type_label, " cover (% points)"),
         y = "Pixel-years", subtitle = "Dashed: mean difference (bias)")

  # map of RAP - field, mean over each pixel's years
  g_map <- d |>
    summarise(diff = mean(diff), .by = c(cell, x, y)) |>
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(read_mask())) |>
    plot_points_conus(color_var = "diff", point_size = 0.15,
                      colorscale = colorscale_diverging(
                        name = "RAP - field\n(% points,\nmean over years)",
                        limits = c(-lim, lim),
                        aesthetics = "colour"),
                      title = NULL,
                      inset_hist = TRUE, inset_drop_zero = FALSE)
  
  wrap_plots(g_scatter, g_hist, g_map, ncol = 1, heights = c(1, 0.7, 1))
}

# one file per cover type, its sets side by side
for (type in cover_types) {
  rows <- which(combos$type == type)
  # one colour scale for the type's maps
  lim <- max(map_dbl(px[rows], \(d) max(abs(quantile(d$diff, c(0.1, 0.9))))))
  plots <- map(rows, \(i) plot_set(px[[i]], type, combos$set[i], lim))

  ggsave(file.path(fig_dir, paste0("rap_vs_field_", type, "_", field_vc,
                                   ".png")),
         wrap_plots(plots, nrow = 1), width = 20, height = 18, dpi = 300)
}