# cover_training_data_maps.R
#
# Maps of the cover training data (field plots + RAP augmentation) for tree,
# shrub, herbaceous and bare ground.
#
# Part 1: pixel maps. Pixel-years at the same pixel are averaged, so each map
# shows one value per 1 km pixel. One figure per group, with a panel for all
# sources combined and one per source, to show each source's extent.
#
# Pixels are drawn as square markers at the cell centre, not as a raster: at
# CONUS scale a 1 km cell is about one device pixel, so isolated cells would
# not be visible. Markers are larger than a cell and overlap where sampling
# is dense, so draw order is randomised to keep high or low values from being
# systematically on top.
#
# Single-source panels use only pixel-years with that source alone. Pixel-years
# with more than one field source (e.g. "fia,ldc") hold their average, so they
# get their own panel rather than being attributed to one source.
#
# Each panel has an inset histogram of its non-zero pixel means, scaled to the
# panel's own tallest bar, so shapes are comparable but counts are not (counts
# are in the panel titles). The share of zeros is printed above the inset.
#
# Part 2: EPA Level III ecoregion summaries, all sources combined (described
# at the start of Part 2).
#
# Inputs:
#   cover_by_pixel_year_all-sources_<vc>.csv - 06_cover_add-rap.R
#   ecoregion_sampling_gap.csv               - 06_cover_add-rap.R
#   EPA_L3_ecoregion_daymet_1000m.tif (load_ecoregion_raster())
#                                            - 01_rasterize_ecoregions.R
#   LCMAP_fracKeep_gte90_1000m.tif           - 03_export_masks.js, via
#                                              04_download_gee_output.R
#   MTBS_fracUnburnedMean_gte90_20yr_2011-2023_1000m.tif
#                                            - 03_export_masks.js, via
#                                              04_download_gee_output.R
#   daymet_conus_snap_1000m.tif              - 00_create_snap_raster.R
#                                              (read_mask())
#
# Outputs (Figures/Cover/training_data/):
#   cover_training_<group>_<vc>.png
#   ecoregion/eco_cover_<vc>.png
#   ecoregion/eco_proportions_<vc>.png
#   ecoregion/eco_density_<vc>.png
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
library(patchwork)
source_functions()
# plot_points_conus(), plot_map_conus(), colorscale_cover(),
# colorscale_biomass() from mapping.R; crs_daymet, align_raster(), read_mask()
# from spatial.R; mean_na() from general.R; load_ecoregion_raster()

vc <- opt$vc  # cover data version

# params ------------------------------------------------------------------

# test_run = TRUE draws at most 2,000 pixels per panel, at low dpi, written to
# _test files. For checking layout and labels without the full render.
test_run <- FALSE

in_file <- file.path(paths$large, "Data_processed/cover_combined",
                     paste0("cover_by_pixel_year_all-sources_", vc, ".csv"))

out_dir <- file.path("Figures/Cover/training_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c(tree        = "cov_tree",
            shrub       = "cov_shrub",
            herbaceous  = "cov_herbaceous",
            bare_ground = "cov_bare_ground")

# panel order and labels; names other than "all" and "multiple" are the
# tokens in `sources`
source_labels <- c(all = "All sources", rap = "RAP", fia = "FIA",
                   ldc = "LDC", lfrdb = "LFRDB",
                   multiple = "Multiple field sources")

trunc_prob   <- 0.99  # upper colour limit: quantile of all-source pixel means
point_size   <- 0.12  # marker size; tune so isolated pixels are visible
point_stroke <- 0     # ggplot's default (0.5) adds ~1 pt to every marker
pal_end      <- 0.8   # drop mako's lightest colours, invisible on white
zero_color   <- "grey70"
fig_width    <- 10
panel_height <- 3.2   # figure height = panel_height * rows + 1
dpi          <- if (test_run) 100 else 600

set.seed(1234)  # draw order

# Part 1: pixel maps ------------------------------------------------------

# read --------------------------------------------------------------------

cover <- read_csv(in_file, show_col_types = FALSE)

stopifnot(all(c("cell", "x", "y", "sources", groups) %in% names(cover)))

src_tokens <- unique(unlist(str_split(unique(cover$sources), ",")))
stopifnot(all(src_tokens %in% names(source_labels)))

message("pixel-years by source:")
count(cover, sources, sort = TRUE) |> as.data.frame() |> print()

# pixel means -------------------------------------------------------------
# Every pixel-year goes in "all", and once more in its source panel.

pixel_means <- cover |>
  mutate(panel = if_else(str_detect(sources, ","), "multiple", sources)) |>
  bind_rows(mutate(cover, panel = "all")) |>
  summarise(across(all_of(unname(groups)), mean_na),
            .by = c(panel, cell, x, y))

message("\npixels with data, by panel and group:")
pixel_means |>
  summarise(across(all_of(unname(groups)), \(z) sum(!is.na(z))), .by = panel) |>
  as.data.frame() |> print()

# maps --------------------------------------------------------------------

for (g in names(groups)) {
  
  d <- pixel_means |>
    filter(!is.na(.data[[groups[[g]]]])) |>
    select(panel, x, y, cover = all_of(groups[[g]]))
  
  # shared across panels, so sources are comparable
  upper <- quantile(d$cover[d$panel == "all"], trunc_prob, names = FALSE)
  
  # facet label with pixel count, in source_labels order
  d <- d |>
    mutate(n = n(), .by = panel) |>
    mutate(panel_label = paste0(source_labels[panel], " (",
                                scales::comma(n), " pixels)"),
           panel_label = fct_reorder(panel_label,
                                     match(panel, names(source_labels))))
  
  if (test_run) d <- slice_sample(d, n = 2000, by = panel)
  
  d <- slice_sample(d, prop = 1)
  
  d_sf <- sf::st_as_sf(d, coords = c("x", "y"), crs = crs(read_mask()))
  
  g_label <- str_to_sentence(str_replace(g, "_", " "))
  
  fig <- plot_points_conus(
    d_sf,
    color_var    = "cover",
    point_size   = point_size,
    point_shape  = 15,
    point_stroke = point_stroke,
    colorscale   = colorscale_cover(name = "Cover (%)",
                                    limits = c(0, upper),
                                    zero_color = zero_color,
                                    pal_end = pal_end,
                                    aesthetics = "colour"),
    facet_var    = "panel_label",
    title        = paste(g_label, "cover, training data"),
    inset_hist   = TRUE
  ) +
    labs(subtitle = paste0("Mean across years at each 1 km pixel. ",
                           "Colour scale truncated at the ",
                           100 * trunc_prob, "th percentile (",
                           round(upper, 1), "%).\nInsets: histogram of the ",
                           "panel's non-zero pixel means, scaled to its ",
                           "tallest bar, with the share of zeros above."))
  
  n_rows <- ceiling(n_distinct(d$panel) / 2)
  
  out_file <- file.path(out_dir, paste0("cover_training_", g, "_", vc,
                                        if (test_run) "_test", ".png"))
  ggsave(out_file, fig, width = fig_width,
         height = panel_height * n_rows + 1, dpi = dpi, bg = "white")
  message("Wrote ", out_file)
}


# Part 2: ecoregion summaries ---------------------------------------------
# All sources combined. Three figures:
#   eco_cover       - mean cover per EPA Level III ecoregion, by group
#   eco_proportions - needle/broad shares of tree, forb/C3/C4 shares of
#                     herbaceous
#   eco_density     - total and unmasked area per sampled pixel, by group
#
# Means are drawn only on unmasked cells (LCMAP and 2011-2023 mean fire
# masks, i.e. the available area in 06_cover_add-rap.R), so masked land is not
# coloured as if it were represented. Density is drawn over the whole
# ecoregion. Sampled pixels include any in masked cells (FIA plots, which skip
# the LCMAP mask, and plots measured before a later fire), as in the density
# check at the end of 06.
#
# Cover means are means of pixel means (repeat years averaged first), so each
# pixel counts once. Proportions follow 05_combine_field_cover.R: components
# are averaged over years, then over pixels, and the proportion is taken once.
# That weights pixels by their cover, and pixels with none drop out. The
# herbaceous proportions use only pixel-years with a C3/C4 split, so FIA-only
# pixel-years are excluded. RAP has no components.

# * params ----------------------------------------------------------------

mask_dir    <- file.path(paths$large, "Data_processed/masks")
eco_out_dir <- file.path(out_dir, "ecoregion")
dir.create(eco_out_dir, recursive = TRUE, showWarnings = FALSE)

min_pixels <- 10     # ecoregion means from fewer pixels are not shown
target_km2 <- 24.28  # target unmasked area per sampled pixel, as in 06
# 4381 x 2733 cells; 5e6 is near-native for these panel sizes
maxcell    <- if (test_run) 1e4 else 5e6

comp_cols <- c("c_needle", "c_broad", "c_forb", "c_c3", "c_c4")

frac_labels <- c(frac_needle = "Needleleaf share of tree",
                 frac_broad  = "Broadleaf share of tree",
                 frac_forb   = "Forb share of herbaceous",
                 frac_c3     = "C3 grass share of herbaceous",
                 frac_c4     = "C4 grass share of herbaceous")

stopifnot(all(comp_cols %in% names(cover)))

# * read ------------------------------------------------------------------

snap <- read_mask()
eco  <- load_ecoregion_raster("L3")  # values are eco_id (`region` in 06)

lcmap_mask <- terra::rast(file.path(mask_dir,
                                    "LCMAP_fracKeep_gte90_1000m.tif")) |>
  align_raster(snap)
fire_mask  <- terra::rast(file.path(
  mask_dir, "MTBS_fracUnburnedMean_gte90_20yr_2011-2023_1000m.tif")) |>
  align_raster(snap)

# ecoregion id on unmasked cells, NA elsewhere
eco_avail <- terra::ifel(lcmap_mask == 1 & fire_mask == 1, eco, NA)

# km2 per ecoregion, as defined in 06
eco_area <- read_csv(file.path(paths$large, "Data_processed/cover_combined",
                               "ecoregion_sampling_gap.csv"),
                     show_col_types = FALSE) |>
  distinct(region, area_total, area_available)
stopifnot(!anyDuplicated(eco_area$region))

# * helpers ---------------------------------------------------------------

#' Paint one value per ecoregion onto an ecoregion raster
#'
#' @param r_eco SpatRaster of ecoregion ids.
#' @param region Ecoregion ids.
#' @param value One value per id. Ids not given, or given NA, become NA.
#' @return SpatRaster of ecoregion values.
eco_to_rast <- function(r_eco, region, value) {
  stopifnot(length(region) == length(value), !anyDuplicated(region))
  terra::subst(r_eco, from = region, to = value, others = NA)
}

#' Map one ecoregion-level variable
#'
#' @param d Data frame with `region` and `value`.
#' @param colorscale ggplot2 fill scale.
#' @param title Panel title.
#' @param r_eco Ecoregion raster to paint onto; `eco_avail` for unmasked cells
#'   only, `eco` for whole ecoregions.
plot_eco <- function(d, colorscale, title, r_eco = eco_avail) {
  r <- eco_to_rast(r_eco, d$region, d$value)
  plot_map_conus(r, colorscale = colorscale, title = title, maxcell = maxcell)
}

#' Share of a total, NA where the total is zero or missing
share <- function(part, total) if_else(total > 0, part / total, NA_real_)

eco_file <- function(name) {
  file.path(eco_out_dir, paste0(name, "_", vc, if (test_run) "_test", ".png"))
}

group_label <- function(g) str_to_sentence(str_replace(g, "_", " "))

# * pixel means with ecoregion --------------------------------------------

px_all <- pixel_means |>
  filter(panel == "all") |>
  select(-panel) |>
  mutate(region = eco[cell][[1]])
stopifnot(!anyNA(px_all$region))

# * cover means -----------------------------------------------------------

eco_cover <- px_all |>
  pivot_longer(all_of(unname(groups)), names_to = "col", values_to = "value") |>
  filter(!is.na(value)) |>
  summarise(value = mean(value), n_pixels = n(), .by = c(region, col)) |>
  mutate(group = names(groups)[match(col, groups)])

# * proportions -----------------------------------------------------------

comp_px <- cover |>
  mutate(
    # tree pair only where both members are present
    tree_ok  = !is.na(c_needle) & !is.na(c_broad),
    c_needle = if_else(tree_ok, c_needle, NA_real_),
    c_broad  = if_else(tree_ok, c_broad,  NA_real_),
    # herbaceous split only where C3/C4 exists; a missing member is zero, as
    # in the herbaceous denominator in 05
    split = !is.na(c_c3) | !is.na(c_c4),
    across(c(c_forb, c_c3, c_c4),
           \(z) if_else(split, coalesce(z, 0), NA_real_))
  ) |>
  summarise(across(all_of(comp_cols), mean_na), .by = cell) |>
  mutate(region = eco[cell][[1]])

eco_frac_wide <- comp_px |>
  summarise(n_tree = sum(c_needle + c_broad > 0, na.rm = TRUE),
            n_herb = sum(c_forb + c_c3 + c_c4 > 0, na.rm = TRUE),
            across(all_of(comp_cols), \(z) mean(z, na.rm = TRUE)),
            .by = region) |>
  mutate(tree_total  = c_needle + c_broad,
         herb_total  = c_forb + c_c3 + c_c4,
         frac_needle = share(c_needle, tree_total),
         frac_broad  = share(c_broad,  tree_total),
         frac_forb   = share(c_forb,   herb_total),
         frac_c3     = share(c_c3,     herb_total),
         frac_c4     = share(c_c4,     herb_total))

# shares within a set sum to one wherever defined
s_tree <- with(eco_frac_wide, frac_needle + frac_broad)
s_herb <- with(eco_frac_wide, frac_forb + frac_c3 + frac_c4)
stopifnot(all(abs(s_tree[!is.na(s_tree)] - 1) < 1e-6),
          all(abs(s_herb[!is.na(s_herb)] - 1) < 1e-6))

eco_frac <- eco_frac_wide |>
  pivot_longer(all_of(names(frac_labels)), names_to = "var",
               values_to = "value") |>
  mutate(n_pixels = if_else(var %in% c("frac_needle", "frac_broad"),
                            n_tree, n_herb)) |>
  select(region, var, value, n_pixels)

# * minimum sample --------------------------------------------------------

eco_cover <- eco_cover |>
  mutate(value = if_else(n_pixels >= min_pixels, value, NA_real_))
eco_frac <- eco_frac |>
  mutate(value = if_else(n_pixels >= min_pixels, value, NA_real_))

message("\necoregions with a mean shown (of ", nrow(eco_area), "):")
bind_rows(select(eco_cover, var = group, value),
          select(eco_frac, var, value)) |>
  summarise(n_shown = sum(!is.na(value)), .by = var) |>
  as.data.frame() |> print()

# * figure: cover ---------------------------------------------------------

eco_caption <- paste0(
  "Mean of pixel means (repeat years averaged), all sources. Shown only on ",
  "cells passing the LCMAP and 2011-2023 mean fire masks.\nEcoregions with ",
  "fewer than ", min_pixels, " sampled pixels are not shown.")

cover_panels <- map(names(groups), \(g) {
  d <- filter(eco_cover, group == g)
  plot_eco(d,
           colorscale = colorscale_cover(
             name = "Cover (%)",
             limits = c(0, max(d$value, na.rm = TRUE)),
             pal_end = pal_end),
           title = group_label(g))
})

fig_cover <- wrap_plots(cover_panels, ncol = 2) +
  plot_annotation(title = "Mean cover by EPA Level III ecoregion, training data",
                  caption = eco_caption)

ggsave(eco_file("eco_cover"), fig_cover, width = 12, height = 7.5,
       dpi = dpi, bg = "white")
message("Wrote ", eco_file("eco_cover"))

# * figure: proportions ---------------------------------------------------

frac_scale <- colorscale_biomass(name = "Proportion", limits = c(0, 1))

frac_panels <- imap(frac_labels, \(lab, v) {
  plot_eco(filter(eco_frac, var == v), frac_scale, lab)
})

# legend goes in the empty cell after the tree pair
fig_frac <- wrap_plots(frac_panels$frac_needle, frac_panels$frac_broad,
                       guide_area(),
                       frac_panels$frac_forb, frac_panels$frac_c3,
                       frac_panels$frac_c4,
                       ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Mean proportions by EPA Level III ecoregion, training data",
    caption = paste0(
      "Components averaged over years, then pixels, and the proportion ",
      "taken once (cover-weighted). Field sources only; herbaceous shares ",
      "exclude FIA-only pixel-years (no C3/C4 split).\n",
      "Shown only on cells passing the LCMAP and 2011-2023 mean fire masks. ",
      "Ecoregions with fewer than ", min_pixels,
      " pixels with cover in the set are not shown."))

ggsave(eco_file("eco_proportions"), fig_frac, width = 15, height = 6.5,
       dpi = dpi, bg = "white")
message("Wrote ", eco_file("eco_proportions"))

# * figure: sampling density ----------------------------------------------

eco_density <- expand_grid(eco_area, group = names(groups)) |>
  left_join(select(eco_cover, region, group, n_pixels),
            by = c("region", "group")) |>
  mutate(n_pixels      = replace_na(n_pixels, 0),
         km2_total     = if_else(n_pixels > 0,
                                 area_total / n_pixels, NA_real_),
         km2_available = if_else(n_pixels > 0 & area_available > 0,
                                 area_available / n_pixels, NA_real_))

message("\necoregions with no sampled pixels, by group:")
eco_density |>
  summarise(n_none = sum(n_pixels == 0), .by = group) |>
  as.data.frame() |> print()

# one log scale across all panels
dens_lim   <- range(c(eco_density$km2_total, eco_density$km2_available),
                    na.rm = TRUE)
dens_scale <- colorscale_biomass(name = "km\u00b2 per\nsampled pixel",
                                 limits = dens_lim, transform = "log10")

dens_panels <- map(names(groups), \(g) {
  d <- filter(eco_density, group == g)
  list(
    plot_eco(transmute(d, region, value = km2_total), dens_scale,
             paste0(group_label(g), ": total area / sampled pixels"),
             r_eco = eco),
    plot_eco(transmute(d, region, value = km2_available), dens_scale,
             paste0(group_label(g), ": unmasked area / sampled pixels"),
             r_eco = eco)
  )
}) |>
  list_flatten()

fig_density <- wrap_plots(dens_panels, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Sampling density by EPA Level III ecoregion, training data",
    caption = paste0(
      "All sources. Unmasked area passes the LCMAP and 2011-2023 mean fire ",
      "masks; the RAP augmentation target was ", target_km2,
      " km\u00b2 of unmasked area per pixel.\nSampled pixels include any in ",
      "masked cells. Blank ecoregions have no sampled pixels."))

ggsave(eco_file("eco_density"), fig_density, width = 11, height = 13,
       dpi = dpi, bg = "white")
message("Wrote ", eco_file("eco_density"))