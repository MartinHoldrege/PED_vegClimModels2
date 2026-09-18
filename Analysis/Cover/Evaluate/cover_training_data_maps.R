# cover_training_data_maps.R
#
# Maps of the cover training data (field plots + RAP augmentation) for tree,
# shrub, herbaceous and bare ground. Pixel-years at the same pixel are
# averaged, so each map shows one value per 1 km pixel. One figure per group,
# with a panel for all sources combined and one per source, to show each
# source's extent.
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
# Input:
#   cover_by_pixel_year_all-sources_<vc>.csv - 06_cover_add-rap.R
#
# Output:
#   Figures/Cover/training_data/cover_training_<group>_<vc>.png
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
# plot_points_conus(), colorscale_cover() from mapping.R; crs_daymet from
# spatial.R; mean_na() from general.R

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
    title        = paste(g_label, "cover, training data")
  ) +
    labs(subtitle = paste0("Mean across years at each 1 km pixel. ",
                           "Colour scale truncated at the ",
                           100 * trunc_prob, "th percentile (",
                           round(upper, 1), "%)."))
  
  n_rows <- ceiling(n_distinct(d$panel) / 2)
  
  out_file <- file.path(out_dir, paste0("cover_training_", g, "_", vc,
                                        if (test_run) "_test", ".png"))
  ggsave(out_file, fig, width = fig_width,
         height = panel_height * n_rows + 1, dpi = dpi, bg = "white")
  message("Wrote ", out_file)
}