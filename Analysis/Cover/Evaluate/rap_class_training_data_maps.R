# class_training_data_maps.R
#
# Maps of the two classifications the models are fit to, from the 2021 RAP
# data:
#   forest:    RAP tree cover > 10%
#   zero tree: > 90% of a cell's natural-land 30 m pixels have < 3% RAP tree
#              cover
# One figure per class, with two panels: the wall-to-wall raster (every cell
# passing the LCMAP and fire masks) and the training sample points (c02,
# 05_rap_training_sample.R). Classes use the same rule as
# 02_fit_classification.R (value > threshold).
#
# Sample points are drawn as square markers at the cell centre, as in
# cover_training_data_maps.R, in random order so neither class is
# systematically on top.
#
# Inputs:
#   RAP_v3_cover_2021_1000m.tif   - 03_rap_sample.js (04_download)
#   RAP_v3_fracZeroTree_lt3_2021-2021_1000m.tif (read_zero_tree_raster())
#                                 - 03_rap_frac-zero-tree.js
#   cover_clim_soils_c02.csv (read_cover_training())
#                                 - 05_rap_training_sample.R
#   daymet_conus_snap_1000m.tif (read_mask())
#                                 - 00_create_snap_raster.R
#
# Outputs (Figures/Cover/training_data/):
#   class_forest_<vc>.png
#   class_zero_tree_<vc>.png
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
library(patchwork)
source_functions()

# params ------------------------------------------------------------------

# test_run = TRUE: coarse raster and at most 2,000 points, at low dpi, written
# to _test files. For checking layout and labels without the full render.
test_run <- FALSE

vc <- "c02"  # RAP sample; the only version with the zero-tree columns

rap_file <- file.path(paths$large, "Data_processed/CoverData/rap",
                      "RAP_v3_cover_2021_1000m.tif")

out_dir <- file.path("Figures/Cover/training_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# one row per class: column in the training data, threshold (%), colour of
# the class (the other class is grey)
classes <- tibble::tribble(
  ~name,       ~label,      ~column,         ~threshold, ~colour,
  "forest",    "forest",    "cov_tree",      10,         "darkgreen",
  "zero_tree", "zero tree", "pct_zero_tree", 90,         "tan3"
)

other_colour <- "grey85"
point_size   <- 0.12  # as in cover_training_data_maps.R
fig_width    <- 14
fig_height   <- 5
dpi          <- if (test_run) 100 else 800
# 4381 x 2733 cells; 5e6 is near-native for these panel sizes
maxcell      <- if (test_run) 1e4 else 1e7

set.seed(1234)  # draw order

# read --------------------------------------------------------------------

snap <- read_mask()

# % per cell: tree cover, and % of natural land with no trees. GEE exports
# come back one cell larger in each direction than the snap grid
pct_r <- list(
  forest    = terra::rast(rap_file)[["tree"]],
  zero_tree = 100 * read_zero_tree_raster("fracZeroTree")
) |>
  map(\(r) align_raster(r, snap))

if (test_run) {
  pct_r <- map(pct_r, \(r) terra::spatSample(r, size = 1e4, method = "regular",
                                             as.raster = TRUE))
}

pts <- read_cover_training(vc)
stopifnot(all(classes$column %in% names(pts)),
          !anyDuplicated(pts$cell))  # one year per cell in c02

# maps --------------------------------------------------------------------

pwalk(classes, \(name, label, column, threshold, colour) {

  class_names <- c(paste("not", label), label)  # 0, 1
  colours <- setNames(c(other_colour, colour), class_names)

  # * raster ----

  r <- terra::as.int(pct_r[[name]] > threshold)
  share_r <- terra::global(r, "mean", na.rm = TRUE)[[1]]
  levels(r) <- data.frame(value = c(0, 1), class = class_names)

  g_rast <- plot_map_conus(
    r,
    colorscale = scale_fill_manual(name = NULL, values = colours,
                                   na.translate = FALSE),
    title = paste0("All cells (", round(100 * share_r, 1), "% ", label, ")"),
    maxcell = maxcell
  )

  # * sample points ----

  d <- pts |>
    filter(!is.na(.data[[column]])) |>
    mutate(class = factor(if_else(.data[[column]] > threshold,
                                  class_names[2], class_names[1]),
                          levels = class_names))
  share_p <- mean(d$class == class_names[2])

  if (test_run) d <- slice_sample(d, n = min(2000, nrow(d)))
  d <- slice_sample(d, prop = 1)  # random draw order

  d_sf <- sf::st_as_sf(d, coords = c("x", "y"), crs = terra::crs(snap))

  g_pts <- plot_points_conus(
    d_sf,
    color_var    = "class",
    point_size   = point_size,
    point_shape  = 15,
    point_stroke = 0,
    colorscale   = scale_colour_manual(name = NULL, values = colours),
    title        = paste0("Training sample, ", vc, " (",
                          scales::comma(nrow(d)), " cells; ",
                          round(100 * share_p, 1), "% ", label, ")")
  ) +
    # larger markers in the legend only
    guides(colour = guide_legend(override.aes = list(size = 3)))

  # * assemble and write ----

  rule <- if (name == "forest") {
    paste0("RAP tree cover > ", threshold, "%")
  } else {
    paste0("> ", threshold, "% of natural-land 30 m pixels with < 3% RAP ",
           "tree cover")
  }

  fig <- g_rast + g_pts +
    plot_annotation(
      title = paste0(str_to_sentence(label), " vs not ", label, ", 2021"),
      subtitle = paste0(str_to_sentence(label), ": ", rule,
                        ". Cells passing the LCMAP and fire masks.")
    )

  out_file <- file.path(out_dir, paste0("class_", name, "_", vc,
                                        if (test_run) "_test", ".png"))
  ggsave(out_file, fig, width = fig_width, height = fig_height, dpi = dpi,
         bg = "white")
  message("Wrote ", out_file)
})

# combined map ------------------------------------------------------------
# Both classes on the raster, to show where the forest and zero-tree layers
# agree: forest and zero tree should rarely overlap (red), and non-forest
# should split into cells with some trees (light green) and none (tan). Cells
# missing either layer (their masks differ) are left out.

combo_names <- c("not forest, not zero tree", "not forest, zero tree",
                 "forest, not zero tree", "forest, zero tree")
combo_colours <- setNames(c("lightgreen", "tan3", "darkgreen", "blue"),
                          combo_names)

th <- setNames(classes$threshold, classes$name)

# 0-3 = 2 * forest + zero tree, the order of combo_names
combo_r <- 2 * (pct_r$forest > th[["forest"]]) +
  (pct_r$zero_tree > th[["zero_tree"]])
combo_r <- terra::as.int(combo_r)
levels(combo_r) <- data.frame(value = 0:3, class = combo_names)

# % of cells in each class:
terra::freq(combo_r) |>
  mutate(percent = 100 * count / sum(count)) |>
  select(class = value, count, percent) |>
  as.data.frame() |> print(digits = 3)

fig <- plot_map_conus(
  combo_r,
  colorscale = scale_fill_manual(name = NULL, values = combo_colours,
                                 na.translate = FALSE),
  title = "Forest and zero-tree classes combined, 2021",
  maxcell = maxcell
) +
  labs(subtitle = paste0("Forest: RAP tree cover > ", th[["forest"]], "%. ",
                         "Zero tree: > ", th[["zero_tree"]], "% of ",
                         "natural-land 30 m pixels with < 3% RAP tree cover.",
                         "\nCells with both layers, passing the LCMAP and ",
                         "fire masks."))

out_file <- file.path(out_dir, paste0("class_forest_zero_tree_", vc,
                                      if (test_run) "_test", ".png"))
# one panel, so half the width of the two-panel figures
ggsave(out_file, fig, width = fig_width / 2 + 1, height = fig_height,
       dpi = dpi, bg = "white")
