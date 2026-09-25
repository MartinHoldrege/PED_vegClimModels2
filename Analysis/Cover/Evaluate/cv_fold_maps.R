# plot_cv_folds.R
#
# Map of the environmental CV folds for one cover model: each training pixel
# coloured by the fold it is held out in. Folds are rebuilt from the spec
# (cluster variables, k, seed) exactly as in 02_fit_classification.R, so the
# model does not need to be fit first.
#
# Model chosen by --cover_type, --cover_model, --vc and --vmc (params.R).
#
# Inputs:
#   cover_clim_soils_<vc>.csv  - 08_combine_cover_and_covariates.R
#   cover_specs                - Functions/models/cover_specs.R
#
# Output:
#   Figures/Folds/folds_<cover_type>_<cover_model>_<vc>-<vmc>.png
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

# spec --------------------------------------------------------------------

spec <- cover_specs[[opt$cover_type]][[opt$cover_model]][[opt$vmc]]

if (is.null(spec)) {
  stop("no spec for cover_specs$", opt$cover_type, "$", opt$cover_model, "$",
       opt$vmc)
}

# data --------------------------------------------------------------------

# same rows as the fitting script, so the folds match
dat <- read_cover_training(opt$vc)

dat <- switch(
  spec$rows,
  all = dat,
  stop("unknown spec$rows: ", spec$rows)
)

source_vars <- unique(c(str_remove(spec$pred_vars, "^log1p_"),
                        spec$cv$cluster_vars))
dat <- dat[complete.cases(dat[c(spec$response, source_vars)]), ]

# folds -------------------------------------------------------------------

clusters <- make_env_clusters(dat,
                              vars = spec$cv$cluster_vars,
                              iter.max = 500,
                              k = spec$cv$k_clusters,
                              seed = 1,
                              group = dat$cell)

folds <- dat |>
  select(cell, x, y) |>
  mutate(fold = factor(clusters$env_cluster)) |>
  # one point per cell; group = cell already put all its years in one fold
  distinct(cell, .keep_all = TRUE) |>
  # random draw order, so no fold is systematically plotted on top
  slice_sample(prop = 1)

# colours -----------------------------------------------------------------

# evenly spaced hues at two lightness levels, shuffled so fold numbers that
# are close don't get similar colours. At large k (e.g. 50) some pairs will
# still look alike
k <- nlevels(folds$fold)
set.seed(1)
pal <- grDevices::hcl(h = seq(15, 375, length.out = k + 1)[-1],
                      c = 90,
                      l = rep(c(45, 75), length.out = k))[sample(k)]
names(pal) <- levels(folds$fold)

# map ---------------------------------------------------------------------

folds_sf <- sf::st_as_sf(folds, coords = c("x", "y"),
                         crs = terra::crs(read_mask()))

g <- plot_points_conus(folds_sf, color_var = "fold", point_size = 0.3,
                       title = paste0(opt$cover_type, "_", opt$cover_model,
                                      " ", opt$vmc, ": CV folds (k = ", k,
                                      ")"),
                       colorscale = scale_colour_manual(name = "Fold",
                                                        values = pal))

# legend only when it's readable
g <- if (k <= 12) {
  g + guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2))
} else {
  g + guides(colour = "none")
}

# fold numbers on the map, so folds can be told apart without colour. Each
# label goes in the fold's best clump: the 100 km bin (map units assumed m)
# with the most of its pixels, weighted by how much of the bin it makes up
if (k <= 10) {
  labels_sf <- folds |>
    mutate(bin_x = round(x / 1e5), bin_y = round(y / 1e5)) |>
    mutate(n_all = n(), .by = c(bin_x, bin_y)) |>
    mutate(n_fold = n(), .by = c(fold, bin_x, bin_y)) |>
    slice_max(n_fold^2 / n_all, n = 1, by = fold, with_ties = FALSE) |>
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(read_mask()))
  
  g <- g +
    geom_sf_label(data = labels_sf, aes(label = fold), size = 3,
                  label.padding = unit(0.15, "lines"))
}

# save --------------------------------------------------------------------

fig_file <- file.path("Figures", "Folds",
                      paste0("folds_", opt$cover_type, "_", opt$cover_model,
                             "_", opt$vc, "-", opt$vmc, ".png"))
dir.create(dirname(fig_file), recursive = TRUE, showWarnings = FALSE)
ggsave(fig_file, g, width = 9, height = 6, dpi = 600)

g