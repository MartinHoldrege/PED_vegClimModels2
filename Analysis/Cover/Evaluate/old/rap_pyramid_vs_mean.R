# rap_pyramid_vs_mean.R
#
# How much does the way RAP is brought to 1 km matter? Compares, for 2021:
#   pyramid  RAP_v3_cover_2021_1000m.tif (03_rap_sample.js): GEE read RAP's own
#            integer pyramids, so values are rounded and come from ~960 m blocks
#            not aligned with the snap grid
#   mean     RAP_v3_cover_2021_from250m_1000m.tif (03_rap_cover_2021_250m_test.js):
#            mean of RAP over each 1 km snap cell (via a 250 m float asset)
#
# Statistics and plots use a random sample of cells with data in both.
# Difference = mean - pyramid (percentage points). "Edge" cells are kept
# cells with a masked or no-data cell among their 8 neighbours (e.g. next to
# cropland, developed land or a burn), where the pyramid footprint is most
# likely to take in land the masks exclude.
#
# Outputs (fig_dir): rap_pyramid_vs_mean_scatter.png, _hist.png,
#   _tree_map.png, and rap_pyramid_vs_mean_summary.csv
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

n_sample <- 1e6
fig_dir <- "Figures/Cover/rap_aggregation"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rap_dir <- file.path(paths$large, "Data_processed/CoverData/rap")
groups <- c("tree", "shrub", "herbaceous", "bare_ground")

# read --------------------------------------------------------------------

snap <- read_mask()
# GEE exports come back one cell larger in each direction than the snap grid
pyramid <- terra::rast(file.path(rap_dir, "RAP_v3_cover_2021_1000m.tif")) |>
  align_raster(snap)
mean_r <- terra::rast(file.path(rap_dir,
                                "RAP_v3_cover_2021_from250m_1000m.tif")) |>
  align_raster(snap)
stopifnot(all(names(pyramid) == groups), all(names(mean_r) == groups))

# kept cells with a masked / no-data neighbour
edge <- terra::focal(is.na(pyramid[["tree"]]), w = 3, fun = "max",
                     na.rm = TRUE) == 1
edge <- terra::mask(edge, pyramid[["tree"]])
names(edge) <- "edge"

# sample ------------------------------------------------------------------

stack <- c(setNames(pyramid, paste0("pyramid_", groups)),
           setNames(mean_r, paste0("mean_", groups)),
           edge)

set.seed(1)
smp <- terra::spatSample(stack, size = n_sample, method = "random",
                         na.rm = TRUE, as.df = TRUE) |>
  as_tibble()

long <- smp |>
  mutate(id = row_number()) |>
  pivot_longer(-c(id, edge), names_to = c("method", "group"),
               names_pattern = "^(pyramid|mean)_(.*)$",
               values_to = "cover") |>
  pivot_wider(names_from = method, values_from = cover) |>
  mutate(diff = mean - pyramid,
         group = factor(group, levels = groups),
         cells = if_else(edge == 1, "edge", "interior"))

# summary -----------------------------------------------------------------

summ <- long |>
  bind_rows(mutate(long, cells = "all")) |>
  summarise(n = n(),
            mean_diff = mean(diff),
            mean_abs_diff = mean(abs(diff)),
            rmsd = sqrt(mean(diff^2)),
            r = cor(mean, pyramid),
            p01 = quantile(diff, 0.01), p05 = quantile(diff, 0.05),
            p25 = quantile(diff, 0.25), p50 = median(diff),
            p75 = quantile(diff, 0.75), p95 = quantile(diff, 0.95),
            p99 = quantile(diff, 0.99),
            # share of values that aren't whole numbers (should be ~0 for
            # pyramid, ~1 for mean)
            frac_nonint_pyramid = mean(pyramid != round(pyramid)),
            frac_nonint_mean = mean(mean != round(mean)),
            .by = c(group, cells)) |>
  arrange(group, cells)

print(summ, n = Inf, width = Inf)
write_csv(summ, file.path(fig_dir, "rap_pyramid_vs_mean_summary.csv"))

# figures -----------------------------------------------------------------

# mean against pyramid, with 1:1
g_scatter <- ggplot(long, aes(pyramid, mean)) +
  geom_bin2d(bins = 100) +
  geom_abline(colour = "red") +
  scale_fill_viridis_c(trans = "log10", name = "Cells") +
  facet_wrap(~ group, scales = "free") +
  labs(x = "Pyramid (current, %)", y = "Mean of 1 km cell (%)",
       title = "RAP 2021 cover: pyramid vs exact 1 km mean",
       subtitle = paste0(format(n_sample, big.mark = ","),
                         " sampled cells; red: 1:1"))
ggsave(file.path(fig_dir, "rap_pyramid_vs_mean_scatter.png"), g_scatter,
       width = 10, height = 9, dpi = 300)

# differences, edge vs interior (central 98%, so the tails don't hide the
# body)
lims <- quantile(long$diff, c(0.01, 0.99))
g_hist <- ggplot(long,
                 aes(diff, after_stat(density), colour = cells)) +
  geom_freqpoly(bins = 200) +
  geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ group, scales = "free_y") +
  coord_cartesian(xlim = lims) +
  labs(x = "Mean - pyramid (% points)", y = "Density", colour = NULL,
       title = "Difference in RAP 2021 cover, edge vs interior cells",
       subtitle = "Edge: kept cell with a masked or no-data neighbour")
ggsave(file.path(fig_dir, "rap_pyramid_vs_mean_hist.png"), g_hist,
       width = 10, height = 7, dpi = 300)

# where the tree cover differences are
diff_tree <- mean_r[["tree"]] - pyramid[["tree"]]
lim <- max(abs(c(raster_quantile(diff_tree, 0.01),
                 raster_quantile(diff_tree, 0.99))))
g_map <- plot_map_conus(diff_tree,
                        colorscale = colorscale_diverging(
                          name = "Mean - pyramid\n(% points)",
                          limits = c(-lim, lim)),
                        title = "RAP 2021 tree cover: exact mean - pyramid",
                        inset_hist = TRUE)
ggsave(file.path(fig_dir, "rap_pyramid_vs_mean_tree_map.png"), g_map,
       width = 9, height = 6, dpi = 300)
