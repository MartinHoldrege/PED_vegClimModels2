# ldc_explore_hit_definitions.R
#
# Exploratory figures comparing any-hit and first-hit cover for the three
# top-level classes. First-hit only sees the topmost layer at a pin, so cover
# beneath a canopy is invisible; any-hit retains it. The difference should be
# largest where there is something above to hide under, which is what the
# second figure tests against tree cover.
#
# Input:  cover_by_plot_year.csv, created in 04_ldc_process.R
# Output: figures in Figures/ldc_hit_definitions/
#
# September, 2026

source("Functions/init.R")

ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")
fig_dir <- "Figures/ldc_hit_definitions"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cover1 <- read_csv(file.path(ldc_processed_dir, "cover_by_plot_year.csv"),
                   show_col_types = FALSE)

groups <- c("tree", "shrub", "herbaceous")

# ---- long format: one row per plot-year per group ------------------------


paired <- cover1 |>
  select(ProjectKey, PlotID, Longitude_NAD83, Latitude_NAD83, year, ah_tree,
         all_of(c(paste0("fh_", groups), paste0("ah_", groups)))) |>
  pivot_longer(matches("^(fh|ah)_(tree|shrub|herbaceous)$"),
               names_to = c("hit", "group"),
               names_pattern = "^(fh|ah)_(.*)$",
               values_to = "cover") |>
  pivot_wider(names_from = hit, values_from = cover) |>
  mutate(diff = ah - fh,
         group = factor(group, levels = groups))

paired2 <- paired |>
  left_join(select(cover1, ProjectKey, PlotID, Longitude_NAD83, 
                   Latitude_NAD83, year, fh_tree, ah_tree, fh_shrub, ah_shrub),
            by = join_by(ProjectKey, PlotID, Longitude_NAD83, Latitude_NAD83, year))

# ---- 1. any-hit vs first-hit --------------------------------------------
p_scatter <- ggplot(paired, aes(fh, ah)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  facet_wrap(~ group) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "First-hit cover (%)", y = "Any-hit cover (%)",
       title = "Any-hit vs first-hit cover",
       subtitle = "Points above the 1:1 line are cover hidden beneath a higher layer") +
  theme_bw()

ggsave(file.path(fig_dir, "ah_vs_fh_scatter.png"), p_scatter,
       width = 10, height = 4, dpi = 300)

# ---- 2. difference against tree cover ------------------------------------
# Tree cover is the candidate explanation for understorey being hidden, so
# only the two understorey classes are shown. Any-hit tree cover is used as
# the predictor since it is the more complete measure of tree presence.
p_diff <- paired2 |>
  filter(group != "tree") |>
  ggplot(aes(fh_tree, diff)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              colour = "red", linewidth = 0.6, se = FALSE) +
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  facet_wrap(~ group) +
  labs(x = "First hit tree cover (%)",
       y = "Any-hit minus first-hit cover (percentage points)",
       title = "Cover hidden by first-hit, against tree cover",
       subtitle = "If trees hide understorey, the difference should rise with tree cover") +
  theme_bw()

ggsave(file.path(fig_dir, "diff_vs_tree_cover.png"), p_diff,
       width = 9, height = 4, dpi = 300)

# ---- 3. herbaceous difference against shrub cover -------------------------
# Shrubs hide herbaceous cover by the same mechanism as trees, one layer
# down, and are far more widespread in these plots.
p_diff_shrub <- paired2 |>
  filter(group == "herbaceous") |>
  ggplot(aes(fh_shrub, diff)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              colour = "red", linewidth = 0.6, se = FALSE) +
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  labs(x = "First-hit shrub cover (%)",
       y = "Any-hit minus first-hit herbaceous cover (percentage points)",
       title = "Herbaceous cover hidden by first-hit, against shrub cover") +
  theme_bw()

ggsave(file.path(fig_dir, "herb_diff_vs_shrub_cover.png"), p_diff_shrub,
       width = 5.5, height = 4, dpi = 300)

# ---- 4. herbaceous difference against combined overstorey ----------------
# Trees and shrubs both sit above the herbaceous layer, so the quantity that
# should predict hidden herbaceous cover is whatever is above it.
p_diff_over <- paired2 |>
  filter(group == "herbaceous") |>
  mutate(fh_over = fh_tree + fh_shrub) |>
  ggplot(aes(fh_over, diff)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              colour = "red", linewidth = 0.6, se = FALSE) +
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  labs(x = "first-hit tree + shrub cover",
       y = "Any-hit minus first-hit herbaceous cover (percentage points)",
       title = "Herbaceous cover hidden by first-hit, against overstory cover") +
  theme_bw()

ggsave(file.path(fig_dir, "herb_diff_vs_overstory.png"), p_diff_over,
       width = 5.5, height = 4, dpi = 300)

# ---- summary numbers -----------------------------------------------------
paired |>
  summarise(median_fh   = round(median(fh), 2),
            median_ah   = round(median(ah), 2),
            median_diff = round(median(diff), 2),
            q90_diff    = round(quantile(diff, 0.9), 2),
            pct_ah_gt_fh = round(100 * mean(diff > 0), 1),
            .by = group) |>
  as.data.frame() |>
  print()

