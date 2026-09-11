# ldc_explore_hit_definitions.R
#
# Exploratory figures comparing any-hit and first-hit cover for the three
# top-level classes. First-hit only sees the topmost layer at a pin, so cover
# beneath a canopy is invisible; any-hit retains it. The difference should be
# largest where there is something above to hide under, which is what the
# second figure tests against tree cover. also compare to fht--first hit,
# unless first hit is tree, then the '2nd' hit for non trees is also counted
# (i.e. understory is seperate from overstory, if tree overstory exists)
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
         all_of(c(paste0("fh_", groups), paste0("ah_", groups), 
                  paste0('fht_', groups)))) |>
  pivot_longer(matches("^(fh|ah|fht)_(tree|shrub|herbaceous)$"),
               names_to = c("hit", "group"),
               names_pattern = "^(fh|ah|fht)_(.*)$",
               values_to = "cover") |>
  pivot_wider(names_from = hit, values_from = cover) |>
  mutate(diff = ah - fh,
         diff_fht = fht - fh,
         group = factor(group, levels = groups))

paired2 <- paired |>
  left_join(select(cover1, ProjectKey, PlotID, Longitude_NAD83, 
                   Latitude_NAD83, year, matches('_(tree|shrub)$')),
            by = join_by(ProjectKey, PlotID, Longitude_NAD83, Latitude_NAD83, year))


# fig params --------------------------------------------------------------

smooth <- function() {
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                        colour = "red", linewidth = 0.6, se = FALSE)
}

caption = 'FH: first pin hit; AH: any pin hit; FHT: FH, unless FH is tree, then shrub/herbaceous is FH after trees)'

base <- function() {
  list(geom_hex(bins = 100),
       coord_equal(xlim = c(0, 100), ylim = c(0, 100)),
       labs(
         caption = caption
       ))
}

# ---- 1. scatterplot vs. first-hit --------------------------------------------

p_scatter <- paired |>
  pivot_longer(c(fht, ah), names_to = "definition", values_to = "cover") |>
  mutate(definition = factor(definition, c("fht", "ah"),
                             c("FHT (through trees)", "AH (any hit)"))) |>
  ggplot(aes(fh, cover)) +
  base() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  facet_grid(definition ~ group) +
  labs(x = "FH cover (%)", y = "Cover under alternative definition (%)",
       title = "How far each definition departs from first hit")

ggsave(file.path(fig_dir, "definitions_vs_fh_scatter.png"), p_scatter,
       width = 10, height = 7, dpi = 300)

# ---- 2. difference against tree cover ------------------------------------
# Tree cover is the candidate explanation for understorey being hidden, so
# only the two understorey classes are shown. Any-hit tree cover is used as
# the predictor since it is the more complete measure of tree presence.
p_diff <- paired2 |>
  filter(group != "tree") |>
  ggplot(aes(fh_tree, diff)) +
  base() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  smooth()+
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  facet_wrap(~ group) +
  labs(x = "FH tree cover (%)",
       y = "AH minus FH cover (percentage points)",
       title = "Cover hidden by first-hit, against tree cover",
       subtitle = "If trees hide understorey, the difference should rise with tree cover") 

ggsave(file.path(fig_dir, "diff_vs_tree_cover.png"), p_diff,
       width = 9, height = 4, dpi = 300)

p_diff_fht <- paired2 |>
  filter(group != "tree") |>
  ggplot(aes(fh_tree, diff_fht)) +
  base() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  smooth()+
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  facet_wrap(~ group) +
  labs(x = "FH tree cover (%)",
       y = "FHT minus FH cover (percentage points)",
       title = "None-tree (FHT) first hit cover, against tree cover",
       subtitle = "If trees hide understorey, the difference should rise with tree cover") 

ggsave(file.path(fig_dir, "diff_fht_vs_tree_cover.png"), p_diff_fht,
       width = 9, height = 4, dpi = 300)

# ---- 3. herbaceous difference against shrub cover -------------------------
# Shrubs hide herbaceous cover by the same mechanism as trees, one layer
# down, and are far more widespread in these plots.
p_diff_shrub <- paired2 |>
  filter(group == "herbaceous") |>
  ggplot(aes(fh_shrub, diff)) +
  base() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  smooth()+
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  labs(x = "FH shrub cover (%)",
       y = "AH minus FH herbaceous cover (percentage points)",
       title = "Herbaceous cover hidden by first-hit, against shrub cover")

ggsave(file.path(fig_dir, "herb_diff_vs_shrub_cover.png"), p_diff_shrub,
       width = 5.5, height = 4, dpi = 300)

# ---- 4. herbaceous difference against combined overstorey ----------------
# Trees and shrubs both sit above the herbaceous layer, so the quantity that
# should predict hidden herbaceous cover is whatever is above it.
p_diff_over <- paired2 |>
  filter(group == "herbaceous") |>
  mutate(fh_over = fh_tree + fh_shrub) |>
  ggplot(aes(fh_over, diff)) +
  base() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  smooth()+
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  labs(x = "first-hit tree + shrub cover",
       y = "Any-hit minus first-hit herbaceous cover (percentage points)",
       title = "Herbaceous cover hidden by first-hit, against overstory cover") 

ggsave(file.path(fig_dir, "herb_diff_vs_overstory.png"), p_diff_over,
       width = 5.5, height = 4, dpi = 300)

# ---- 5. share of the FH-to-AH shortfall recovered by FHT ------------------
# If this is near 1, FHT is effectively AH and the extra definition buys
# nothing. If near 0, trees are not what hides understorey and FHT is not
# worth carrying.
p_recovered <- paired2 |>
  filter(group != "tree", diff > 0) |>
  mutate(recovered = diff_fht / diff) |>
  ggplot(aes(fh_tree, recovered)) +
  geom_hex(bins = 100) +
  smooth() +
  scale_fill_viridis_c(trans = "log10", name = "plot-years") +
  facet_wrap(~ group) +
  labs(x = "FH tree cover (%)", y = "(FHT - FH) / (AH - FH)",
       title = "Share of the first-hit shortfall recovered by FHT",
       caption = caption)

ggsave(file.path(fig_dir, "fht_share_recovered.png"), p_recovered,
       width = 9, height = 4, dpi = 300)

# ---- summary numbers -----------------------------------------------------
paired |>
  summarise(median_fh = round(median(fh), 2),
            median_fht = round(median(fht), 2),
            median_ah = round(median(ah), 2),
            pct_fht_gt_fh = round(100 * mean(diff_fht > 0), 1),
            pct_ah_gt_fh = round(100 * mean(diff > 0), 1),
            .by = group) |>
  as.data.frame() |> print()

