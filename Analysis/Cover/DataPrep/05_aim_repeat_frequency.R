# 05_aim_repeat_frequency.R
#
# Repeat frequency of AIM sampling, used as the temporal density target for
# the RAP augmentation.
#
# Counted two ways:
#   per plot   - (PlotID, daymet cell); PlotID alone is reused across
#                field offices
#   per pixel  - the Daymet 1 km snap grid, which is the unit augmentation
#                operates on. Distinct plots in one pixel count as repeats of
#                that pixel, so this runs higher than the per-plot figure.
#
# Pixel-years, not rows: two plots sampled in the same pixel in the same year
# are one pixel-year, since they are averaged. RAP can only add pixel-years,
# so n_years per pixel is the quantity the target must match.
#
# Input:  cover_by_plot_year.csv, daymet_conus_snap_1000m.tif
# Output: printed summary; aim_repeat_target.csv
#
# given very low temporal reameasurement, likely won't worry about
# getting getting repeats from RAP data (i.e. output of this script not used)
#
# September, 2026


# dependencies ------------------------------------------------------------


source("Functions/init.R")
source_functions()

# read in data ------------------------------------------------------------


ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")

ldc <- read_csv(file.path(ldc_processed_dir, "cover_by_plot_year.csv"),
                show_col_types = FALSE)

aim <- ldc |> filter(ProjectKey == "BLM_AIM")

message("AIM plot-years: ", nrow(aim),
        " | years ", min(aim$year), "-", max(aim$year))

# ---- assign the 1 km snap grid cell --------------------------------------
snap <- read_mask()

aim_vect <- terra::vect(aim, geom = c("Longitude_NAD83", 'Latitude_NAD83'), 
                        crs = "EPSG:4269") |> 
  terra::project(terra::crs(snap))
aim$cell <- terra::extract(snap, aim_vect, cells = TRUE)$cell

n_off <- sum(is.na(aim$cell))

aim <- aim |> filter(!is.na(cell))

# ---- per plot ------------------------------------------------------------
repeats_plot <- aim |>
  summarise(n_years = n_distinct(year),
            span    = max(year) - min(year),
            .by = c(PlotID, cell))

# ---- per pixel -----------------------------------------------------------
repeats_px <- aim |>
  summarise(n_years = n_distinct(year),
            n_plots = n_distinct(paste(ProjectKey, PlotID)),
            n_rows  = n(),
            span    = max(year) - min(year),
            .by = cell)

# ---- summary -------------------------------------------------------------
message("\nyears per plot:")
count(repeats_plot, n_years) |>
  mutate(pct = round(100 * n / sum(n), 1)) |> as.data.frame() |> print()

message("\nyears per 1 km pixel:")
count(repeats_px, n_years) |>
  mutate(pct = round(100 * n / sum(n), 1)) |> as.data.frame() |> print()

message("\nplots per 1 km pixel:")
count(repeats_px, n_plots) |>
  mutate(pct = round(100 * n / sum(n), 1)) |> as.data.frame() |> print()

target <- bind_rows(
  plot  = select(repeats_plot, n_years),
  pixel = select(repeats_px,   n_years),
  .id = "unit"
) |>
  summarise(n_units      = n(),
            mean_years   = mean(n_years),
            median_years = median(n_years),
            q75_years    = quantile(n_years, 0.75),
            pct_repeat   = 100 * mean(n_years > 1),
            .by = unit) |>
  mutate(across(where(is.numeric), \(z) round(z, 3)))

message("\ntemporal density target:")
as.data.frame(target) |> print()

message("\nAIM plot-years by year:")
count(aim, year) |> arrange(year) |> as.data.frame() |> print()

write_csv(target, file.path(ldc_processed_dir, "aim_repeat_target.csv"))