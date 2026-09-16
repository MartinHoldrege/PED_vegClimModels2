# 05_combine_field_cover.R
#
# Combine the LDC, LFRDB and FIA plot-year cover tables into one pixel-year
# table on the Daymet 1 km snap grid.
#
# Cover definitions differ by source:
#   LDC    fht_*  first pin hit, seeing through trees
#   LFRDB  cov_*  LANDFIRE overlap-adjusted, else the crew's lifeform call
#   FIA    cov_*  P2VEG aerial cover (tally + non-tally for trees)
# All are percent cover and are combined directly. `sources` records which
# contributed to each pixel-year so the difference can be modelled.
#
# Component covers are averaged across plot-years within a pixel, and the
# fractions computed once from those averages. Averaging the fractions
# instead would give a plot with trace herbaceous cover the same weight as
# one where the layer dominates.
#
# Herbaceous split: frac_forb + frac_c3 + frac_c4 = 1, so modelled herbaceous
# cover can be disaggregated. FIA has no pathway split; where it is the sole
# source in a pixel-year the denominator falls back to forb plus unsplit
# graminoid cover and the C3/C4 fractions are NA.
#
# Tree split: frac_needle + frac_broad = 1. FIA's is a basal-area fraction
# over tally species only, applied to a tree cover that includes non-tally
# species; it is converted to a cover-scale component here so it can be
# averaged with the other sources.
#
# Bare ground is LDC only.
#
# Inputs:
#   cover_by_plot_year.csv             - 04_ldc_process.R
#   lfrdb_cover_by_plot.csv            - 02_lfrdb_process.R
#   fia_cover_by_plot_year<suffix>.csv - 02_FIA_combine.R
#   daymet_conus_snap_1000m.tif        - 00_create_snap_raster.R (read_mask())
#   LCMAP_fracKeep_gte90_1000m.tif     - 03_export_masks.js, downloaded in
#                                        04_download_gee_output.R
#   MTBS_fracUnburned_gte90_20yr_2000-2024_1000m.tif
#                                      - 03_export_masks.js, downloaded in
#                                        04_download_gee_output.R
#
# Output:
#   field_cover_by_pixel_year.csv
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
# using some functions in general.R and spatial.R

# params ------------------------------------------------------------------

ldc_dir   <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")
lfrdb_dir <- file.path(paths$large, "Data_processed/LANDFIRE_LFRDB")
fia_dir   <- file.path(paths$large, "Data_processed/FIA")
out_dir   <- file.path(paths$large, "Data_processed/cover_combined")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fia_suffix <- "_fia1"

snap <- read_mask()

cov_cols  <- c("cov_tree", "cov_shrub", "cov_herbaceous", "cov_bare_ground")
comp_cols <- c("c_forb", "c_c3", "c_c4", "c_gram", "c_needle", "c_broad") # component cover
frac_cols <- c("frac_forb", "frac_c3", "frac_c4", "frac_needle", "frac_broad")

# read and standardise each source ----------------------------------------

# ---- LDC -----------------------------------------------------------------
# ah_grass_unknown is left out of the components: grass of unknown pathway
# would otherwise inflate the forb fraction.
ldc <- read_csv(file.path(ldc_dir, "cover_by_plot_year.csv"),
                show_col_types = FALSE) |>
  assign_cell("Longitude_NAD83", "Latitude_NAD83", "EPSG:4269") |>
  drop_offgrid() |>
  transmute(
    source = "ldc",
    cell, year,
    cov_tree        = fht_tree,
    cov_shrub       = fht_shrub,
    cov_herbaceous  = fht_herbaceous,
    cov_bare_ground = fh_bare_ground,
    c_forb   = ah_forb,
    c_c3     = ah_grass_c3,
    c_c4     = ah_grass_c4,
    c_gram   = NA_real_,     # unsplit graminoid cover; FIA only
    c_needle = ah_tree_needle,
    c_broad  = ah_tree_broad
  )

# ---- LFRDB ---------------------------------------------------------------
lfrdb <- read_csv(file.path(lfrdb_dir, "lfrdb_cover_by_plot.csv"),
                  show_col_types = FALSE) |>
  assign_cell("Longitude_wgs84", "Latitude_wgs84", "EPSG:4326") |>
  drop_offgrid() |>
  transmute(
    source = "lfrdb",
    cell, year,
    cov_tree        = cov_tree,
    cov_shrub       = cov_shrub,
    cov_herbaceous  = cov_herbaceous,
    cov_bare_ground = NA_real_,
    c_forb   = sum_forb,
    c_c3     = sum_grass_c3,
    c_c4     = sum_grass_c4,
    c_gram   = NA_real_,
    c_needle = sum_tree_needle,
    c_broad  = sum_tree_broad
  )

# ---- FIA -----------------------------------------------------------------
# Graminoid cover cannot be split by pathway, so it is carried as c_gram and
# used only where no C3/C4 source is present in the pixel-year. The
# needle/broad components are reconstructed from the basal-area fraction
# against total tree cover, putting them on a cover scale.
fia <- read_csv(file.path(fia_dir,
                          paste0("fia_cover_by_plot_year", fia_suffix, ".csv")),
                show_col_types = FALSE) |>
  assign_cell("Longitude_NAD83", "Latitude_NAD83", "EPSG:4269") |>
  drop_offgrid() |>
  transmute(
    source = "fia",
    cell, year,
    cov_tree        = cov_tree,
    cov_shrub       = cov_shrub,
    cov_herbaceous  = cov_herbaceous,
    cov_bare_ground = NA_real_,
    c_forb   = cov_forb,
    c_c3     = NA_real_,
    c_c4     = NA_real_,
    c_gram   = cov_graminoid,
    c_needle = cov_tree * frac_needle,
    c_broad  = cov_tree * frac_broad
  )


# combine to pixel-year ---------------------------------------------------

field_long <- bind_rows(ldc, lfrdb, fia)
stopifnot(all(c("source", "cell", "year", cov_cols, comp_cols) %in%
                names(field_long)))

message("\nplot-years by source:")
count(field_long, source) |> as.data.frame() |> print()

field_px <- field_long |>
  summarise(across(all_of(c(cov_cols, comp_cols)), mean_na),
            n_plots = n(),
            sources = paste(sort(unique(source)), collapse = ","),
            .by = c(cell, year)) |>
  mutate(
    # Where a pathway split exists the herbaceous denominator is
    # forb + C3 + C4; otherwise it falls back to forb plus unsplit graminoid
    # cover. In a mixed pixel-year carrying both, FIA's graminoid cover drops
    # out of the denominator while its forb cover stays in - check the
    # `sources` tally for how often that arises.
    herb_den    = if_else(!is.na(c_c3) | !is.na(c_c4),
                          coalesce(c_forb, 0) + coalesce(c_c3, 0) +
                            coalesce(c_c4, 0),
                          coalesce(c_forb, 0) + coalesce(c_gram, 0)),
    frac_forb   = safe_frac(c_forb, herb_den),
    frac_c3     = safe_frac(c_c3,   herb_den),
    frac_c4     = safe_frac(c_c4,   herb_den),
    frac_needle = safe_frac(c_needle, c_needle + c_broad),
    frac_broad  = safe_frac(c_broad,  c_needle + c_broad)
  ) |>
  select(-herb_den) |> 
  filter(year >= 2000 & year < 2025)

# cell centroid, for joining to other gridded data and for mapping
xy <- terra::xyFromCell(snap, field_px$cell)
field_px$x <- xy[, "x"]
field_px$y <- xy[, "y"]


# mask flags --------------------------------------------------------------
# Flags only, no filtering: whether to apply the LCMAP mask to FIA is decided
# downstream, and FIA already excludes developed and agricultural conditions.
# TRUE means the pixel would be masked out.
#
# The fire mask is per-year (less than 10% of the cell burned in the
# preceding 20 years), so each pixel-year is checked against its own band.
# Pixel-years outside the mask's year range get NA.

mask_dir <- file.path(paths$large, "Data_processed/masks")

lcmap_mask <- terra::rast(file.path(mask_dir, "LCMAP_fracKeep_gte90_1000m.tif"))
fire_mask  <- terra::rast(file.path(mask_dir,
                                    "MTBS_fracUnburned_gte90_20yr_2000-2024_1000m.tif"))

lcmap_mask <- align_raster(lcmap_mask, snap)
fire_mask <- align_raster(fire_mask, snap)

# LCMAP: one layer, so a straight cell lookup
field_px$masked_by_lcmap <- !as.logical(lcmap_mask[field_px$cell][[1]])

# Fire: pick the band matching each pixel-year's year
fire_years <- as.integer(str_remove(names(fire_mask), "^year_"))

fire_vals <- terra::extract(fire_mask, field_px$cell)
band_idx  <- match(field_px$year, fire_years)
stopifnot(all(!is.na(band_idx))) # suggests more years of data provided than years fire masks made for
field_px$masked_by_fire <- if_else(
  is.na(band_idx), NA,
  !as.logical(fire_vals[cbind(seq_len(nrow(fire_vals)), band_idx)])
)

message("\nmask flags:")
field_px |>
  count(masked_by_lcmap, masked_by_fire) |>
  mutate(pct = round(100 * n / sum(n), 2)) |>
  as.data.frame() |> print()

message("\nby source:")
field_px |>
  summarise(n = n(),
            pct_lcmap = round(100 * mean(masked_by_lcmap, na.rm = TRUE), 1),
            pct_fire  = round(100 * mean(masked_by_fire,  na.rm = TRUE), 1),
            .by = sources) |>
  arrange(desc(n)) |> as.data.frame() |> print()

# checks ------------------------------------------------------------------

message("\npixel-years: ", nrow(field_px),
        " | unique pixels: ", n_distinct(field_px$cell),
        " | from plot-years: ", nrow(field_long))

message("\nsource combinations per pixel-year:")
count(field_px, sources, sort = TRUE) |>
  mutate(pct = round(100 * n / sum(n), 2)) |> as.data.frame() |> print()

message("\nnon-missing counts by column:")
field_px |> summarise(across(all_of(c(cov_cols, comp_cols, frac_cols)),
                             \(z) sum(!is.na(z)))) |> glimpse()

message("\nmedian values:")
field_px |> summarise(across(all_of(c(cov_cols, frac_cols)),
                             \(z) round(median(z, na.rm = TRUE), 3))) |> glimpse()

# covers are percentages, fractions are 0-1
stopifnot(all(map_lgl(cov_cols, \(cc) {
  z <- field_px[[cc]]; all(is.na(z) | (z >= 0 & z <= 100))
})))
stopifnot(all(map_lgl(frac_cols, \(cc) {
  z <- field_px[[cc]]; all(is.na(z) | (z >= 0 & z <= 1))
})))

# fraction sets sum to one wherever any member is present
s_herb <- rowSums(as.matrix(field_px[c("frac_forb", "frac_c3", "frac_c4")]),
                  na.rm = FALSE)
s_herb <- s_herb[!is.na(s_herb)]
s_tree <- rowSums(as.matrix(field_px[c("frac_needle", "frac_broad")]),
                  na.rm = FALSE)
s_tree <- s_tree[!is.na(s_tree)]
stopifnot(all(s_herb == 0 | abs(s_herb - 1) < 1e-6),
          all(s_tree == 0 | abs(s_tree - 1) < 1e-6))

# duplicate pixel-years would break the modelling unit
stopifnot(!any(duplicated(field_px[c("cell", "year")])))

write_csv(field_px, file.path(out_dir, "field_cover_by_pixel_year.csv"))