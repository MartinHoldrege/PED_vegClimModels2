# Part 1:
#
# Observed sampling density against a target density, by EPA Level III
# ecoregion, for each cover group.
#
# Observations are the unmasked field pixel-years. A pixel-year from FIA
# alone is masked only by fire: FIA already excludes developed and
# agricultural conditions, so the LCMAP mask would drop plots twice. Any
# pixel-year involving another source is masked by both.
#
# Available area is the count of snap cells passing both the LCMAP mask and
# the decade-mean fire mask. Note the two are not the same fire criterion as
# the observation flags: observations use the per-year mask (that year's
# 20-year window), area uses the mean burned fraction over 2010-2023. Area is
# a property of the ecoregion, not of a year.
#
# Target density is one sampled pixel per ~2,500 ha, the FIA base grid.
#
# Inputs:
#   field_cover_by_pixel_year.csv    - 05_combine_field_cover.R
#   EPA_L3_ecoregion_daymet_1000m.tif, EPA_L3_ecoregion_lookup.csv
#                                    - 01_rasterize_ecoregions.R
#   LCMAP_fracKeep_gte90_1000m.tif   - 03_export_masks.js, via
#                                      04_download_gee_output.R
#   MTBS_fracUnburnedMean_gte90_20yr_2010-2023_1000m.tif
#                                    - 03_export_masks.js, via
#                                      04_download_gee_output.R
#   daymet_conus_snap_1000m.tif      - 00_create_snap_raster.R (read_mask())
#
# Output:
#   ecoregion_sampling_gap.csv
# 
# Part 2: sample rap pixels to augment the field data,
# in ecoregions where field data is sparse (see full description
# at begin of part 2 below)
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

vc <- opt$vc # 'c01' # cover output version (code not yet setup to do anything, 
# different depending on different versions)

# Part 1 ------------------------------------------------------------------

# params ------------------------------------------------------------------

in_dir   <- file.path(paths$large, "Data_processed/cover_combined")
mask_dir <- file.path(paths$large, "Data_processed/masks")
out_dir  <- in_dir


target_density_km <- 24.28  # one sampled pixel per this much available area
# based on FIA based grid which is 1 sample per ~6,000 acres

# cover groups and the column each is observed in
groups <- c(tree = "cov_tree", shrub = "cov_shrub",
            herbaceous = "cov_herbaceous", bare_ground = "cov_bare_ground")


# read in data ------------------------------------------------------------

snap <- read_mask()

field <- read_csv(file.path(in_dir, "field_cover_by_pixel_year.csv"),
                  show_col_types = FALSE)

eco        <- load_ecoregion_raster("L3")
eco_lookup <- load_ecoregion_lookup("L3")


lcmap_mask <- terra::rast(file.path(mask_dir, "LCMAP_fracKeep_gte90_1000m.tif")) |>
  # GEE exports come back one cell larger in each direction than the snap grid
  align_raster(snap)

fire_mask  <- terra::rast(file.path(
  mask_dir, "MTBS_fracUnburnedMean_gte90_20yr_2011-2023_1000m.tif")) |>
  align_raster(snap)

# available area by ecoregion ---------------------------------------------

names(lcmap_mask) <- "lcmap_keep"
names(fire_mask)  <- "fire_keep"

cells <- c(eco, lcmap_mask, fire_mask, cellSize(snap, unit = 'km'),
           snap) |>
  terra::as.data.frame(cells = TRUE, na.rm = FALSE) |>
  as_tibble() 

# checking the snap cell_id is the same as just the cellnumbers
stopifnot(!any(cells$cell != cells$cell_id, na.rm = TRUE),
         !any(!is.na(cells$cell_id) & is.na(cells$cell)))

cells <- cells |> 
  filter(!is.na(cell_id)) |> 
  select(-cell_id)

area_eco <- cells |>
  summarise(
    area_total     = sum(area),
    area_lcmap     = sum(area[lcmap_keep == 1], na.rm = TRUE),
    area_fire      = sum(area[fire_keep == 1], na.rm = TRUE),
    area_available = sum(area[lcmap_keep == 1 & fire_keep == 1], na.rm = TRUE),
    .by = region
  ) |>
  mutate(
    frac_available    = area_available / area_total,
    n_target          = area_available / target_density_km
  ) |> 
  filter(!is.na(region))


# observed pixels by ecoregion and group ----------------------------------

field_eco <- field |>
  left_join(select(cells, cell, region, lcmap_keep, fire_keep), by = "cell") |>
  filter(!is.na(region))

n_no_region <- nrow(field) - nrow(field_eco)
if (n_no_region > 0) {
  # when run this flag didn't get thrown
  warning("pixel-years with no ecoregion, dropped: ", n_no_region)
}

# for the final dataset FIA-only pixel-years are masked by fire alone; everything else by both
# but here deciding to only count the fia pixels after masking by both
# so the actually sampled and target area refers to the same potential 
# pixel area (only makes about a ~2% difference)
field_use <- field_eco |>
  mutate(keep = !masked_by_fire & !masked_by_lcmap & 
           # for this purpose also filtering by the mean unburned 
           # mask which is used to estimate 'available' area
           fire_keep == 1) |>
  filter(keep)

message("pixel-years: ", nrow(field_eco), 
        "\nunmasked pixel-years: ", 
        nrow(distinct(field_use[c('cell', 'year')])),
        '\nunmasked pixels: ', 
        length(unique(field_use$cell)))

# repeats
repeats <- field_use |> 
  summarise(n = n(), .by = c(sources, cell)) |> 
  summarise(perc_gt1 = mean(n >1)*100,
            perc_2 = mean(n == 2)*100,
            .by = 'sources') 

# proportion of fia cells that got remeasured--this 
# is the target repeat frequency when sampling rap cells
target_repeats <- repeats$perc_gt1[repeats$sources == 'fia'] |> 
  round(1)

#' Observed pixel and pixel-year counts for one cover group
#'
#' @param col Column holding that group's cover.
#' @param label Group name for the output.
observed_group <- function(col, label, data = field_use) {
  data |>
    filter(!is.na(.data[[col]])) |>
    summarise(n_pixels      = n_distinct(cell),
              n_pixel_years = n(),
              .by = region) |>
    mutate(group = label)
}

observed <- imap(groups, \(col, label) observed_group(col, label)) |>
  list_rbind()


# gap ---------------------------------------------------------------------
# how many missing observations based on the target observation density

gap <- expand_grid(region = area_eco$region, group = names(groups)) |>
  left_join(area_eco, by = "region") |>
  left_join(observed, by = c("region", "group")) |>
  mutate(
    n_pixels       = replace_na(n_pixels, 0),
    n_pixel_years  = replace_na(n_pixel_years, 0),
    n_gap          = pmax(n_target - n_pixels, 0),
    pct_of_target  = 100 * n_pixels / n_target
  ) |>
  left_join(eco_lookup, by = c("region" = "eco_id")) |>
  select(region, eco_code, eco_name, group, frac_available,
         area_total, area_available,
         n_pixels, n_pixel_years, n_target, n_gap, pct_of_target) |>
  arrange(group, desc(n_gap))

# * check for fia only ----------------------------------------------------
# understanding what the observational density in the fia data
# --it's above the target density in a few ecoregions, because some regions
# sample more densely
field_use_fia <- field_use |> 
  filter(str_detect(sources, 'fia'))
observed_fia <- imap(groups, \(col, label) observed_group(col, label,data = field_use_fia)) |>
  list_rbind()

gap_fia <- expand_grid(region = area_eco$region, group = names(groups)) |>
  left_join(area_eco, by = "region") |>
  left_join(observed_fia, by = c("region", "group")) |>
  mutate(
    n_pixels       = replace_na(n_pixels, 0),
    n_pixel_years  = replace_na(n_pixel_years, 0),
    n_gap          = pmax(n_target - n_pixels, 0),
    n_target_theory = area_total/target_density_km,
    pct_of_target  = 100 * n_pixels / n_target,
    pct_of_target_theory  = 100 * n_pixels / n_target_theory
  ) |>
  left_join(eco_lookup, by = c("region" = "eco_id")) |>
  select(region, eco_code, eco_name, group, frac_available,
         area_total, area_available,
         n_pixels, n_pixel_years, n_target, n_gap, pct_of_target, 
         matches('theory')) |>
  arrange(group, desc(n_gap))
  
gap_fia |> 
  filter(group == 'tree', pct_of_target > 100) |> 
  pull(eco_name)

# checks ------------------------------------------------------------------

message("\necoregions: ", n_distinct(gap$region),
        " | with a lookup entry: ", sum(!is.na(unique(gap$eco_code))))

message("\navailable fraction of ecoregion area:")
area_eco |>
  summarise(min = round(min(frac_available), 3),
            median = round(median(frac_available), 3),
            max = round(max(frac_available), 3)) |>
  as.data.frame() |> print()

message("\nCONUS totals by group:")
gap |>
  summarise(n_pixels = sum(n_pixels),
            n_gap    = round(sum(n_gap)),
            .by = group) |>
  as.data.frame() |> print()

message("\necoregions already at target, by group:")
gap |>
  summarise(n_eco = n(),
            n_at_target = sum(pct_of_target >= 100),
            .by = group) |>
  as.data.frame() |> print()

# 
message("\nlargest gaps (tree):")
gap |> filter(group == "tree") |>
  select(eco_name, frac_available, n_pixels, n_target, n_gap) |>
  head(15) |> as.data.frame() |> print()

stopifnot(all(!is.na(gap$n_target)))

write_csv(gap, file.path(out_dir, "ecoregion_sampling_gap.csv"))


# Part 2: RAP augmentation ---------------------------------------------------
# Fill the per-ecoregion gap with RAP pixels, subject to:
#   - only cells in the available area (LCMAP + decade-mean fire masks), which is
#     the same area the target was computed from
#   - only cells on the thinned (every 5th) RAP grid, i.e. cells with RAP
#     data. RAP is already masked by LCMAP and the per-year fire mask, so a
#     non-NA value means that cell-year passed both.
#   - never a cell that appears anywhere in the observational data
#
# Groups are nested: one ordered random sample of size max(n_gap) is drawn
# per ecoregion, and group g takes the first n_gap[g] of it. So the pixels
# filling a small gap are a subset of those filling a larger one, and a
# single RAP row can carry several groups.
#
# Each selected pixel gets one random year, and a `target_repeats` share of
# them get a second. Years are drawn only from years where that cell has RAP
# data. The gap is in pixels, so a second year is extra rather than counting
# against it.
#
# RAP understorey is not usable under tree canopy, so cov_shrub and
# cov_herbaceous are set to NA in any RAP pixel-year with tree cover at or
# above `tree_cutoff`. That means shrub and herbaceous gaps will be
# under-filled in forested ecoregions - reported below.

set.seed(1234)

rap_dir     <- file.path(paths$large, "Data_processed/CoverData/rap")
rap_years   <- 2011:2023
tree_cutoff <- 10   # percent; RAP understorey invalid at or above this

# RAP file per group; band names are year_YYYY
# files created in 03_rap_sample.js. these have been already
# masked (i.e. 90% unburned in prior 20 years, and not ag, developed or water)
rap_files <- c(
  tree        = "RAP_v3_cover-tree_2011-2023_thin5_1000m.tif",
  shrub       = "RAP_v3_cover-shrub_2011-2023_thin5_1000m.tif",
  herbaceous  = "RAP_v3_cover-herbaceous_2011-2023_thin5_1000m.tif",
  bare_ground = "RAP_v3_cover-bare_ground_2011-2023_thin5_1000m.tif"
)
stopifnot(all(file.exists(file.path(rap_dir, rap_files))))


# * observational data, with the site-year filters only -------------------
# This is the observational half of the final dataset: the decade-mean fire
# restriction used for the gap calculation is NOT applied here, since a plot
# measured before a fire is still a valid observation of that year.
field_final <- field_eco |>
  mutate(keep = if_else(sources == "fia",
                        # fia data was pre-filtered by ground condition codes (e.g. avoiding developed)
                        !masked_by_fire, 
                        !masked_by_fire & !masked_by_lcmap)) |>
  filter(keep) |>
  select(-keep, -lcmap_keep, -fire_keep)

message("\nobservational pixel-years in the final dataset: ", nrow(field_final))


# * read RAP, wide by year ------------------------------------------------

#' Read one RAP cover raster as a cell-by-year table
#'
#' Rows are dropped only where every year is NA, so the result is the thinned
#' grid. Values are percent cover.
#'
#' @param file File name within rap_dir.
read_rap_wide <- function(file) {
  terra::rast(file.path(rap_dir, file)) |>
    align_raster(snap) |>
    terra::as.data.frame(cells = TRUE, na.rm = NA) |>
    as_tibble() |>
    rename(cell = cell)
}

rap_wide <- map(rap_files, read_rap_wide)

stopifnot(all(map_lgl(rap_wide, \(d) identical(sort(d$cell),
                                               sort(rap_wide$tree$cell)))))

message("RAP cells on the thinned grid: ", nrow(rap_wide$tree))


# * candidate cells -------------------------------------------------------

observed_cells <- unique(field_final$cell)

candidates <- cells |>
  filter(!is.na(region)) |>
  filter(cell %in% rap_wide$tree$cell) |>
  filter(!cell %in% observed_cells) |>
  select(cell, region)

message("candidate RAP cells: ", nrow(candidates),
        " | in ", n_distinct(candidates$region), " ecoregion(s)")


# * how many cells to draw per ecoregion ----------------------------------
# n_gap is fractional (area / target density); rounded up so a partial gap
# still gets a pixel.
gap_int <- gap |>
  mutate(n_gap_int = pmax(ceiling(n_gap), 0)) |>
  select(region, group, n_gap_int)

eco_draw <- gap_int |>
  summarise(n_draw = max(n_gap_int), .by = region) |>
  filter(n_draw > 0)

pool <- candidates |> summarise(n_pool = n(), .by = region)

eco_draw <- eco_draw |>
  left_join(pool, by = "region") |>
  mutate(n_pool = replace_na(n_pool, 0),
         n_draw_actual = pmin(n_draw, n_pool),
         shortfall = n_draw - n_draw_actual)

message("\necoregions needing RAP: ", nrow(eco_draw),
        " | with too few candidate cells: ", sum(eco_draw$shortfall > 0))

if (any(eco_draw$shortfall > 0)) {
  eco_draw |> filter(shortfall > 0) |>
    left_join(select(eco_lookup, eco_id, eco_name), by = c("region" = "eco_id")) |>
    arrange(desc(shortfall)) |>
    select(eco_name, n_draw, n_pool, shortfall) |>
    head(15) |> as.data.frame() |> print()
}


# * draw the nested sample ------------------------------------------------
# `rank` is the position in the ecoregion's ordered draw; group g keeps the
# pixels with rank <= its own gap, which gives the nesting.
sampled_cells <- candidates |>
  # inner join b/ a couple ecoregions don't require any draws
  inner_join(select(eco_draw, region, n_draw_actual), by = "region") |>
  filter(n_draw_actual > 0) |>
  slice_sample(prop = 1) |> # random reordering of rows
  mutate(rank = row_number(), .by = region) |>
  filter(rank <= n_draw_actual) |>
  select(cell, region, rank)

message("RAP cells drawn: ", nrow(sampled_cells))


# * assign years ----------------------------------------------------------
# Years are drawn only from years where the cell has RAP data, so a cell
# masked by fire in some years can still be used in others.
years_available <- rap_wide$tree |>
  filter(cell %in% sampled_cells$cell) |>
  pivot_longer(-cell, names_to = "year", values_to = "cover") |>
  filter(!is.na(cover)) |>
  mutate(year = as.integer(str_remove(year, "^year_"))) |>
  select(cell, year)

#' Draw one or two years for a cell
#'
#' @param yrs Years with RAP data for that cell.
#' @param second TRUE if a second year should be drawn.
draw_years <- function(yrs, second) {
  n <- if (second && length(yrs) > 1) 2 else 1
  # see ?sample
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  resample(x = yrs, size = n)
}

sampled_year <- years_available |>
  summarise(yrs = list(year), .by = cell) |>
  mutate(second = runif(n()) < target_repeats / 100,
         year = map2(yrs, second, draw_years)) |>
  select(cell, year) |>
  unnest(year)

message("RAP pixel-years: ", nrow(sampled_year),
        " | % pixels with two years: ",
        mean(count(sampled_year, cell)$n > 1)*100,
        " (target ", target_repeats, "%)")

# RAP pixel-years: 144447 | % pixels with two years: 53.3603006752453 (target 54.1%)
# * attach cover values ---------------------------------------------------

rap_long <- imap(rap_wide, \(d, g) {
  d |>
    filter(cell %in% sampled_cells$cell) |>
    pivot_longer(-cell, names_to = "year", values_to = "cover") |>
    mutate(year = as.integer(str_remove(year, "^year_")),
           group = g) |>
    filter(!is.na(cover))
}) |>
  list_rbind()

rap_px <- sampled_year |>
  left_join(sampled_cells, by = "cell") |>
  left_join(rap_long, by = c("cell", "year"),
            relationship = "many-to-many") |>
  pivot_wider(names_from = group, values_from = cover, names_prefix = "cov_")

# a group only keeps the pixels its own gap reached
rap_px <- rap_px |>
  left_join(pivot_wider(gap_int, names_from = group, values_from = n_gap_int,
                        names_prefix = "gap_"),
            by = "region") |>
  mutate(
    cov_tree        = if_else(rank <= gap_tree,        cov_tree,        NA_real_),
    cov_shrub       = if_else(rank <= gap_shrub,       cov_shrub,       NA_real_),
    cov_herbaceous  = if_else(rank <= gap_herbaceous,  cov_herbaceous,  NA_real_),
    cov_bare_ground = if_else(rank <= gap_bare_ground, cov_bare_ground, NA_real_)
  )

# RAP understorey is not usable under tree canopy
n_under_tree <- sum(rap_px$cov_tree >= tree_cutoff, na.rm = TRUE)

rap_px <- rap_px |>
  mutate(under_tree = !is.na(cov_tree) & cov_tree >= tree_cutoff,
         cov_shrub      = if_else(under_tree, NA_real_, cov_shrub),
         cov_herbaceous = if_else(under_tree, NA_real_, cov_herbaceous))

message("\nRAP pixel-years with tree cover >= ", tree_cutoff,
        "%, understorey dropped: ", n_under_tree,
        " (", round(100 * n_under_tree / nrow(rap_px), 1), "%)")


# * combine with the observational data -----------------------------------

xy <- terra::xyFromCell(snap, rap_px$cell)

rap_out <- rap_px |>
  transmute(cell, year, region,
            x = xy[, "x"], y = xy[, "y"],
            sources = "rap", n_plots = NA_integer_,
            cov_tree, cov_shrub, cov_herbaceous, cov_bare_ground)

cover_final <- field_final |> 
  select(-masked_by_fire, -masked_by_lcmap) |> 
  bind_rows(rap_out) |>
  left_join(eco_lookup, by = c('region' = 'eco_id')) |> 
  select(-region) |> 
  arrange(eco_code, cell, year) 

test <- cover_final |> 
  select(cell, year, x, y, eco_code, eco_name, sources) |> 
  is.na() |> 
  sum()

stopifnot(test == 0)

# * checks ----------------------------------------------------------------

message("\nfinal dataset: ", nrow(cover_final), " pixel-years")
count(cover_final, sources, sort = TRUE) |>
  mutate(pct = round(100 * n / sum(n), 1)) |> as.data.frame() |> print()

message("\nnon-missing cover by group and source:")
cover_final |>
  summarise(across(all_of(unname(groups)), \(z) sum(!is.na(z))), .by = sources) |>
  as.data.frame() |> print()

# RAP must never share a pixel-year with field data
overlap <- inner_join(select(field_final, cell, year),
                      select(rap_out, cell, year), by = c("cell", "year"))
stopifnot(nrow(overlap) == 0)

# and never a pixel that has field data at all
stopifnot(!any(rap_out$cell %in% field_final$cell))

# nesting: a group's RAP pixels must be a subset of any larger group's
rap_cells_by_group <- map(unname(groups), \(col) {
  rap_out |> filter(!is.na(.data[[col]])) |> distinct(region, cell)
}) |> set_names(names(groups))

# getting the cell id's of group with the most cells in a region
rap_cells_largest <- bind_rows(rap_cells_by_group, .id = 'group') |> 
  mutate(n = n(), .by = c(region, group)) |> 
  filter(group == group[n == max(n)][1]) |> 
  select(region, cell)

# check
walk(rap_cells_by_group, function(df) {
  stopifnot(all(df$cell %in% rap_cells_largest$cell))
})

# covers stay in range
stopifnot(all(map_lgl(unname(groups), \(cc) {
  z <- cover_final[[cc]]; all(is.na(z) | (z >= 0 & z <= 100))
})))

# no duplicate pixel-years
stopifnot(!any(duplicated(cover_final[c("cell", "year")])))

message("\nresulting density against target, by group:")
density_achieved <- cover_final |>
  filter(sources != "rap" | TRUE) |>
  summarise(across(all_of(unname(groups)),
                   \(z) n_distinct(cell[!is.na(z)])), .by = eco_name) |>
  pivot_longer(-eco_name, names_to = "col", values_to = "n_pixels_final") |>
  mutate(group = names(groups)[match(col, unname(groups))]) |>
  left_join(select(gap, eco_name, eco_code, group, n_target), by = c("eco_name", "group")) 
density_achieved |> 
  summarize(across(c(n_target, n_pixels_final), .fns = sum), .by = group) |> 
  print()

# saving ------------------------------------------------------------------

write_csv(cover_final, file.path(out_dir, paste0("cover_by_pixel_year_all-sources_", 
                                                 vc, ".csv")))

write_csv(density_achieved, 
          paste0('diagnostics_pixel-n_achieved-vs-target_', vc, '.csv'))
