# 03_ldc_process.R
#
# Compile LDC line-point intercept batches into plot-visit level cover.
#
# Absolute covers use the FIRST (topmost) hit at each pin, so classes are
# mutually exclusive and sum to at most 100%.
#   tree, shrub, herbaceous, bare_ground, litter, rock, ...
#
# Component covers use ANY hit anywhere in the pin column, deduplicated per
# pin, so each is independently bounded at 100% but they do NOT sum to the
# corresponding first-hit total (a pin can hold a C3 grass and a forb).
#   tree_needle, tree_broad, grass_c3, grass_c4, forb
#
# Output is one row per PrimaryKey (plot-visit). Aggregation to plot-year
# happens downstream.
#
# Inputs:
#   code_class.rds  - (SpeciesKey, code) -> growth-habit class, from 02_
#   code_traits.rds - code -> leaf_type, pathway; required for the any-hit
#                     components (see contract below)
#
# Outputs:
#   cover_by_visit.csv     - one row per plot visit (PrimaryKey)
#   cover_by_plot_year.csv - one row per plot-year, with x/y coordinates
#
# September, 2026

source("Functions/init.R")

rerun <- FALSE # rerun chunks if output file already created
ldc_dir <- file.path(paths$large, "Data_raw/LandscapeDataCommonsDat")
lpi_dir <- file.path(ldc_dir, "lpi_batches")
ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")

# Must match 02_ldc_build_code_lookup.R exactly.
layer_order  <- c("TopCanopy", paste0("Lower", 1:10), "SoilSurface")
no_hit_codes <- c("none", "")

hdr_all  <- readRDS(file.path(ldc_dir, "header_raw.rds"))

hdr_key  <- hdr_all |> 
  as_tibble() |> 
  distinct(PrimaryKey, SpeciesKey) 

code_class <- readRDS(file.path(ldc_processed_dir, "code_class.rds")) |>
  as_tibble() |> 
  select(SpeciesKey, code, class)

# ---- species traits for the any-hit components ---------------------------
# Built separately: 02c_ldc_graminoid_pathway.R (C3/C4) and
# 02d_ldc_tree_leaf_type.R (needle/broad). Joined on `code` alone - both are
# resolved from ScientificName, which is code-level rather than list-level.
pathway_path   <- file.path(ldc_processed_dir, "graminoid_pathway.csv")
leaf_type_path <- file.path(ldc_processed_dir, "tree_leaf_type.csv")

stopifnot(file.exists(pathway_path), file.exists(leaf_type_path))

graminoid_pathway <- read_csv(pathway_path) |>
  select(code, pathway) 
tree_leaf_type <- read_csv(leaf_type_path) |>
  select(code, leaf_type) |> 
  distinct()

stopifnot(!any(duplicated(graminoid_pathway$code)),
          !any(duplicated(tree_leaf_type$code)))

code_traits <- full_join(graminoid_pathway, tree_leaf_type, by = "code")
stopifnot(!any(duplicated(code_traits$code)))

# ---- functions -----------------------------------------------------------

#' Long code -> component group lookup for the any-hit covers, first hit covers,
#' and first hit cover, of the understory (i.e. after first hit trees have 
#' been removed)
#'
#' Groups deliberately overlap: a code contributes to its top-level class
#' (tree / shrub / herbaceous) and to its component group (tree_needle,
#' grass_c4, ...). The compile step joins many-to-many on this.
#'
#' @param code_class (SpeciesKey, code, class) growth-habit classification.
#' @param code_traits code -> leaf_type, pathway.
#' @return Tibble with columns SpeciesKey, code, group.
build_code_group <- function(code_class, code_traits) {
  x <- left_join(code_class, code_traits, by = "code")
  
  bind_rows(
    # top-level classes; the any-hit counterparts of fh_tree / fh_shrub /
    # fh_herbaceous, so either definition can be used at model time
    x |> filter(class == "tree")  |> transmute(SpeciesKey, code, group = "tree"),
    x |> filter(class == "shrub") |> transmute(SpeciesKey, code, group = "shrub"),
    x |> filter(class %in% c("graminoid", "forb", "herbaceous", 
                             "nonwoody_unknown")) |>
      transmute(SpeciesKey, code, group = "herbaceous"),
    # components, used only for the proportions - always any-hit
    x |> filter(class == "tree", leaf_type == "needle") |>
      transmute(SpeciesKey, code, group = "tree_needle"),
    x |> filter(class == "tree", leaf_type == "broad") |>
      transmute(SpeciesKey, code, group = "tree_broad"),
    x |> filter(class == "tree", !leaf_type %in% c("needle", "broad")) |>
      transmute(SpeciesKey, code, group = "tree_unknown"),
    x |> filter(class == "graminoid", pathway == "C3") |>
      transmute(SpeciesKey, code, group = "grass_c3"),
    x |> filter(class == "graminoid", pathway == "C4") |>
      transmute(SpeciesKey, code, group = "grass_c4"),
    x |> filter(class == "graminoid", !pathway %in% c("C3", "C4")) |>
      transmute(SpeciesKey, code, group = "grass_unknown"),
    x |> filter(class == "forb") |> transmute(SpeciesKey, code, group = "forb")
  ) |>
    distinct()
}

#' Plot-visit cover from one batch of tall LPI data
#'
#' @param path Path to a cached LPI batch.
#' @param code_class (SpeciesKey, code, class) growth-habit classification.
#' @param code_group Long, overlapping (SpeciesKey, code, group) lookup.
#' @return One row per PrimaryKey: n_pins, fh_* first-hit covers (%),
#'   ah_* any-hit component covers (%).
summarise_lpi <- function(path, code_class, code_group) {
  
  lpi <- readRDS(path) |>
    mutate(code = str_trim(code)) |>
    left_join(hdr_key, by = "PrimaryKey")
  
  # Denominator from ALL rows, before dropping no-hit sentinels, so that a
  # pin recording only "None" still counts as a sampled pin.
  pins <- lpi |>
    distinct(PrimaryKey, LineKey, PointNbr) |>
    count(PrimaryKey, name = "n_pins")
  
  hits <- lpi |>
    filter(!is.na(code), !str_to_lower(code) %in% no_hit_codes) |>
    mutate(rank = match(layer, layer_order))
  
  # A few rows carry a corrupt layer label (e.g. "LowerNA"). They are kept for the
  # any-hit covers, which ignore depth, but cannot be ranked so they are
  # excluded from the first-hit calculation. Warn if this is more than a
  # handful of rows.
  n_bad_layer <- sum(is.na(hits$rank))
  if (n_bad_layer > 0) {
    message(basename(path), ": ", n_bad_layer, " rows with unrankable layer (",
            paste(unique(hits$layer[is.na(hits$rank)]), collapse = ", "), ")")
  }
  
  classify_herbaceous <- function(df) {
    df |> 
      mutate(class = case_when(
        class %in% c("graminoid", "forb", "herbaceous",
                     "nonwoody_unknown")                ~ "herbaceous",
        is.na(class)                                    ~ "unclassified",
        .default = class
    ))
  }
  
  # --- first hit: absolute, mutually exclusive ---------------------------
  fh <- hits |>
    filter(!is.na(rank)) |>
    slice_min(rank, n = 1, with_ties = FALSE,
              by = c(PrimaryKey, LineKey, PointNbr)) |>
    left_join(code_class, by = c("SpeciesKey", "code")) |>
    classify_herbaceous() |> 
    count(PrimaryKey, class, name = "n_hit") |>
    left_join(pins, by = "PrimaryKey") |>
    mutate(cover = 100 * n_hit / n_pins) |>
    select(PrimaryKey, class, cover) |>
    pivot_wider(names_from = class, values_from = cover,
                names_prefix = "fh_", values_fill = 0)
  
  # --- first hit, seeing through trees -----------------------------------
  # A pin contributes its first hit, and if that hit is a tree, also its
  # first non-tree VEGETATION hit. Nothing deeper counts, and surface classes
  # are never picked up on the second pass - bare ground stays strictly
  # first-hit, so fht_bare_ground equals fh_bare_ground.
  #
  # Not mutually exclusive: a pin with a tree over a shrub counts once for
  # tree and once for shrub, so fht_* can sum above 100%.
  veg_classes <- c("tree", "shrub", "herbaceous", "vine", "woody_unknown")
  
  ranked <- hits |>
    filter(!is.na(rank)) |>
    left_join(code_class, by = c("SpeciesKey", "code")) |>
    classify_herbaceous()
  
  fht <- bind_rows(
    # the first hit itself, whatever it is
    ranked |>
      slice_min(rank, n = 1, with_ties = FALSE,
                by = c(PrimaryKey, LineKey, PointNbr)),
    # and, only where that first hit was a tree, the first non-tree
    # vegetation hit beneath it
    ranked |>
      filter(any(class == "tree" & rank == min(rank)),
             class %in% setdiff(veg_classes, "tree"),
             .by = c(PrimaryKey, LineKey, PointNbr)) |>
      slice_min(rank, n = 1, with_ties = FALSE,
                by = c(PrimaryKey, LineKey, PointNbr))
  ) |>
    distinct(PrimaryKey, LineKey, PointNbr, class) |>
    count(PrimaryKey, class, name = "n_hit") |>
    left_join(pins, by = "PrimaryKey") |>
    mutate(cover = 100 * n_hit / n_pins) |>
    select(PrimaryKey, class, cover) |>
    pivot_wider(names_from = class, values_from = cover,
                names_prefix = "fht_", values_fill = 0)
  
  # --- any hit: components, overlapping ----------------------------------
  ah <- hits |>
    inner_join(code_group, by = c("SpeciesKey", "code"),
               relationship = "many-to-many") |>
    distinct(PrimaryKey, LineKey, PointNbr, group) |>
    count(PrimaryKey, group, name = "n_hit") |>
    left_join(pins, by = "PrimaryKey") |>
    mutate(cover = 100 * n_hit / n_pins) |>
    select(PrimaryKey, group, cover) |>
    pivot_wider(names_from = group, values_from = cover,
                names_prefix = "ah_", values_fill = 0)
  
  pins |>
    left_join(fh, by = "PrimaryKey") |>
    left_join(ah, by = "PrimaryKey") |> 
    left_join(fht, by = "PrimaryKey")
}

summarise_lpi_safe <- function(path, code_class, code_group) {
  out <- try(summarise_lpi(path, code_class, code_group), silent = TRUE)
  if (inherits(out, "try-error")) {

    message("SKIPPED ", basename(path), ": ",
            conditionMessage(attr(out, "condition")))
    return(NULL)
  }
  if (nrow(out) == 0) return(NULL)
  out
}


# ---- run -----------------------------------------------------------------
code_group <- build_code_group(code_class, code_traits)

p_cover_visit <- file.path(ldc_processed_dir, "cover_by_visit.csv")

if(!file.exists(p_cover_visit) | rerun) {
  cover_visit <- list.files(lpi_dir, pattern = "^lpi_batch-\\d+\\.rds$",
                            full.names = TRUE) |>
    map(\(f) summarise_lpi_safe(f, code_class, code_group)) |>
    list_rbind() |>
    # a group absent from a batch yields no column there, so fill after binding
    mutate(across(starts_with(c("fh_", "ah_", 'fht_')), \(x) replace_na(x, 0)))
  write_csv(cover_visit, p_cover_visit)
} else {
  cover_visit <- read_csv(p_cover_visit) |> 
    mutate(across(starts_with(c("fh_", "ah_", 'fht_')), \(x) replace_na(x, 0)))
}


# ---- checks --------------------------------------------------------------
# First-hit classes are mutually exclusive, so they must sum to 100%.
fh_total <- cover_visit |>
  transmute(PrimaryKey,
            total = rowSums(across(starts_with("fh_")), na.rm = TRUE))

message("first-hit totals: min ", round(min(fh_total$total), 2),
        ", max ", round(max(fh_total$total), 2))
stopifnot(all(fh_total$total <= 100 + 1e-6))
message("visits below 100%: ", sum(fh_total$total < 100 - 1e-6))

cover_visit |>
  summarise(across(starts_with(c("fh_", "ah_")),
                   \(x) round(median(x), 2))) |>
  glimpse()

message("visits: ", nrow(cover_visit),
        " | median pins: ", median(cover_visit$n_pins))


# ---- aggregate to plot-year ---------------------------------------------
# Plot identity is the rounded coordinate pair, plotid and project
#
# Where a plot-year has more than one visit (~0.5% of plot-years) covers are
# averaged. Most such groups are NWERN research sites where it's likely several distinct
# plots share one recorded coordinate.

coord_digits <- 4   # ~11 m; the plot-identity tolerance

visits <- cover_visit |>
  left_join(
    hdr_all |>
      as_tibble() |>
      select(PrimaryKey, ProjectKey, PlotID, DateVisited,
             Latitude_NAD83, Longitude_NAD83),
    by = "PrimaryKey"
  ) |>
  mutate(DateVisited = as_date(DateVisited),
         year = year(DateVisited),
         Longitude_NAD83 = round(Longitude_NAD83, coord_digits),
         Latitude_NAD83 = round(Latitude_NAD83,  coord_digits))

n_dropped <- sum(is.na(visits$year) | is.na(visits$Latitude_NAD83) | is.na(visits$Longitude_NAD83))
if (n_dropped > 0) {
  message("dropping ", n_dropped, " visits with no date or coordinates")
}
visits <- visits |> filter(!is.na(year), !is.na(Latitude_NAD83), !is.na(Longitude_NAD83))
cover_cols <- str_subset(names(visits), "^(fh|ah|fht)_")

plot_year <- visits |>
  summarise(
    across(all_of(cover_cols), \(z) mean(z, na.rm = TRUE)),
    n_visits   = n(),
    n_pins     = mean(n_pins),
    .by = c(ProjectKey, PlotID, Longitude_NAD83, Latitude_NAD83, year)
  ) |>
  arrange(ProjectKey, PlotID, year)

message("plot-years: ", nrow(plot_year),
        " | from visits: ", nrow(visits),
        " | with >1 visit: ", sum(plot_year$n_visits > 1))

# first-hit classes stay mutually exclusive under averaging
fh_total_py <- rowSums(plot_year[str_subset(cover_cols, "^fh_")], na.rm = TRUE)
stopifnot(all(fh_total_py <= 100 + 1e-6))

write_csv(plot_year,   file.path(ldc_processed_dir, "cover_by_plot_year.csv"))
