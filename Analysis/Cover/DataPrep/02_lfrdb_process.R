# 02_lfrdb_process.R
#
# Summarise LFRDB cover into the same functional groups as the LDC pipeline.
#
# HOW LFRDB COVER RELATES TO THE LDC COVERS
#
#   src_*  the crew's own lifeform call, "accounting for most overlap"
#          (dtStands Source*Cov). A direct field observation rather than a
#          derived quantity, but often missing.
#   sum_*  species LFAbsCov summed within group. Overlap NOT accounted for,
#          may exceed 100%. The only option for the component groups.
#   adj_*  LANDFIRE's overlap-adjusted version of sum_*, bounded 0-100%
#          (dtStands LF*CovAdj). Method undocumented; the check below tests
#          whether it is the standard random-overlap transform, in which case
#          it is a deterministic function of sum_* and adds no information.
# the raw data folders contain the data dictionary description pdfs 
#
# Summation is unreliable for absolute cover but fine for the proportions:
# the overlap bias falls on numerator and denominator alike, so C3/C4/forb
# and needle/broad fractions are usable from summed species cover.
#
# LFRDB has NO bare ground. Source/LFNvasCov are non-vascular plant cover,
# not exposed soil.
#
# created rds files in 01_lfrdb_read.R
# Inputs:  lfrdb_stands.rds, lfrdb_species.rds, Osborne Table S3 .xls
# Outputs: lfrdb_cover_by_plot.csv, lfrdb_graminoid_pathway.csv,
#          lfrdb_tree_leaf_type.csv
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source("Functions/data/plant_traits.R")


# params ------------------------------------------------------------------

# Group columns the output must always carry, so the schema is stable even
# if a group is absent from the data.
component_groups <- c("tree_needle", "tree_broad", "tree_unknown",
                      "grass_c3", "grass_c4", "grass_unknown", "forb")
toplevel_groups  <- c("tree", "shrub", "herbaceous")
all_groups       <- c(toplevel_groups, component_groups)


year_min <- 2000

# read in data ------------------------------------------------------------


lfrdb_processed_dir <- file.path(paths$large, "Data_processed/LANDFIRE_LFRDB")
osborne_xls <- file.path(paths$large, "Data_raw/plant_spp_info",
                         "Osborne_etal_2014_newphytologist_tables3_c4-genera.xls")

stands  <- readRDS(file.path(lfrdb_processed_dir, "lfrdb_stands.rds")) |> 
  mutate(YYYY = as.numeric(YYYY)) |> 
  filter(YYYY >= year_min)
species <- readRDS(file.path(lfrdb_processed_dir, "lfrdb_species.rds")) |> 
  filter(EventID %in% stands$EventID)

osborne <- read_osborne(osborne_xls)
known_genera <- unique(c(osborne$genus, non_poaceae_pathway$genus, 
                         .needle_genera))


# ---- map LFRDB lifeform codes to the LDC classes -------------------------
# Lifeform (dtSpecies): F = forb, G = graminoid, H = herb, N = nonvascular,
# S = shrub, T = tree, V = vine. "H" is herbaceous of unspecified type and
# cannot be split into forb or graminoid - it counts toward herbaceous but
# not toward the forb / C3 / C4 fractions.
#
# LFRDB has no subshrub or succulent code; both fall under S. That matches
# the LDC decision to lump subshrubs and succulents with shrubs, so no
# reconciliation is needed.
lifeform_class <- tibble::tribble(
  ~Lifeform, ~class,
  "T",       "tree",
  "S",       "shrub",
  "G",       "graminoid",
  "F",       "forb",
  "H",       "herbaceous",
  "N",       "nonvascular",
  "V",       "vine"
)

sp <- species |>
  filter(!is.na(LFAbsCov)) |>
  left_join(lifeform_class, by = "Lifeform") |>
  mutate(genus = clean_genus(SciName, known_genera))

# An unmapped Lifeform means lifeform_class above needs extending; those rows
# are silently excluded from every group. <1% of cover, so not trying to solve
unmapped <- sp |> filter(is.na(class)) |> count(Lifeform, sort = TRUE)
if (nrow(unmapped) > 0) {
  message("\nUNMAPPED Lifeform codes - extend lifeform_class:")
  print(as.data.frame(unmapped))
}

# "H" cannot be split into forb or graminoid. If it carries a large share of
# herbaceous cover, the forb and C3/C4 fractions rest on a small subset.
herb_share <- sp |>
  filter(class %in% c("graminoid", "forb", "herbaceous")) |>
  summarise(cover = sum(LFAbsCov, na.rm = TRUE), .by = class) |>
  mutate(pct = round(100 * cover / sum(cover), 1))
message("\nherbaceous cover by lifeform class ('herbaceous' = unsplittable H):")
print(as.data.frame(herb_share))


# ---- graminoid pathway ---------------------------------------------------

gram_base <- sp |>
  filter(class == "graminoid") |>
  summarise(cover_sum = sum(LFAbsCov, na.rm = TRUE), n_rows = n(),
            .by = c(Item, genus)) |>
  arrange(is.na(genus), desc(cover_sum)) |>
  # One genus per code
  # in case one Item maps to multiple genera (but doesn't look to happen in this dataset)
  slice_head(n = 1, by = Item) 


gram_lookup <- gram_base |> 
  left_join(assign_pathway(gram_base$genus, osborne), by = 'genus')

check_unique(gram_lookup, Item, label = "gram_lookup")

flag_genera(gram_lookup, "c3_by_default",
            "Graminoid genera assigned C3 by elimination - scan for C4 misses:",
            weight_col = "cover_sum", n = 50)

flag_genera(gram_lookup, "no_scientific_name",
            "Graminoid codes with no usable genus:",
            weight_col = "cover_sum")


# ---- tree leaf type ------------------------------------------------------
tree_base <- sp |>
  filter(class == "tree") |>
  summarise(cover_sum = sum(LFAbsCov, na.rm = TRUE), n_rows = n(),
            .by = c(Item, genus)) |>
  arrange(is.na(genus), desc(cover_sum)) |>
  slice_head(n = 1, by = Item)

tree_lookup <- tree_base |> 
  left_join(assign_leaf_type(tree_base$genus), by = 'genus')
check_unique(tree_lookup, Item, label = "tree_lookup")

flag_genera(tree_lookup, "broad_by_default",
            "Tree genera assigned broadleaf by elimination - scan for conifers:",
            source_col = "leaf_type_source", weight_col = "cover_sum", n = 50)

# ---- long code -> group lookup, matching the LDC groups ------------------

# A species code can carry different Lifeform assignments across contributing
# programs (114 codes, 2.3% of cover). Resolve to one class per code by the
# number of source programs using it, breaking ties on plot count: a code
# recorded by many independent programs as one lifeform is more reliable than
# one program's heavier cover totals (e.g. Typha is graminoid in 49 programs
# but forb in 4, which hold most of its cover). NA sorted last so it only
# wins where no program assigned a lifeform.
sp_groups <- sp |>
  left_join(select(stands, EventID, SourceID), by = "EventID") |>
  summarise(n_sources = n_distinct(SourceID),
            n_plots   = n_distinct(EventID),
            .by = c(Item, class)) |>
  arrange(is.na(class), desc(n_sources), desc(n_plots)) |>
  slice_head(n = 1, by = Item) |>
  select(Item, class) |>
  left_join(select(gram_lookup, Item, pathway),   by = "Item") |>
  left_join(select(tree_lookup, Item, leaf_type), by = "Item")

# a code carrying two classes across plots would double-count on the join
check_unique(sp_groups, Item, label = "sp_groups")

code_group <- bind_rows(
  sp_groups |> filter(class == "tree")  |> transmute(Item, group = "tree"),
  sp_groups |> filter(class == "shrub") |> transmute(Item, group = "shrub"),
  sp_groups |> filter(class %in% c("graminoid", "forb", "herbaceous")) |>
    transmute(Item, group = "herbaceous"),
  sp_groups |> filter(class == "tree", leaf_type == "needle") |>
    transmute(Item, group = "tree_needle"),
  sp_groups |> filter(class == "tree", leaf_type == "broad") |>
    transmute(Item, group = "tree_broad"),
  sp_groups |> filter(class == "tree", !leaf_type %in% c("needle", "broad")) |>
    transmute(Item, group = "tree_unknown"),
  sp_groups |> filter(class == "graminoid", pathway == "C3") |>
    transmute(Item, group = "grass_c3"),
  sp_groups |> filter(class == "graminoid", pathway == "C4") |>
    transmute(Item, group = "grass_c4"),
  sp_groups |> filter(class == "graminoid", !pathway %in% c("C3", "C4")) |>
    transmute(Item, group = "grass_unknown"),
  sp_groups |> filter(class == "forb") |> transmute(Item, group = "forb")
) |>
  distinct()

# a group name not in all_groups would produce an unexpected output column
bad_group <- setdiff(unique(code_group$group), all_groups)
if (length(bad_group) > 0) {
  stop("code_group has unexpected group(s): ", paste(bad_group, collapse = ", "))
}

message("\ncodes per group:")
count(code_group, group) |> as.data.frame() |> print()


# ---- sum_: summed species cover by group ---------------------------------

# these can sum to >100 per group, b/ summing individual species
# cover estimates
sum_cover <- sp |>
  inner_join(code_group, by = "Item", relationship = "many-to-many") |>
  summarise(cover = sum(LFAbsCov, na.rm = TRUE), .by = c(EventID, group)) |>
  pivot_wider(names_from = group, values_from = cover,
              names_prefix = "sum_", values_fill = 0) |>
  ensure_cols(paste0("sum_", all_groups))

# ---- adj_ and src_: stand-level lifeform cover ---------------------------

adj_cover <- stands |>
  transmute(EventID,
            adj_tree       = LFTreeCovAdj,
            adj_shrub      = LFShrubCovAdj,
            adj_herbaceous = LFHerbCovAdj)

src_cover <- stands |>
  transmute(EventID,
            src_tree       = SourceTreeCov,
            src_shrub      = SourceShrubCov,
            src_herbaceous = SourceHerbCov)


# ---- assemble one row per plot -------------------------------------------

cover_plot <- stands |>
  transmute(EventID, region, SourceID, Protocol, LFVersion, Type, Purpose,
            LocMeth,
            Longitude_wgs84 = Long, Latitude_wgs84 = Lat,
            year = YYYY) |>
  filter(!is.na(Latitude_wgs84), !is.na(Longitude_wgs84)) |> 
  left_join(src_cover, by = "EventID") |>
  left_join(adj_cover, by = "EventID") |>
  left_join(sum_cover, by = "EventID") |>
  mutate(across(starts_with("sum_"), \(z) replace_na(z, 0)))

check_unique(cover_plot, EventID, label = "cover_plot")

# ---- resolved top-level cover: adj_ preferred, src_ as fallback ----------
# adj_ is derived from species-level LFAbsCov and is continuous; src_ is the
# crew's lifeform call, which for most source programs is an ordinal cover
# class with a midpoint assigned (MIDNR1, the largest contributor, has four
# distinct values across 143k plots). adj_ is therefore preferred where both
# exist, but the two are strongly complementary in coverage, so src_ fills
# the plots adj_ is missing.
#
# src_ is not documented as bounded and contains values above 100; both are
# clamped to [0, 100]. The *_origin columns record which field each value
# came from, for use as a source effect - the two definitions diverge at
# high cover, where adj_ saturates and src_ does not.

#' Clamp a cover percentage to [0, 100]
#' @param x Numeric cover values.
clamp_cover <- function(x) pmin(pmax(x, 0), 100)

cover_plot <- cover_plot |>
  mutate(across(all_of(c(paste0("adj_", toplevel_groups),
                         paste0("src_", toplevel_groups))), clamp_cover))

for (g in toplevel_groups) {
  cover_plot[[paste0("cov_", g)]] <- coalesce(cover_plot[[paste0("adj_", g)]],
                                              cover_plot[[paste0("src_", g)]])
  cover_plot[[paste0("cov_", g, "_origin")]] <- case_when(
    !is.na(cover_plot[[paste0("adj_", g)]]) ~ "adj",
    !is.na(cover_plot[[paste0("src_", g)]]) ~ "src",
    .default = NA_character_
  )
}


# ---- drop src_ cover recorded as coarse ordinal classes ------------------
# Some programs recorded cover in a few ordinal classes, to which LANDFIRE
# assigned midpoints (MIDNR1, 143k plots, has four distinct tree values).
# adj_ is continuous, so only src_-origin values are at risk. Null those for
# the affected class, leaving the plot's other classes intact.
for (g in toplevel_groups) {
  cv <- paste0("cov_", g); og <- paste0("cov_", g, "_origin")
  
  bad <- cover_plot |>
    filter(.data[[og]] == "src", .data[[cv]] > 0) |>
    summarise(n_vals = n_distinct(.data[[cv]]), n = n(), .by = SourceID) |>
    # discarding if only having very coarse 
    filter(n >= 100, n_vals <6) |>
    pull(SourceID)
  
  drop <- cover_plot$SourceID %in% bad & cover_plot[[og]] %in% "src"
  cover_plot[[cv]][drop] <- NA_real_
  cover_plot[[og]][drop] <- NA_character_
  
  message(g, ": nulled ", sum(drop), " ordinal src_ values from ",
          length(bad), " source(s)")
}


# ---- component groups: zero only where the parent class was assessed -----
# sum_ comes from summing species records, so an absent group yields no row
# and pivot_wider fills 0. That is correct where the class WAS assessed and
# simply had no species, and wrong where it was never assessed. Set the
# latter back to NA, using the parent class as the indicator of assessment.
cover_plot <- cover_plot |>
  mutate(
    across(all_of(paste0("sum_", c("tree", "tree_needle", "tree_broad",
                                   "tree_unknown"))),
           \(z) if_else(is.na(cov_tree), NA_real_, z)),
    across(all_of("sum_shrub"),
           \(z) if_else(is.na(cov_shrub), NA_real_, z)),
    across(all_of(paste0("sum_", c("herbaceous", "forb", "grass_c3",
                                   "grass_c4", "grass_unknown"))),
           \(z) if_else(is.na(cov_herbaceous), NA_real_, z))
  )


# ---- checks --------------------------------------------------------------
cov_cols <- paste0("cov_", toplevel_groups)

stopifnot(all(map_lgl(cov_cols, \(cc) {
  z <- cover_plot[[cc]]
  all(is.na(z) | (z >= 0 & z <= 100))
})))

# origin matrix
map_dfr(toplevel_groups, \(g) {
  count(cover_plot, origin = .data[[paste0("cov_", g, "_origin")]]) |>
    mutate(group = g)
}) |>
  pivot_wider(names_from = origin, values_from = n, values_fill = 0) |>
  as.data.frame() |> print()

# "\nNA counts by column after assessment masking:"
cover_plot |>
  summarise(across(starts_with(c("cov_", "sum_")), \(z) sum(is.na(z)))) |>
  glimpse()

message("\nplots retained: ", nrow(cover_plot))

test <- cover_plot |> 
  select(Latitude_wgs84, Longitude_wgs84, year) |> 
  duplicated() |> 
  sum()

if(test > 0) stop('year-location duplicates exist')

# ---- outputs -------------------------------------------------------------
write_csv(cover_plot, file.path(lfrdb_processed_dir, "lfrdb_cover_by_plot.csv"))
# write_csv(select(gram_lookup, Item, genus, pathway, pathway_source, cover_sum),
#           file.path(lfrdb_processed_dir, "lfrdb_graminoid_pathway.csv"))
# write_csv(select(tree_lookup, Item, genus, leaf_type, leaf_type_source, cover_sum),
#           file.path(lfrdb_processed_dir, "lfrdb_tree_leaf_type.csv"))