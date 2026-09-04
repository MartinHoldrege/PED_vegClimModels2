# 02_ldc_build_code_lookup.R
#
# Build the code -> plant attribute lookup needed to classify LPI hits.
#
# Fetches the LDC geoSpecies table (growth habit, duration) for all plots,
# tallies every code observed in the cached LPI batches, and attaches a
# growth-habit class to each observed code.
#
# Codes are only guaranteed unique within a species list, so four lookups
# are built and applied in order (first non-missing wins):
#   sp_class       per-(SpeciesKey, code) - authoritative species-list habit
#   global_habit   per-code              - fallback where a list is missing
#                                          an entry for a code it recorded
#   manual_codes   per-code, in-script   - hand classifications for AIM surface
#                                          codes and coarse functional groups
#   generic_lookup per-code              - terradactyl::generic.species growth
#                                          habits for generic/cryptogam codes
#
# Manual code classifications live in the `manual_codes` table below rather
# than in an external file, so this script is the single source of truth.
# Extend that table using the still-missing report it prints and writes.
#
# Outputs:
#   species_attr.rds      - distinct (SpeciesKey, code) plant attributes,
#                           including ScientificName
#   code_class.rds        - every observed (SpeciesKey, code) with its class,
#                           source (list/global/manual/generic), and
#                           ScientificName (for downstream C3/C4 + leaf type)
#   ldc_missing_codes.csv - report of codes still matching nothing, for
#                           triage; output only, never read back
#
# geoSpecies fetching is cached per batch and resumable, as in 01_ldc_fetch.R.
#
# September, 2026

source("Functions/init.R")

# terradactyl::generic.species provides authoritative growth habits for the
# generic / non-plant codes (cryptogams, generic plant lifeforms).
stopifnot(requireNamespace("terradactyl", quietly = TRUE))

ldc_dir <- file.path(paths$large, "Data_raw/LandscapeDataCommonsDat")
lpi_dir <- file.path(ldc_dir, "lpi_batches")
sp_dir  <- file.path(ldc_dir, "geospecies_batches")
ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")

out_csv <- file.path(ldc_processed_dir, "ldc_missing_codes.csv")

# Tallying codes across all LPI batches is slow. When rerun_code_tally is FALSE
# the cached result is read from disk instead (set TRUE to recompute, e.g. after
# the LPI batches change).
rerun_code_tally <- FALSE   # re-tally after the None-sentinel fix; set FALSE once the cache is refreshed
code_tally_cache <- file.path(ldc_processed_dir, "code_tally.rds")

dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ldc_processed_dir, showWarnings = FALSE, recursive = TRUE)

hdr_all <- readRDS(file.path(ldc_dir, "header_raw.rds"))


# ---- helper functions ----------------------------------------------------

replace_with_na <- function(df) {
  mutate(df, across(where(is.character), 
         .fns = \(x) na_if(x, "<NA>")))
}
#' Standardise a growth-habit string
#'
#' Lowercases and strips spaces and hyphens so that variants like
#' "Sub-shrub" and "Forb/herb" compare equal across species lists.
#' NA is returned unchanged.
#'
#' @param x Character vector of raw GrowthHabit or GrowthHabitSub values.
clean_strings <- function(x) {
  x |>
    str_to_lower() |>
    str_remove_all("[\\s-]") |>
    str_replace_all("forb/herb", "forb") |>
    # blank -> NA so habit_class() can fall back to GrowthHabit (e.g. codes
    # with an empty GrowthHabitSub but GrowthHabit == "Nonvascular")
    na_if("")
}

# Classes that habit_class() can return, in the order used to resolve codes
# that appear as more than one class across species lists. Earlier entries
# win. Woody beats non-woody and tree beats shrub; the modelled classes all
# beat the non-target and unknown ones, which only win when nothing else is
# available. Every class habit_class() returns must appear here.
habit_precedence <- c("tree", "shrub", "graminoid", "forb",
                      "vine", "woody_unknown", "nonwoody_unknown",
                      "moss", "lichen", "nonvascular")

#' Collapse a cleaned growth habit to a modelling class
#'
#' Subshrubs and succulents lump with shrubs, sedges and rushes with
#' graminoids. Non-target lifeforms get explicit classes rather than NA so
#' that what is excluded is counted rather than silently dropped. Where
#' GrowthHabitSub is missing, GrowthHabit still distinguishes woody from
#' non-woody, which is not enough to assign a class but is enough to know
#' the code is a plant.
#'
#' @param sub Cleaned GrowthHabitSub.
#' @param habit Cleaned GrowthHabit, used only where `sub` is missing.
habit_class <- function(sub, habit) {
  case_when(
    sub == "tree"                                  ~ "tree",
    sub %in% c("shrub", "subshrub", "succulent")   ~ "shrub",
    sub %in% c("graminoid", "sedge", "rush")       ~ "graminoid",
    sub == "forb"                                  ~ "forb",
    sub == "vine"                                  ~ "vine",
    sub %in% c("moss", "liverwort", "bryophyte",
               "hornwort")                         ~ "nonvascular",
    sub %in% c("lichen", "lichenous")              ~ "nonvascular",
    sub %in% c("nonvascular", "alga", "fungus")    ~ "nonvascular",
    is.na(sub) & habit == "woody"                  ~ "woody_unknown",
    is.na(sub) & habit == "nonwoody"               ~ "nonwoody_unknown",
    is.na(sub) & habit == "nonvascular"            ~ "nonvascular",
    .default = NA_character_
  )
}

#' Flag rows where GrowthHabit contradicts GrowthHabitSub
#'
#' GrowthHabit is documented as derived from GrowthHabitSub, so mismatches
#' are data-entry errors in the source species list rather than real
#' disagreement about the plant.
#'
#' @param sub,habit Cleaned GrowthHabitSub and GrowthHabit.
habit_inconsistent <- function(sub, habit) {
  expected <- case_when(
    sub %in% c("tree", "shrub", "subshrub", "succulent")      ~ "woody",
    sub %in% c("forb", "graminoid", "sedge", "rush")          ~ "nonwoody",
    sub %in% c("nonvascular", "moss", "lichen", "lichenous",
               "liverwort", "bryophyte", "hornwort",
               "alga", "fungus")                              ~ "nonvascular",
    .default = NA_character_
  )
  !is.na(expected) & !is.na(habit) & expected != habit & 
    # this is an acceptable case
    !(habit == 'nonwoody' & expected == 'nonvascular')
}

#' Build a one-row-per-code growth-habit lookup
#'
#' Collapses per-species-list attributes into a single global lookup,
#' resolving codes that appear as more than one class by `precedence`.
#' Most conflicts are a multi-valued PLANTS habit (e.g. "Forb/herb,
#' Subshrub") collapsed differently by different lists, so one documented
#' rule is preferable to a per-list majority vote.
#'
#' @param x Data frame with columns code and class.
#' @param precedence Ordered class vector; earlier classes win conflicts.
#' @param overrides Optional data frame with columns code and class, applied
#'   after the precedence rule.
#' @return One row per code: code, class, n_class, classes (all classes seen,
#'   for auditing), resolved_by ("single", "precedence" or "override").
build_global_habit <- function(x, precedence = habit_precedence,
                               overrides = NULL) {
  
  stopifnot(all(c("code", "class") %in% names(x)))
  
  unknown <- setdiff(na.omit(unique(x$class)), precedence)
  if (length(unknown) > 0) {
    stop("class(es) missing from precedence: ", paste(unknown, collapse = ", "))
  }
  
  out <- x |>
    filter(!is.na(code), !is.na(class)) |>
    distinct(code, class) |>
    summarise(n_class = n_distinct(class),
              classes = paste(sort(unique(class)), collapse = "/"),
              class   = precedence[min(match(class, precedence))],
              .by = code) |>
    mutate(resolved_by = if_else(n_class == 1, "single", "precedence"))
  
  if (!is.null(overrides)) {
    stopifnot(all(c("code", "class") %in% names(overrides)),
              !any(duplicated(overrides$code)))
    out <- out |>
      left_join(select(overrides, code, class_override = class), by = "code") |>
      mutate(resolved_by = if_else(!is.na(class_override), "override", resolved_by),
             class       = coalesce(class_override, class)) |>
      select(-class_override)
  }
  
  stopifnot(!any(duplicated(out$code)))
  out
}

#' Read the distinct species attributes from one cached geoSpecies batch
#'
#' @param path Path to a cached geoSpecies batch. Returns NULL for empty
#'   batches (plots whose project has no species list).
distinct_species <- function(path) {
  df <- readRDS(path)
  if (!isTRUE(nrow(df) > 0)) return(NULL)
  
  df |>
    mutate(across(c(GrowthHabit, GrowthHabitSub, Duration, ScientificName),
                  .f = \(x) na_if(na_if(x, "<NA>"), ""))) |>
    distinct(SpeciesKey, Species, ScientificName,
             GrowthHabit, GrowthHabitSub, Duration)
}

#' Fetch and cache geoSpecies for one batch of plots
#'
#' @param keys Character vector of PrimaryKeys.
#' @param batch_id Batch number used in the file name.
fetch_sp_batch <- function(keys, batch_id) {
  out <- file.path(sp_dir, sprintf("sp_batch-%05d.rds", as.integer(batch_id)))
  if (file.exists(out)) return(invisible(NULL))
  
  x <- try(trex::fetch_ldc(
    data_type = "geoSpecies",
    query_parameters = list("PrimaryKey" = list("=" = keys)),
    take = 10000, delay = 500, timeout = 1800
  ), silent = TRUE)
  
  if (inherits(x, "try-error")) {
    message("FAILED batch ", batch_id)
    return(invisible(NULL))
  }
  
  tmp <- paste0(out, ".partial")
  saveRDS(x, tmp)
  file.rename(tmp, out)
  invisible(NULL)
}

#' Tally codes in one LPI batch, split by first-hit status
#'
#' @param path Path to a cached LPI batch.
# No-hit sentinels: codes meaning "nothing recorded" rather than a real hit.
# Dropped BEFORE first-hit ranking so a bare pin's true first hit is its
# SoilSurface code, not a "None" placeholder sitting in TopCanopy. Matched
# case-insensitively (both "None" and "NONE" occur). "N"/"NA" are deliberately
# left out - they do not occur in the data and could be real codes elsewhere.
no_hit_codes <- c("none", "")

tally_codes <- function(path) {
  lpi <- readRDS(path) |>
    mutate(code = str_trim(code)) |>
    filter(!is.na(code), !str_to_lower(code) %in% no_hit_codes) |>
    left_join(hdr_key, by = "PrimaryKey")
  
  first <- lpi |>
    mutate(rank = match(layer, layer_order)) |>
    filter(!is.na(rank)) |>
    slice_min(rank, n = 1, with_ties = FALSE,
              by = c(PrimaryKey, LineKey, PointNbr))
  
  full_join(
    lpi   |> count(SpeciesKey, code, name = "n_hits"),
    first |> count(SpeciesKey, code, name = "n_first"),
    by = c("SpeciesKey", "code")
  )
}


# ---- geoSpecies for all plots, cached, same batching as LPI --------------
all_keys <- sort(unique(hdr_all$PrimaryKey))
iwalk(split(all_keys, ceiling(seq_along(all_keys) / 200)), fetch_sp_batch)

sp_attr0 <- list.files(sp_dir, pattern = "^sp_batch-\\d+\\.rds$", full.names = TRUE) |>
  map(distinct_species) |>
  list_rbind() |>
  distinct()

# where a (SpeciesKey, code) appears both with and without attributes, keep
# only the populated rows
sp_attr <- sp_attr0 |>
  group_by(SpeciesKey, Species) |>
  filter(!(length(Species) > 1 & is.na(GrowthHabit) & any(!is.na(GrowthHabit)))) |>
  ungroup()

saveRDS(sp_attr, file.path(ldc_dir, "species_attr.rds"), compress = "xz")

# ---- tally observed LPI codes -------------------------------------------
# LPI layers, top to bottom. Include Lower1..Lower10 (not just Lower5): the
# data contains deeper layers (e.g. Lower6, Lower7) that would otherwise get
# rank NA and be silently dropped from the first-hit calculation.
layer_order <- c("TopCanopy", paste0("Lower", 1:10), "SoilSurface")

hdr_key <- hdr_all |> as_tibble() |> distinct(PrimaryKey, SpeciesKey)

if (!rerun_code_tally && file.exists(code_tally_cache)) {
  message("code_tally: reading cached ", code_tally_cache)
  code_tally <- readRDS(code_tally_cache)
} else {
  message("code_tally: tallying LPI batches (slow) ...")
  code_tally <- list.files(lpi_dir, pattern = "^lpi_batch-\\d+\\.rds$", full.names = TRUE) |>
    map(tally_codes) |>
    list_rbind() |>
    summarise(n_hits  = sum(n_hits, na.rm = TRUE),
              n_first = sum(n_first, na.rm = TRUE),
              .by = c(SpeciesKey, code))
  saveRDS(code_tally, code_tally_cache, compress = "xz")
  message("code_tally: cached to ", code_tally_cache)
}


# ---- growth-habit lookups ------------------------------------------------
sp_lookup <- sp_attr |>
  select(SpeciesKey, code = Species, GrowthHabit, GrowthHabitSub) |>
  mutate(sub   = clean_strings(GrowthHabitSub),
         habit = clean_strings(GrowthHabit),
         class = habit_class(sub, habit),
         bad   = habit_inconsistent(sub, habit))

message("rows where GrowthHabit contradicts GrowthHabitSub: ", sum(sp_lookup$bad))

# per-list lookup; a few lists still disagree internally, resolved the same way
sp_class <- sp_lookup |>
  filter(!is.na(class)) |>
  distinct(SpeciesKey, code, class) |>
  slice_min(match(class, habit_precedence), n = 1, with_ties = FALSE,
            by = c(SpeciesKey, code))

global_habit <- build_global_habit(sp_lookup)

print(count(global_habit, resolved_by))
print(count(global_habit, class))

# every conflict pattern and how it was resolved, for the record
global_habit |>
  filter(n_class > 1) |>
  count(classes, class, sort = TRUE) |> 
  print()

# ---- manual code table --------
# Hand classifications for codes that neither lookup resolves: non-plant
# surface codes and the coarse functional-group codes from plots recorded
# without a species list. Extend this as new codes show up in the
# still-missing report at the end. Manual entries only fill codes left
# unclassified by the list/global lookups - they never override a species
# classification.
#
# The coarse functional groups can only be classed at the top level (NOT the
# finer needleleaf/broadleaf or C3/C4 splits, which need species-level data).
# "herbaceous" = graminoid or forb, unresolved, handled explicitly downstream.
manual_codes <- tibble::tribble(
  ~code,                             ~class,        ~notes,
  # coarse functional groups (plots recorded without a species list)
  "Trees",                           "tree",        "coarse; no needleleaf/broadleaf split",
  "Shrubs",                          "shrub",       NA_character_,
  "Perennial grasses",               "graminoid",   "coarse; no C3/C4 split",
  "Annual plants",                   "herbaceous",  "coarse; grass/forb unknown",
  "Sub-shrubs and perennial forbs",  "forb",        NA_character_,
  "Plant base",                      NA_character_, "hit artifact; assign or drop",
  'DEAD SHRUB',                       "dead",        NA_character_,
  'DEAD BUSH',                       "dead",        NA_character_,
  'DEAD TREE',                       "dead",        NA_character_,
  # non-plant surface codes NOT in terradactyl::generic.species. Rock cover
  # (rock, cobble, gravel, bedrock, stone) is grouped as bare_ground
  # The 2XXX / cryptogam generics (2ALGA, CY,
  # M, LC, ...) are handled by the terradactyl tier below, not here.
  "S",                               "bare_ground", "soil",
  "R",                               "bare_ground", "rock",
  "GR",                              "bare_ground", "gravel",
  "CB",                              "bare_ground", "cobble",
  "ST",                              "bare_ground", "stone",
  "BR",                              "bare_ground", "bedrock",
  "RF",                              "bare_ground", "rock fragments",
  "WA",                              "water",       "water",
  "L",                               "litter",      NA_character_,
  "HL",                              "litter",      "herbaceous litter",
  "WL",                              "litter",      "woody litter",
  "EL",                              "litter",      "embedded litter (AIM data form)",
  "D",                               "litter",      "duff / decomposed litter (AIM data form)",
  # None / NONE are no-hit sentinels dropped in tally_codes, so they never
  # reach this table.
  # PC/LM/BY occur only as SoilSurface hits (PC, LM from NWERN DIMA; BY from
  # TerrADat) but their meanings are unconfirmed; left NA pending review so
  # they are tracked, not miscoded.
  "PC",                              NA_character_, "SoilSurface (NWERN DIMA); confirm meaning",
  "LM",                              NA_character_, "SoilSurface (NWERN DIMA); confirm meaning",
  "BY",                              NA_character_, "SoilSurface (TerrADat); confirm meaning"
)
manual_codes <- manual_codes |> mutate(code = str_trim(code))
stopifnot(!any(duplicated(manual_codes$code)))

# LPI codes are only str_trim()'d, never clean_strings()'d, so manual codes
# must match the raw (trimmed) code exactly. Warn on any that match nothing -
# usually a typo or a case/spacing difference from the observed code.
unseen <- setdiff(manual_codes$code, code_tally$code)
if (length(unseen) > 0) {
  message("manual_codes not seen in the LPI tally (typo/case/spacing?): ",
          paste(unseen, collapse = ", "))
}


# ---- terradactyl generic-code lookup -------------------------------------
# terradactyl::generic.species gives growth habits for generic / non-plant
# codes (cryptogams, generic plant lifeforms) under both a 2-prefixed `Code`
# and an unprefixed `Prefix`. Run these through the same habit_class() used
# for the species data, then key on both forms. Applied AFTER manual_codes so
# the collision Prefix "S" (= Shrub here) can never override S = bare_ground.
generic_species <- local({
  e <- new.env()
  utils::data("generic.species", package = "terradactyl", envir = e)
  e$generic.species
})

generic_class <- generic_species |>
  transmute(Code, Prefix,
            class = habit_class(clean_strings(GrowthHabitSub),
                                clean_strings(GrowthHabit))) |>
  filter(!is.na(class))

generic_lookup <- bind_rows(
  transmute(generic_class, code = Code,   class),
  transmute(generic_class, code = Prefix, class)
) |>
  filter(!is.na(code)) |>
  distinct(code, class) |>
  slice_min(match(class, habit_precedence), n = 1, with_ties = FALSE, by = code)


# One ScientificName per (SpeciesKey, code), carried into code_class.rds for
# downstream species-level traits (C3/C4 pathway, needle/broad leaf). A few
# pairs carry more than one name; keep the first. Non-plant / functional /
# unmatched codes have no name and get NA.
sp_names <- sp_attr |>
  filter(!is.na(ScientificName)) |>
  distinct(SpeciesKey, code = Species, ScientificName) |>
  slice_head(n = 1, by = c(SpeciesKey, code))


# ---- attach a class to every observed (SpeciesKey, code) -----------------
# Precedence (first non-missing wins): species-list habit, global per-code
# habit, the manual table, then the terradactyl generic lookup. manual before
# generic so our surface codes (S = bare_ground, etc.) win any Prefix clash.
code_joined <- code_tally |>
  left_join(sp_class, by = c("SpeciesKey", "code")) |>
  left_join(select(global_habit, code, class_global = class), by = "code") |>
  left_join(select(manual_codes, code, class_manual = class), by = "code") |>
  left_join(rename(generic_lookup, class_generic = class), by = "code") |>
  mutate(source = case_when(!is.na(class)         ~ "list",
                            !is.na(class_global)  ~ "global",
                            !is.na(class_manual)  ~ "manual",
                            !is.na(class_generic) ~ "generic",
                            .default              = "unmatched"),
         class  = coalesce(class, class_global, class_manual, class_generic)) |>
  select(-class_global, -class_manual, -class_generic) |>
  left_join(sp_names, by = c("SpeciesKey", "code")) |> 
  replace_with_na()

code_joined |>
  summarise(n_first = sum(n_first, na.rm = TRUE), .by = source) |>
  mutate(pct = 100 * n_first / sum(n_first)) |>
  print()


# ---- still-missing codes -------------------------------------------------
# Codes matching none of the four lookups. Collapsed to a global `code` key
# (still-missing codes are dominated by surface/unknown codes that are global)
# and sorted by first-hit frequency so the biggest gaps triage first. Add any
# that should be classified to `manual_codes` above and re-run; the remainder
# are genuinely unclassifiable species or unknown codes.
missing_codes <- code_joined |>
  filter(is.na(class)) |>
  summarise(n_hits  = sum(n_hits, na.rm = TRUE),
            n_first = sum(n_first, na.rm = TRUE),
            n_lists = n_distinct(SpeciesKey),
            .by = code) |>
  arrange(desc(n_first)) |> 
  replace_with_na()

message("\n(SpeciesKey, code) pairs: ", nrow(code_joined),
        " | classified: ", sum(!is.na(code_joined$class)),
        " | still-missing codes: ", nrow(missing_codes))
message("first-hits still unclassified: ",
        round(100 * sum(missing_codes$n_first, na.rm = TRUE) /
                sum(code_joined$n_first, na.rm = TRUE), 2), "%")

print(missing_codes |> slice_max(n_first, n = 50))


# ---- outputs -------------------------------------------------------------
# Final per-(SpeciesKey, code) classification, with ScientificName for the
# downstream C3/C4 and leaf-type trait build, for LPI hit classing.
saveRDS(code_joined, file.path(ldc_processed_dir, "code_class.rds"), compress = "xz")

# Report of still-missing codes. Output only - never read back; manual
# classifications live in `manual_codes` above, the single source of truth.
write_csv(missing_codes, out_csv)

