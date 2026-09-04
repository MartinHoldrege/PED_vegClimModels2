# 02_ldc_build_code_lookup.R
#
# Build the code -> plant attribute lookup needed to classify LPI hits.
#
# Fetches the LDC geoSpecies table (growth habit, duration) for all plots,
# then tallies every code observed in the cached LPI batches and joins the
# two on (SpeciesKey, code). Codes are only unique within a species list,
# so SpeciesKey is part of the join key.
#
# Outputs:
#   species_attr.rds      - distinct (SpeciesKey, code) plant attributes
#   ldc_code_residue.csv  - codes with no match, for hand assignment
#
# The residue should be mostly non-plant surface codes (soil, rock, litter,
# water, moss, crust). Rerunning refreshes hit counts without overwriting
# assignments already entered in the csv.
#
# geoSpecies fetching is cached per batch and resumable, as in 01_ldc_fetch.R.
# 
# September, 2026

source("Functions/init.R")

ldc_dir <- file.path(paths$large, "Data_raw/LandscapeDataCommonsDat")
lpi_dir <- file.path(ldc_dir, "lpi_batches")
ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")

out_csv <- file.path(ldc_processed_dir, "/ldc_code_residue.csv")

hdr_all <- readRDS(file.path(ldc_dir, "header_raw.rds"))

# ---- geoSpecies for all plots, cached, same batching as LPI -------------
sp_dir <- file.path(ldc_dir, "geospecies_batches")
dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ldc_processed_dir, showWarnings = FALSE, recursive = TRUE)

#' Fetch and cache geoSpecies for one batch of plots
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

all_keys <- sort(unique(hdr_all$PrimaryKey))
iwalk(split(all_keys, ceiling(seq_along(all_keys) / 200)), fetch_sp_batch)

# ---- distinct code attributes -------------------------------------------

distinct_species <- function(path) {
  df <- readRDS(path)
  if(!isTRUE(nrow(df) > 0)) return(NULL)
  
  df |> 
    mutate(across(c(GrowthHabit, GrowthHabitSub, Duration),
                  .f = \(x) na_if(x, "<NA>"))) |> 
    distinct(SpeciesKey, Species, GrowthHabit,
              GrowthHabitSub, Duration)
    
}
sp_attr0 <- list.files(sp_dir, pattern = "^sp_batch-\\d+\\.rds$", full.names = TRUE) |>
  map(distinct_species) |>
  list_rbind() |>
  distinct()


sp_attr <- sp_attr0 |>
  group_by(SpeciesKey, Species) |> 
  # remove NAs when duplicates have not NAs
  filter(!(length(Species) > 1 & is.na(GrowthHabit) & any(!is.na(GrowthHabit)))) |> 
  ungroup()

# a code should have one growth habit within a species list; check
sp_conflict <- sp_attr |>
  summarise(n = n_distinct(GrowthHabitSub), .by = c(SpeciesKey, Species)) |>
  filter(n > 1)

if(nrow(sp_conflict) > 0) {
  message("codes with conflicting growth habit: ", nrow(sp_conflict))
}
  
saveRDS(sp_attr, file.path(ldc_dir, "species_attr.rds"), compress = "xz")

# ---- observed LPI codes, with SpeciesKey attached ------------------------
layer_order <- c("TopCanopy", "Lower1", "Lower2", "Lower3",
                 "Lower4", "Lower5", "SoilSurface")

hdr_key <- hdr_all |> as_tibble() |> distinct(PrimaryKey, SpeciesKey)

#' Tally codes in one LPI batch, split by first-hit status
#' @param path Path to a cached LPI batch.
tally_codes <- function(path) {
  lpi <- readRDS(path) |>
    mutate(code = str_trim(code)) |>
    filter(!is.na(code), code != "") |>
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

code_tally <- list.files(lpi_dir, pattern = "^lpi_batch-\\d+\\.rds$", full.names = TRUE) |>
  map(tally_codes) |>
  list_rbind() |>
  summarise(n_hits  = sum(n_hits, na.rm = TRUE),
            n_first = sum(n_first, na.rm = TRUE),
            .by = c(SpeciesKey, code))

# ---- join and isolate the residue ---------------------------------------
sp_code <- sp_attr |>
  count(Species, GrowthHabitSub, name = "n_lists") |>
  slice_max(n_lists, n = 1, with_ties = FALSE, by = Species) |>
  select(code = Species, GrowthHabitSub)

code_joined <- code_tally |>
  summarise(n_hits = sum(n_hits, na.rm = TRUE),
            n_first = sum(n_first, na.rm = TRUE),
            .by = code) |>
  left_join(sp_code, by = "code")

residue <- code_joined |>
  filter(is.na(GrowthHabitSub)) |>
  arrange(desc(n_first))

message("\ncodes total: ", nrow(code_joined),
        " | matched: ", sum(!is.na(code_joined$GrowthHabitSub)),
        " | residue: ", nrow(residue))
message("first-hits unmatched: ",
        round(100 * sum(residue$n_first, na.rm = TRUE) /
                sum(code_joined$n_first, na.rm = TRUE), 2), "%")

print(residue |> select(SpeciesKey, code, n_hits, n_first) |> slice_max(n_first, n = 50), n = 50)

# ---- write residue for hand assignment, preserving prior work ------------
blank <- residue |>
  select(SpeciesKey, code, n_hits, n_first) |>
  mutate(class = NA_character_, notes = NA_character_)

if (file.exists(out_csv)) {
  prior <- read_csv(out_csv, show_col_types = FALSE) |>
    select(SpeciesKey, code, class, notes)
  blank <- blank |> select(-c(class, notes)) |>
    left_join(prior, by = c("SpeciesKey", "code"))
}

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write_csv(blank, out_csv)