# 03_ldc_graminoid_pathway.R
#
# Assign C3/C4 to graminoid codes.
#
# Poaceae genera come from Osborne et al. (2014) New Phytologist 204:441-446,
# Supporting Information Table S3 (genus-level). Osborne covers grasses only,
# so Cyperaceae and Juncaceae genera are supplied from a small supplementary
# table (Bruhl & Wilson 2007, Aliso 23:99-148). Note, Osborne misses
# a bunch of new genera but they look to be all c3, so defaulting them
# to c3 in the code is working fine.
#
# Genera containing both C3 and C4 species cannot be assigned at genus level
# and are left NA. The script flags any that actually occur in the data, with
# hit counts, so they can be resolved at species level only if they matter.
#
# Input:  code_class.rds (must carry ScientificName), Osborne Table S3 .xls
# (file created 02_ldc_build_code_lookup)
# Output: graminoid_pathway.csv (code, genus, pathway, pathway_source)
#
# September, 2026

source("Functions/init.R")

ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")
osborne_xls <- file.path(paths$large, "Data_raw/plant_spp_info",
                         "Osborne_etal_2014_newphytologist_tables3_c4-genera.xls")

code_class <- readRDS(file.path(ldc_processed_dir, "code_class.rds")) |> 
  as_tibble()
stopifnot("ScientificName" %in% names(code_class))


# ---- Osborne genus-level pathways (Poaceae) ------------------------------
# Header sits on row 8; rows above are the caption and explanatory notes.
osborne <- readxl::read_excel(osborne_xls, sheet = "Table S1", skip = 7) |>
  rlang::set_names(c("genus_authority", "pathway_raw", "species_sampled",
                     "references", "notes")) |>
  filter(!is.na(genus_authority)) |>
  mutate(genus = word(genus_authority, 1),
         # Only unambiguous single-pathway genera are usable at genus level.
         # Everything else ("C4 & C3", "unknown", ...) stays NA and is flagged.
         pathway = case_when(str_trim(pathway_raw) == "C3" ~ "C3",
                             str_trim(pathway_raw) == "C4" ~ "C4",
                             .default = 'mixed')) |>
  select(genus, pathway, pathway_raw)

stopifnot(nrow(osborne) == 708, !any(duplicated(osborne$genus)))

table(osborne$pathway)

osborne$pathway[osborne$genus == 'Aristida'] <- 'C4' #in north america common sp are all c4
  
# ---- non-Poaceae graminoids ----------------------------------------------
# Osborne covers Poaceae only. Cyperaceae pathways from Bruhl & Wilson (2007)
# Aliso 23:99-148, Table 1. Listed below is every C4 and C3/C4 genus in that
# table, in table order. "mixed" means the genus contains both pathways and
# cannot be assigned at genus level. All other genera in Table 1 are C3, so a
# genus matching neither source is defaulted to C3 - Juncaceae included.
non_poaceae <- tibble::tribble(
  ~genus,             ~pathway_np, ~np_note,
  "Abildgaardia",     "mixed",     "C3/C4; = Fimbristylis in World Checklist",
  "Alinula",          "C4",        NA_character_,
  "Ascolepis",        "C4",        NA_character_,
  "Bulbostylis",      "C4",        NA_character_,
  "Crosslandia",      "C4",        NA_character_,
  "Cyperus",          "mixed",     "C3/C4; most N. American spp. C4",
  "Eleocharis",       "mixed",     "C3/C4",
  "Fimbristylis",     "mixed",     "C3/C4",
  "Kyllinga",         "C4",        "= Cyperus in World Checklist",
  "Lipocarpha",       "C4",        NA_character_,
  "Nelmesia",         "C4",        NA_character_,
  "Nemum",            "C4",        NA_character_,
  "Pycreus",          "C4",        "= Cyperus in World Checklist",
  "Queenslandiella",  "C4",        "= Cyperus in World Checklist",
  "Rhynchospora",     "mixed",     "C3/C4",
  "Sphaerocyperus",   "C4",        NA_character_,
  "Volkiella",         "C4",       NA_character_
)



# ---- genus cleaning ------------------------------------------------------
# ScientificName is not consistently formatted. Three failure modes matter,
# because an unrecognised genus is silently defaulted to C3:
#   1. authority text leaks in       ("Juncus L.", "Vulpia C.C.")
#   2. epithet runs into the genus   ("Bromustectorum", "Poasecunda")
#   3. placeholders instead of names ("Unknown", "Generic", "#N/A")

# Strings carrying no botanical information; set to NA rather than defaulted.
bad_genus <- c("Unknown", "Generic", "Perennial", "Check", "#N/A", "NA")

# Common names appearing in the ScientificName field.
common_names <- tibble::tribble(
  ~genus,   ~genus_fixed,
  "Fescue", "Festuca"
)

# C4 genera absent from Osborne: synonyms, segregate genera, and misspellings.
# Without these they fall through to the C3 default, which is a silent error.
c4_extra <- tibble::tribble(
  ~genus,        ~pathway_np, ~np_note,
  "Dasyochloa",  "C4",        "= Erioneuron",
  "Pleuraphis",  "C4",        "= Hilaria",
  "Monroa",      "C4",        "= Munroa",
  "Vilfa",       "C4",        "= Sporobolus",
  "Chaetochloa", "C4",        "= Setaria",
  "Erianthus",   "C4",        "= Saccharum",
  "Sorgastrum",  "C4",        "misspelling of Sorghastrum",
  "Hopia",       "C4",        "segregate of Panicum",
  "Hypogynium",  "C4",        "Andropogoneae"
)

non_poaceae <- bind_rows(non_poaceae, c4_extra)
stopifnot(!any(duplicated(non_poaceae$genus)))

known_genera <- unique(c(osborne$genus, non_poaceae$genus))

#' Split a genus name that has run together with its epithet
#'
#' Returns the longest known genus the string starts with, or the string
#' unchanged if none matches. Can fail in some cases
#'
#' @param x Character vector of candidate genus strings.
#' @param genera Character vector of known genus names.
fix_runtogether <- function(x, genera = known_genera) {
  map_chr(x, \(s) {
    if (is.na(s) || s %in% genera) return(s)
    hit <- genera[str_starts(s, fixed(genera))]
    if (length(hit) == 0) return(s)
    hit[which.max(nchar(hit))]
  })
}

# testing
if(FALSE) {
  known_genera <- c("Poa", "Bromus", "Festuca", "Elymus", "Eragrostis",
                    "Pseudoroegneria", "Setaria", "Sporobolus", "Carex",
                    "Stipa", "Bromidium", "Panicularia")
  
  x <- c("Bromustectorum", "Poasecunda", "Elymuselymoides",
         "Festucaidahoensis", "Pseudoroegneriaspicata",
         "EragrostiscilianensisAll", "SetariaP", "Sporoboluscryptandrus",
         "Carexfilifolia", "Stipacomata",
         "Poa", "Bromidium", "Panicularia", "Poaceae", "Xyz", NA)
  
  fix_runtogether(x) # fails for 'poaceae'
}

# ---- graminoid codes in the data -----------------------------------------
graminoids_raw <- code_class |>
  filter(class == "graminoid") |>
  mutate(
    # leading alphabetic run only: drops authorities and punctuation
    genus = word(ScientificName, 1) |> str_remove("[^A-Za-z].*$"),
    genus = if_else(genus %in% bad_genus | genus == "", NA_character_, genus)
  ) |>
  left_join(common_names, by = "genus") |>
  mutate(genus = coalesce(genus_fixed, genus),
         genus = fix_runtogether(genus)) |>
  select(-genus_fixed) |>
  summarise(n_hits  = sum(n_hits,  na.rm = TRUE),
            n_first = sum(n_first, na.rm = TRUE),
            .by = c(code, ScientificName, genus))

graminoids <- graminoids_raw |>
  summarise(n_hits = sum(n_hits, na.rm = TRUE),
            n_first = sum(n_first, na.rm = TRUE),
            .by = c(code, genus)) |>
  # genus with the most hits wins; NA sorted last so it only wins if it's
  # the only option (b/ can be mis-spellings etc. of the genus)
  arrange(is.na(genus), desc(n_hits)) |>
  slice_head(n = 1, by = code) |> 
  # totals must cover all rows for the code, not just the winning genus
  left_join(summarise(graminoids_raw,
                      n_hits_tot  = sum(n_hits,  na.rm = TRUE),
                      n_first_tot = sum(n_first, na.rm = TRUE),
                      .by = code),
            by = "code") |>
  transmute(code, genus, n_hits = n_hits_tot, n_first = n_first_tot)

pathway_lookup <- graminoids |>
  left_join(osborne, by = "genus") |>
  left_join(non_poaceae, by = "genus") |>
  mutate(
    pathway_source = case_when(
      !is.na(pathway_raw) ~ "osborne",
      !is.na(pathway_np)  ~ "bruhl_wilson",
      is.na(genus)        ~ "no_scientific_name",
      .default            = "c3_by_default"
    ),
    pathway = case_when(
      !is.na(pathway)    ~ pathway,
      !is.na(pathway_np) ~ pathway_np,
      is.na(genus)       ~ NA_character_,
      # Osborne is complete for Poaceae and Bruhl & Wilson lists every C4 and
      # C3/C4 Cyperaceae genus, so an unlisted genus is C3 by elimination.
      .default           = "C3"
    )
  )

# ---- flags ---------------------------------------------------------------
#' Print a flagged subset of the pathway lookup, by genus and hit count
#'
#' @param x pathway_lookup.
#' @param src pathway_source value to report on.
#' @param msg Message printed above the table.
#' @param n Maximum genera to print.
flag_genera <- function(x, src, msg, n = 30) {
  out <- x |>
    filter(pathway_source == src) |>
    summarise(n_codes = n(),
              n_hits  = sum(n_hits),
              n_first = sum(n_first),
              .by = genus) |>
    arrange(desc(n_hits))
  
  if (nrow(out) > 0) {
    message("\n", msg)
    print(as.data.frame(head(out, n)))
  }
  invisible(out)
}

flag_genera(pathway_lookup, "osborne_mixed",
            paste("MIXED C3/C4 Poaceae genera present - resolve at species",
                  "level only if the hit counts warrant it:"))

flag_genera(pathway_lookup, "non_poaceae_mixed",
            "MIXED C3/C4 Cyperaceae genera present:")

flag_genera(pathway_lookup, "unmatched_genus",
            paste("Graminoid genera in neither table - check for name",
                  "mismatches or missing non-Poaceae genera:"))

no_name <- pathway_lookup |> filter(pathway_source == "no_scientific_name")
if (nrow(no_name) > 0) {
  message("\nGraminoid codes with no ScientificName (", nrow(no_name),
          " codes, ", sum(no_name$n_hits), " hits):")
  print(as.data.frame(head(arrange(no_name, desc(n_hits)), 20)))
}


# ---- summary and output --------------------------------------------------
pathway_lookup |>
  summarise(n_codes = n(), n_hits = sum(n_hits), .by = pathway_source) |>
  mutate(pct_hits = round(100 * n_hits / sum(n_hits), 2)) |>
  arrange(desc(n_hits)) |>
  as.data.frame() |>
  print()

pathway_lookup |>
  summarise(n_codes = n(), n_hits = sum(n_hits), .by = pathway) |>
  mutate(pct_hits = round(100 * n_hits / sum(n_hits), 2)) |>
  as.data.frame() |>
  print()

write_csv(select(pathway_lookup, code, genus, pathway, pathway_source),
        file.path(ldc_processed_dir, "graminoid_pathway.csv"))
