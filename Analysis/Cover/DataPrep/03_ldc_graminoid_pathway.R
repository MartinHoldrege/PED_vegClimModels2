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
# relies heavily on Functions/data/plant_traits.R
source('Functions/data/plant_traits.R')

ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")
osborne_xls <- file.path(paths$large, "Data_raw/plant_spp_info",
                         "Osborne_etal_2014_newphytologist_tables3_c4-genera.xls")

code_class <- readRDS(file.path(ldc_processed_dir, "code_class.rds")) |> 
  as_tibble()
stopifnot("ScientificName" %in% names(code_class))


# ---- Osborne genus-level pathways (Poaceae) ------------------------------
osborne <- read_osborne(osborne_xls)

# ---- genus cleaning ------------------------------------------------------
# ScientificName is not consistently formatted. 
non_poaceae <- non_poaceae_pathway # from plant_traits.R
known_genera <- unique(c(osborne$genus, non_poaceae$genus))

# ---- graminoid codes in the data -----------------------------------------
graminoids_raw <- code_class |>
  filter(class == "graminoid") |>
  mutate(genus = clean_genus(scientific_name = ScientificName, 
                             genera = known_genera)) |>
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
  left_join(assign_pathway(graminoids$genus, osborne), by = "genus")

# ---- flags ---------------------------------------------------------------
# mixed genera are unresolvable at genus level and fall to grass_unknown
pathway_lookup |>
  filter(pathway == "mixed") |>
  summarise(n_codes = n(), n_hits = sum(n_hits), .by = c(genus, pathway_source)) |>
  arrange(desc(n_hits)) |> as.data.frame() |> print()

flag_genera(pathway_lookup, "c3_by_default",
            "Genera assigned C3 by elimination - scan for missed C4 genera:",
            n = 50)


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
