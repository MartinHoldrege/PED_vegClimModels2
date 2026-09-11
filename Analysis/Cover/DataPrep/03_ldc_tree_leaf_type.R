# 03_ldc_tree_leaf_type.R
#
# Assign needle-leaved vs broad-leaved to tree codes.
#
# Conifers (Pinophyta) are needle-leaved. Juniperus, Calocedrus, Thuja and
# Hesperocyparis are scale-leaved, and defined 
# as needle-leaved here
#
#
# Everything else with a genus is broad-leaved by elimination. Codes with no ScientificName cannot be assigned and are NA.
#
# Input:  code_class.rds (with ScientificName)
# Output: tree_leaf_type.csv (code, genus, leaf_type, leaf_type_source)
#
# September, 2026

source("Functions/init.R")
source("Functions/data/plant_traits.R")

ldc_processed_dir <- file.path(paths$large, "Data_processed/LandscapeDataCommonsDat")
code_class <- readRDS(file.path(ldc_processed_dir, "code_class.rds"))


trees <- code_class |>
  filter(class == "tree") |>
  mutate(
    genus = word(ScientificName, 1) |> str_remove("[^A-Za-z].*$"),
    genus = if_else(genus %in% bad_genus | genus == "", 
                    NA_character_, genus)
  ) |>
  summarise(n_hits  = sum(n_hits,  na.rm = TRUE),
            n_first = sum(n_first, na.rm = TRUE),
            .by = c(code, ScientificName, genus))

leaf_lookup <- trees |>
  left_join(assign_leaf_type(trees$genus), by = "genus")

# ---- audit ---------------------------------------------------------------
leaf_lookup |>
  summarise(n_codes = n(), n_hits = sum(n_hits), .by = leaf_type_source) |>
  mutate(pct = round(100 * n_hits / sum(n_hits), 2)) |>
  arrange(desc(n_hits)) |> as.data.frame() |> print()

# genera defaulted to broadleaf - scan for any conifer missed above
leaf_lookup |>
  filter(leaf_type_source == "broad_by_default") |>
  summarise(n_hits = sum(n_hits), .by = genus) |>
  arrange(desc(n_hits)) |> as.data.frame() |> print()

write_csv(select(leaf_lookup, code, genus, leaf_type, leaf_type_source),
        file.path(ldc_processed_dir, "tree_leaf_type.csv"))
