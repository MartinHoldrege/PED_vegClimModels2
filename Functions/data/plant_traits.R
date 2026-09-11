# Functions/plant_traits.R
#
# Shared plant-trait lookups and helpers, used by both the LDC and LFRDB
# pipelines. Everything here operates on a scientific name or a genus, so it
# is independent of which survey the data came from.
#
# Sources:
#   Osborne et al. (2014) New Phytologist 204:441-446, Table S3 - Poaceae
#     genus-level C3/C4
#   Bruhl & Wilson (2007) Aliso 23:99-148, Table 1 - Cyperaceae C3/C4
#
# September, 2026

# ---- genus cleaning ------------------------------------------------------
# Scientific-name fields are not consistently formatted in either source.
# Three failure modes matter, because an unrecognised genus is silently
# defaulted to C3:
#   1. authority text leaks in       ("Juncus L.", "Vulpia C.C.")
#   2. epithet runs into the genus   ("Bromustectorum", "Poasecunda")
#   3. placeholders instead of names ("Unknown", "Generic", "#N/A")

# Strings carrying no botanical information; set to NA rather than defaulted.
bad_genus <- c("Unknown", "Generic", "Perennial", "Check", "#N/A", "NA")

# Common names appearing in the scientific-name field.
common_names <- tibble::tribble(
  ~genus,   ~genus_fixed,
  "Fescue", "Festuca"
)

# Conifer genera occurring in North America, native and planted. Juniperus,
# Calocedrus, Thuja and Hesperocyparis are scale-leaved and are defined as
# needle-leaved here.
.needle_genera <- c(
  "Abies", "Calocedrus",  "Cedrus", "Chamaecyparis", "Cupressus",
  "Hesperocyparis", "Juniperus", "Larix", "Picea", "Pinus", "Pseudotsuga",
  "Sequoia", "Sequoiadendron", "Taxodium", "Taxus", "Thuja", "Torreya",
  "Tsuga", "Callitropsis", "Cryptomeria", "Metasequoia", "Platycladus",
  "Podocarpus", "Araucaria", 'Pine', 'Pseudolarix', 'Afrocarpus'
)

# Cyperaceae C4 and C3/C4 genera, Bruhl & Wilson (2007) Table 1, in table
# order, plus C4 genera absent from Osborne (synonyms, segregates,
# misspellings). "mixed" means the genus contains both pathways and cannot be
# assigned at genus level. All other Cyperaceae genera in Table 1 are C3, so a
# genus matching neither source is C3 by elimination - Juncaceae included.
non_poaceae_pathway <- tibble::tribble(
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
  "Volkiella",        "C4",        NA_character_,
  "Dasyochloa",       "C4",        "= Erioneuron",
  "Pleuraphis",       "C4",        "= Hilaria",
  "Monroa",           "C4",        "= Munroa",
  "Vilfa",            "C4",        "= Sporobolus",
  "Chaetochloa",      "C4",        "= Setaria",
  "Erianthus",        "C4",        "= Saccharum",
  "Sorgastrum",       "C4",        "misspelling of Sorghastrum",
  "Hopia",            "C4",        "segregate of Panicum",
  "Hypogynium",       "C4",        "Andropogoneae"
)
stopifnot(!any(duplicated(non_poaceae_pathway$genus)))
stopifnot(all(non_poaceae_pathway$pathway_np %in% c("C3", "C4", "mixed")))


# ---- small utilities -----------------------------------------------------

#' Stop if a data frame has duplicate values in the given columns
#'
#' @param x A data frame.
#' @param ... Columns expected to be unique together (tidyselect).
#' @param label Name used in the error message.
check_unique <- function(x, ..., label = deparse(substitute(x))) {
  n_dup <- x |> count(...) |> filter(n > 1) |> nrow()
  if (n_dup > 0) {
    stop(label, ": ", n_dup, " duplicated key(s). ",
         "Deduplicate before joining, or the join will fan out.")
  }
  invisible(x)
}

#' Stop if expected columns are missing
#'
#' @param x A data frame.
#' @param cols Character vector of required column names.
#' @param label Name used in the error message.
check_cols <- function(x, cols, label = deparse(substitute(x))) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    stop(label, ": missing expected column(s): ",
         paste(missing, collapse = ", "))
  }
  invisible(x)
}

#' Add any missing columns as zero, then order them
#'
#' Guarantees a stable output schema when a group is absent from the data, so
#' downstream code can rely on the column set regardless of the input.
#'
#' @param x A data frame.
#' @param cols Character vector of columns that must exist.
#' @param fill Value used for columns that have to be created.
ensure_cols <- function(x, cols, fill = 0) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    message("adding ", length(missing), " absent column(s) as ", fill, ": ",
            paste(missing, collapse = ", "))
    x[missing] <- fill
  }
  x
}

#' Cache a computed object to disk, recomputing only if absent
#'
#' Writes to a `.partial` file first so an interrupted write leaves no
#' truncated cache.
#'
#' @param path Path to the .rds cache file.
#' @param fn Zero-argument function returning the object to cache.
fetch_cached <- function(path, fn, rerun = FALSE) {
  if (file.exists(path) & !rerun) return(readRDS(path))
  
  message("computing: ", basename(path))
  x <- fn()
  stopifnot(is.data.frame(x), nrow(x) > 0)
  
  tmp <- paste0(path, ".partial")
  saveRDS(x, tmp, compress = "xz")
  file.rename(tmp, path)
  x
}


# ---- trait sources -------------------------------------------------------

#' Read the Osborne et al. (2014) genus-level C3/C4 table
#'
#' Header sits on row 8; rows above are the caption and explanatory notes.
#' Aristida is recorded as mixed by Osborne on the strength of Old World
#' species, but all North American species are C4 and it is overridden here.
#'
#' @param path Path to the Osborne Table S3 .xls file.
#' @return Tibble with columns genus, pathway ("C3"/"C4"/"mixed"), pathway_raw.
read_osborne <- function(path) {
  stopifnot(file.exists(path))
  
  out <- readxl::read_excel(path, sheet = "Table S1", skip = 7) |>
    rlang::set_names(c("genus_authority", "pathway_raw", "species_sampled",
                       "references", "notes")) |>
    filter(!is.na(genus_authority)) |>
    mutate(genus = word(genus_authority, 1),
           pathway = case_when(str_trim(pathway_raw) == "C3" ~ "C3",
                               str_trim(pathway_raw) == "C4" ~ "C4",
                               .default = "mixed")) |>
    select(genus, pathway, pathway_raw)
  
  # 708 genera in Table S3; a different count means the file or the header
  # offset has changed and the parse is wrong.
  if (nrow(out) != 708) {
    stop("read_osborne: expected 708 genera, got ", nrow(out),
         ". Check the sheet name and `skip` against the file.")
  }
  stopifnot(!any(duplicated(out$genus)))
  
  out$pathway[out$genus == "Aristida"] <- "C4" #in north america common sp are all c4
  out
}


# ---- genus extraction ----------------------------------------------------

#' Split a genus name that has run together with its epithet
#'
#' Returns the longest known genus the string starts with, or the string
#' unchanged if none matches. Fails where a genuine genus absent from
#' `genera` happens to start with a known genus name ("Poaceae" -> "Poa"), so
#' inspect what it changes (see clean_genus, which reports this).
#'
#' @param x Character vector of candidate genus strings.
#' @param genera Character vector of known genus names.
fix_runtogether <- function(x, genera) {
  stopifnot(is.character(genera), length(genera) > 0)
  
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
  
  fix_runtogether(x, known_genera) # fails for 'poaceae'
}

#' Extract a clean genus from a scientific-name string
#'
#' Takes the leading alphabetic run (dropping authorities and punctuation),
#' maps common names, nulls placeholders, and splits run-together names.
#'
#' @param scientific_name Character vector of scientific names.
#' @param genera Known genus names, used by fix_runtogether().
#' @param report If TRUE, message the distinct strings that fix_runtogether
#'   altered, so a bad split can be spotted.
#' @return Character vector of genera, NA where nothing usable was found.
clean_genus <- function(scientific_name, genera, report = TRUE) {
  raw <- word(scientific_name, 1) |> str_remove("[^A-Za-z].*$")
  raw <- if_else(raw %in% bad_genus | raw == "" | is.na(raw),
                 NA_character_, raw)
  
  idx <- match(raw, common_names$genus)
  mapped <- if_else(is.na(idx), raw, common_names$genus_fixed[idx])
  
  out <- fix_runtogether(mapped, genera)
  
  if (report) {
    changed <- tibble(from = mapped, to = out) |>
      filter(!is.na(from), from != to) |>
      distinct() |>
      arrange(from)
    if (nrow(changed) > 0) {
      message("clean_genus: split ", nrow(changed),
              " run-together name(s); check for bad splits:")
      print(as.data.frame(changed))
    }
  }
  
  out
}


# ---- trait assignment ----------------------------------------------------

#' Assign C3/C4 photosynthetic pathway from genus
#'
#' Osborne is complete for Poaceae and Bruhl & Wilson lists every C4 and
#' C3/C4 Cyperaceae genus, so a genus in neither is C3 by elimination. That
#' also absorbs post-Osborne segregate genera, which in North America are
#' Pooideae and therefore C3.
#'
#' @param genus Character vector of cleaned genus names.
#' @param osborne Output of read_osborne().
#' @return Tibble, one row per input, with columns pathway
#'   ("C3"/"C4"/"mixed"/NA) and pathway_source ("osborne"/"bruhl_wilson"/
#'   "no_scientific_name"/"c3_by_default").
assign_pathway <- function(genus, osborne) {
  check_cols(osborne, c("genus", "pathway", "pathway_raw"), "osborne")
  
  out <- tibble(genus = unique(genus)) |>
    left_join(osborne, by = "genus") |>
    left_join(non_poaceae_pathway, by = "genus") |>
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
        .default           = "C3"
      )
    ) |>
    select(genus, pathway, pathway_source)
  
  # a genus in both sources would silently take the Osborne value; the two
  # tables are meant to be disjoint (Poaceae vs Cyperaceae)
  overlap <- intersect(osborne$genus, non_poaceae_pathway$genus)
  if (length(overlap) > 0) {
    stop("assign_pathway: genus in both trait tables: ",
         paste(overlap, collapse = ", "))
  }
  
  out
}

#' Assign needle-leaved vs broad-leaved from genus
#'
#' Conifer genera are a short closed list in North America, so a genus not on
#' it is broadleaf by elimination. Palms and other monocot "trees" fall to
#' broad: fronds are broad laminae and functionally closer to broadleaf
#' evergreen than to conifers.
#'
#' @param genus Character vector of cleaned genus names.
#' @return Tibble, one row per input, with columns leaf_type and
#'   leaf_type_source.
assign_leaf_type <- function(genus) {
  out <- tibble(genus = unique(genus)) |>
    mutate(
      leaf_type = case_when(
        is.na(genus)             ~ NA_character_,
        genus %in% .needle_genera ~ "needle",
        .default                 = "broad"
      ),
      leaf_type_source = case_when(
        is.na(genus)             ~ "no_scientific_name",
        genus %in% .needle_genera ~ "conifer_list",
        .default                 = "broad_by_default"
      )
    ) |>
    select(genus, leaf_type, leaf_type_source)

  out
}


# ---- reporting -----------------------------------------------------------

#' Print a flagged subset of a trait lookup, by genus and weight
#'
#' @param x A lookup carrying `genus`, a source column, and a weight column.
#' @param src Value of `source_col` to report on.
#' @param msg Message printed above the table.
#' @param source_col Name of the source column.
#' @param weight_col Name of the column to sum and sort by.
#' @param n Maximum genera to print.
flag_genera <- function(x, src, msg, source_col = "pathway_source",
                        weight_col = "n_hits", n = 30) {
  check_cols(x, c("genus", source_col, weight_col), "flag_genera input")
  
  out <- x |>
    filter(.data[[source_col]] == src) |>
    summarise(n_codes = n(),
              weight  = sum(.data[[weight_col]], na.rm = TRUE),
              .by = genus) |>
    arrange(desc(weight))
  
  if (nrow(out) > 0) {
    message("\n", msg)
    print(as.data.frame(head(out, n)))
  }
  invisible(out)
}