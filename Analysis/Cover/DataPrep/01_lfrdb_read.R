# 01_lfrdb_read.R
#
# Read the six regional LANDFIRE Reference Database (LFRDB) Access files and
# combine them into two cached tables.
#
# LFRDB is plot-level field data compiled by LANDFIRE from ~300 source
# programs, NOT a LANDFIRE raster product. LF 2.0.0 (2016 Remap) is the most
# recent release. 
#
# Tables used (see LFRemap_PUBLIC_LFRDB_DataDictionary.pdf):
#   dtStands  lifeform cover and height within the sampled unit
#   dtSpecies species-level cover within the sampled unit
#   dtPoints  location (Lat, Long, LocMeth)
#   dtVisits  year, source program, protocol
#
# Requires the mdbr package, which needs mdbtools installed on the system.
#
# Outputs:
#   lfrdb_stands.rds  - one row per EventID, lifeform cover + location + visit
#   lfrdb_species.rds - one row per EventID x species
#
# September, 2026


# dependencies ------------------------------------------------------------


source("Functions/init.R")
source("Functions/plant_traits.R")

stopifnot(requireNamespace("mdbr", quietly = TRUE))

# params ------------------------------------------------------------------


lfrdb_dir <- file.path(paths$large, "Data_raw/LANDFIRE_LFRDB")
lfrdb_processed_dir <- file.path(paths$large, "Data_processed/LANDFIRE_LFRDB")
dir.create(lfrdb_processed_dir, showWarnings = FALSE, recursive = TRUE)
rerun <- FALSE
regions <- c("NC", "NE", "NW", "SC", "SE", "SW")


# functions ---------------------------------------------------------------


#' Path to one region's Access database
#' @param region Two-letter region code.
lfrdb_path <- function(region) {
  file.path(lfrdb_dir,
            sprintf("%s_Public_LFRDB_LF2.0.0", region),
            sprintf("%s_Public_LFRDB_LF2.0.0.accdb", region))
}

missing_db <- regions[!file.exists(map_chr(regions, lfrdb_path))]
if (length(missing_db) > 0) {
  stop("missing LFRDB file(s) for region(s): ",
       paste(missing_db, collapse = ", "))
}

#' Read one table from one region, tagging the region
#' @param region Two-letter region code.
#' @param table Table name, e.g. "dtStands".
read_lfrdb <- function(region, table) {
  message("reading ", region, " / ", table)
  mdbr::read_mdb(lfrdb_path(region), table = table) |>
    as_tibble() |>
    mutate(region = region)
}

#' Read and combine one table across all six regions
#'
#' Regions are separate databases with no shared EventIDs, so a plain bind is
#' correct. The check below confirms that.
#'
#' @param table Table name.
read_lfrdb_all <- function(table) {
  out <- map(regions, read_lfrdb, table = table) |> list_rbind()
  stopifnot(nrow(out) > 0)
  out
}


# read in data ------------------------------------------------------------


# ---- stands: lifeform cover, joined to location and visit ----------------
stands <- fetch_cached(
  file.path(lfrdb_processed_dir, "lfrdb_stands.rds"),
  \() {
    dt_stands <- read_lfrdb_all("dtStands")
    dt_points <- read_lfrdb_all("dtPoints") |> select(-region)
    dt_visits <- read_lfrdb_all("dtVisits") |> select(-region)
    
    # EventID is assigned by LANDFIRE and must be unique across regions;
    # if it is not, the joins below will fan out silently.
    check_unique(dt_stands, EventID, label = "dtStands")
    check_unique(dt_points, EventID, label = "dtPoints")
    check_unique(dt_visits, EventID, label = "dtVisits")
    
    out <- dt_stands |>
      left_join(dt_points, by = "EventID") |>
      left_join(dt_visits, by = "EventID")
    
    stopifnot(nrow(out) == nrow(dt_stands))
    out
  },
  rerun = rerun
)

# ---- species: one row per EventID x species ------------------------------
species <- fetch_cached(
  file.path(lfrdb_processed_dir, "lfrdb_species.rds"),
  \() read_lfrdb_all("dtSpecies"),
  rerun = rerun
)

# ---- structural checks ---------------------------------------------------
check_cols(stands,
           c("EventID", "LFTreeCov", "LFTreeCovAdj", "LFConiferTreeCov",
             "LFShrubCov", "LFShrubCovAdj", "LFHerbCov", "LFHerbCovAdj",
             "SourceTreeCov", "SourceShrubCov", "SourceHerbCov",
             "Lat", "Long", "LocMeth", "YYYY", "SourceID", "Protocol"),
           "stands")

check_cols(species,
           c("EventID", "Item", "SciName", "Lifeform", "Duration",
             "LFAbsCov", "LFRelCov"),
           "species")

# LFAbsCov must be numeric. dtExotics stores categorical infestation levels
# ("P", "L", "M", "H") in a similarly named field; if dtSpecies ever does the
# same, every downstream sum is wrong.
if (!is.numeric(species$LFAbsCov)) {
  stop("species$LFAbsCov is ", class(species$LFAbsCov)[1],
       ", expected numeric. Check for categorical cover codes.")
}

# a species should appear once per plot; repeats would double-count on sum
check_unique(species, EventID, Item, label = "species (EventID x Item)")

# ---- audit ---------------------------------------------------------------
message("plots with species rows: ", n_distinct(species$EventID),
        " (", round(100 * n_distinct(species$EventID) / nrow(stands), 1), "%)")

# LocMeth: G = GPS, M = digitised in office, X = unknown. Plots with unknown
# locations cannot be joined to gridded covariates.
message("\nLocMeth (G = GPS, M = office, X = unknown):")
count(stands, LocMeth) |> mutate(pct = round(100 * n / sum(n), 1)) |>
  as.data.frame() |> print()

# Lifeform codes present, for the group mapping in 06_. An unexpected code
# there means 06_'s lifeform_class table needs extending.
message("\nLifeform codes (F forb, G graminoid, H herb, N nonvascular, ",
        "S shrub, T tree, V vine):")
count(species, Lifeform, sort = TRUE) |> as.data.frame() |> print()

message("\nsampling years:")
count(stands, YYYY) |> arrange(YYYY) |> as.data.frame() |> print()

# ---- provenance: does LFRDB overlap the LDC plots? -----------------------
# Text search of the source metadata found no AIM / NRI / LMF programs, but
# a contribution could be relabelled. Nearest-neighbour distance is the
# stronger test: a real overlap puts LFRDB plots within metres of LDC plots.
# file created in "Analysis/Cover/DataPrep/04_ldc_process.R", so 
# this check can only be done after 04_ldc_process.R has run (and it's optional)
# --it shows that the plots are only ~48 plots within 50m, so essentially
# no overlap)

ldc_cover <- file.path(paths$large,
                       "Data_processed/LandscapeDataCommonsDat",
                       "cover_by_plot_year.csv")

if (file.exists(ldc_cover)) {
  ldc_xy <- read_csv(ldc_cover, show_col_types = FALSE) |>
    distinct(Longitude_NAD83, Latitude_NAD83) |>
    sf::st_as_sf(coords = c("Longitude_NAD83", "Latitude_NAD83"), crs = 4269) |>
    sf::st_transform(5070)
  
  lf_xy <- stands |>
    filter(!is.na(Long), !is.na(Lat)) |>
    distinct(SourceID, Long, Lat)
  
  lf_sf  <- sf::st_as_sf(lf_xy, coords = c("Long", "Lat"), crs = 4326) |>
    sf::st_transform(5070)
  
  idx <- sf::st_nearest_feature(lf_sf, ldc_xy)
  lf_xy$dist_m <- as.numeric(sf::st_distance(lf_sf, ldc_xy[idx, ], by_element = TRUE))
  
  message("\nLFRDB plots near an LDC plot (possible shared provenance):")
  lf_xy |> summarise(n = n(),
                     within_50m  = sum(dist_m < 50),
                     within_500m = sum(dist_m < 500)) |>
    as.data.frame() |> print()
  
  near <- lf_xy |> filter(dist_m < 50) |> count(SourceID, sort = TRUE)
  if (nrow(near) > 0) {
    message("SourceIDs with co-located plots - check whether these are the ",
            "same network under another label:")
    print(as.data.frame(head(near, 20)))
  }
} else {
  message("\nskipping provenance check (needs cover_by_plot_year.csv and nngeo)")
}