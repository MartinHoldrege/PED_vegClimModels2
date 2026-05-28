# Script started May 2026
#
# Purpose: Download Jornada Basin LTER seasonal aboveground biomass (AGB)
# data from EDI. This data set provides per-quadrat, per-species, per-season
# AGB estimates (g/m^2) from 1989-ongoing at 15 sites spanning 5 Chihuahuan
# Desert ecosystem types, including three shrubland types (creosotebush-,
# mesquite-, and tarbush-dominated).
#
# Note: this is the raw seasonal AGB package The sites
# are referred to as "NPP study sites" because they are part of the
# long-term NPP study, but the data here are biomass values.
#
# Package: knb-lter-jrn.210011001
# https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-jrn.210011001.4
#
# Saves all data entities to jornada_agb/.
# Run once (or rerun to refresh if a new package revision is released).

# dependencies ------------------------------------------------------------

library(EDIutils)
source("Functions/init.R")

# params ------------------------------------------------------------------
# https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-jrn.210011001.4
scope <- "knb-lter-jrn"
identifier <- 210011001

out_dir <- file.path(paths$large, "Data_raw", "BiomassDataSources", "jornada_agb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# resolve latest package revision -----------------------------------------

revisions <- list_data_package_revisions(scope = scope, identifier = identifier)
latest_rev <- max(as.integer(revisions))
package_id <- paste(scope, identifier, latest_rev, sep = ".")
message("Using package: ", package_id)

# list entities -----------------------------------------------------------

entities <- read_data_entity_names(packageId = package_id)
print(entities)

# download ----------------------------------------------------------------

for (i in seq_len(nrow(entities))) {
  entity_name <- entities$entityName[i]
  entity_id   <- entities$entityId[i]
  
  slug <- entity_name |>
    tolower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
  dest_name <- paste0(slug, ".csv")
  
  dest <- file.path(out_dir, dest_name)
  
  message("Downloading: ", dest_name)
  raw <- read_data_entity(packageId = package_id, entityId = entity_id)
  data <- read_csv(raw)
  write_csv(data, dest)
}

eml <- read_metadata(packageId = package_id)
xml2::write_xml(eml, file.path(out_dir, "eml_metadata.xml"))

# save provenance --------------------------------------------------------

# record exactly which version was downloaded, for reproducibility
writeLines(
  c(paste0("Downloaded: ", Sys.time()),
    paste0("Package ID: ", package_id),
    paste0("Source: https://portal.edirepository.org/nis/metadataviewer?packageid=",
           package_id)),
  con = file.path(out_dir, "provenance.txt")
)

message("Done. Files in: ", out_dir)
