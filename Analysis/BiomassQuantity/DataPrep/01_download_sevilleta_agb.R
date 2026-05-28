# Script started May 2026
#
# Purpose: Download Sevilleta LTER core-site seasonal biomass and NPP data
# from EDI. 
#
# Note: this package contains both seasonal AGB and derived seasonal/annual
# NPP entities. The AGB-per-quadrat file is what's of primary interest here;
# the NPP entities are derived products.
#
# Package: knb-lter-sev.182
# https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sev.182
#
# Saves all data entities to sevilleta_agb/.
# Run once (or rerun to refresh if a new package revision is released).

# dependencies ------------------------------------------------------------

library(EDIutils)
source("Functions/init.R")

# params ------------------------------------------------------------------

scope <- "knb-lter-sev"
identifier <- 182

out_dir <- file.path(paths$large, "Data_raw", "BiomassDataSources", "sevilleta_agb")
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
  
  dest_name <- entity_name
  stopifnot(str_detect(entity_name, '.csv$'))
  
  dest <- file.path(out_dir, dest_name)
  
  message("Downloading: ", dest_name)
  raw <- read_data_entity(packageId = package_id, entityId = entity_id)
  data <- read_csv(raw)
  write_csv(data, dest)
}

# save EML metadata -------------------------------------------------------

# pulls the full EML XML for the package, which often contains site-level
# geographic coverage (lat/lon per site) not present in the data files
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