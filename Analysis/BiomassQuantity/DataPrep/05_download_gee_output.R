# Martin Holdrege
# Script started April 2026
#
# Purpose: Download GeoTIFFs exported to Google Drive from GEE
# (exported in Analysis/GEE/03_export_tifs.js)
# Only downloads if the Drive version is newer than the local copy.

# dependencies ------------------------------------------------------------

library(googledrive)

source("Functions/init.R")
source("Functions/read_write.R")

# params ------------------------------------------------------------------

drive_folder <- "PED_vegClimModels2"
base_local <- file.path(paths$large, "Data_processed/")

# local destinations per file type
dest <- list(
  lcmap = file.path(base_local, "masks"),
  fire  = file.path(base_local, "masks"),
  rap_cover = file.path(base_local, "CoverData/rap/"),
  rap_biomass = file.path(base_local, "BiomassQuantityData/rap")
)

# create directories if needed
walk(unique(unlist(dest)), dir.create, recursive = FALSE, showWarnings = FALSE)

# patterns and their local destinations
file_specs <- tribble(
  ~pattern,          ~local_dir,
  "^LCMAP",          dest$lcmap,
  "^MTBS",           dest$fire,
  "^RAP_v\\d_cover",   dest$rap_cover,
  "^RAP_v\\d_herb",    dest$rap_biomass,
  "^RAP_v\\d_fracNotForest", dest$rap_cover
)

# list files on Drive -----------------------------------------------------

# drive_auth() # run to setup authentication
drive_files <- drive_ls(path = drive_folder, pattern = "\\.tif$") %>%
  mutate(modifiedTime = map_chr(drive_resource, function(x) x$modifiedTime)) %>%
  # if duplicates, keep only the newest
  group_by(name) %>%
  filter(modifiedTime == max(modifiedTime)) %>%
  ungroup()

# download ----------------------------------------------------------------

# match files to destinations and download
for (i in seq_len(nrow(file_specs))) {
  matched <- drive_files %>%
    filter(str_detect(name, file_specs$pattern[i]))
  
  if (nrow(matched) == 0) {
    message("No files matched pattern: ", file_specs$pattern[i])
    next
  }
  
  for (j in seq_len(nrow(matched))) {
    download_if_newer(matched[j, ], file_specs$local_dir[i])
  }
}
