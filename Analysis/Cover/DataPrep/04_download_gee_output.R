# Download the mask and RAP cover GeoTIFFs from Google Drive.
#
# Exported in:
#   04_export_masks.js   LCMAP_*, MTBS_*
#   05_rap_cover.js      RAP_v3_cover-*
#
# Only downloads if the Drive version is newer than the local copy.

library(googledrive)

source("Functions/init.R")
source("Functions/read_write.R")

drive_folder <- "PED_vegClimModels2"
base_local <- file.path(paths$large, "Data_processed")

file_specs <- tribble(
  ~pattern,             ~local_dir,
  "^(LCMAP|MTBS).*gt",      file.path(base_local, "masks"),
  "^RAP_v\\d_cover.*thin",    file.path(base_local, "CoverData/rap")
)

walk(file_specs$local_dir, dir.create, recursive = TRUE, showWarnings = FALSE)

# drive_auth() # run once to set up authentication
drive_files <- drive_ls(path = drive_folder, pattern = "\\.tif$") |>
  mutate(modifiedTime = map_chr(drive_resource, \(x) x$modifiedTime)) |>
  slice_max(modifiedTime, n = 1, by = name)

for (i in seq_len(nrow(file_specs))) {
  matched <- filter(drive_files, str_detect(name, file_specs$pattern[i]))
  
  if (nrow(matched) == 0) {
    message("no files matched: ", file_specs$pattern[i])
    next
  }
  
  walk(seq_len(nrow(matched)),
       \(j) download_if_newer(matched[j, ], file_specs$local_dir[i]))
}