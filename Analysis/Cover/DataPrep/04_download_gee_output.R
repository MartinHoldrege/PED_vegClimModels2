# Download the mask and RAP cover GeoTIFFs from Google Drive.
#
# Exported in:
#   03_export_masks.js         LCMAP_*, MTBS_*
#   01_rap_cover.js            RAP_v3_cover-*
#   03_rap_frac-zero-tree.js   RAP_v3_fracZeroTree_*
#
# Only downloads if the Drive version is newer than the local copy.

library(googledrive)

source("Functions/init.R")
source("Functions/read_write.R")

drive_folder <- "PED_vegClimModels2"
base_local <- file.path(paths$large, "Data_processed")

file_specs <- tribble(
  ~pattern,                 ~local_dir,
  "^(LCMAP|MTBS).*gt.*tif$",      file.path(base_local, "masks"),
  "^RAP_v\\d_cover-.*tif$",       file.path(base_local, "CoverData/rap"),
  "^RAP_v\\d_fracZeroTree.*tif$", file.path(base_local, "CoverData/rap"),
  "^RAP_v\\d_zeroTree_.*cell-years.*csv$", file.path(base_local, "CoverData/rap")
)

walk(file_specs$local_dir, dir.create, recursive = TRUE, showWarnings = FALSE)

# drive_auth() # run once to set up authentication
drive_files <- drive_ls(path = drive_folder) |>
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