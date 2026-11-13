# 05_write_combined_all.R
#
# Writes all layers of one version of cover_specs$combined$all
# (read_combined_all()) to a single multi-layer raster, for sharing.
#
# Output (path from combined_path()):
#   Data_processed/CoverData/Predictions/combined_all_<vmc>_<scenario>.tif
#
# November, 2026

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

vmc <- "m21.5" # opt$vmc
scenario <- "current"

r <- read_combined_all(vmc, scenario)
out_file <- combined_path("all", vmc, "prediction", scenario)
terra::writeRaster(r, out_file, overwrite = TRUE, datatype = "FLT4S")
message("Wrote ", out_file, "\n  layers: ", paste(names(r), collapse = ", "))

if(FALSE) {
  # to move a copy locally for sharing
  local <- combined_path("all", vmc, "prediction", scenario, root = '.')
  dir.create(dirname(local))
  file.copy(out_file, local)
}
