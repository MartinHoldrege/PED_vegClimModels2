# Script started May 2026
#
# Purpose: Download Wright (2015) sagebrush fuels data from the Forest
# Service Research Data Archive. 
#
# Package: RDS-2015-0016
# https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0016
#
# Saves the zip and its extracted contents to wright_sagebrush_fuels/.

# dependencies ------------------------------------------------------------

source("Functions/init.R")

# params ------------------------------------------------------------------

# https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0016
rds_id <- "RDS-2015-0016"
zip_url <- paste0("https://www.fs.usda.gov/rds/archive/products/",
                  rds_id, "/", rds_id, ".zip")

out_dir <- file.path(paths$large, "Data_raw", "BiomassDataSources",
                     "wright_sagebrush_fuels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zip_path <- file.path(out_dir, paste0(rds_id, ".zip"))

# download ----------------------------------------------------------------

if (file.exists(zip_path)) {
  message("Skipping (exists): ", basename(zip_path))
} else {
  message("Downloading: ", zip_url)
  download.file(zip_url, zip_path, mode = "wb")
}

# extract -----------------------------------------------------------------

message("Extracting to: ", out_dir)
unzip(zip_path, exdir = out_dir, overwrite = TRUE)

# list what we got --------------------------------------------------------

extracted <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
message("Files in package:")
print(extracted)

# save provenance ---------------------------------------------------------

writeLines(
  c(paste0("Downloaded: ", Sys.time()),
    paste0("Package ID: ", rds_id),
    paste0("Source: ", zip_url)),
  con = file.path(out_dir, "provenance.txt")
)

message("Done. Files in: ", out_dir)