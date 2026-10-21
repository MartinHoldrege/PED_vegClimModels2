# Download MACAv2-LIVNEH monthly climate projections (pr, tasmax, tasmin)
# for all 20 GCMs from the NKN THREDDS server.
#
# Periods downloaded:
#   historical: 1970-1989, 1990-2005
#   rcp45, rcp85: 2066-2085, 2086-2099
#
# Files are written to
#   paths$large/Data_raw/macaClimateProjections/Data
# which is where 01_summarise_projected_climate_data.R reads them from.
# Re-running skips files that already exist, so an interrupted run can be
# restarted.
# Note: if on HPC run this script as a job using slurm/download_maca.sh

source("config/paths.R") # provides paths
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(xml2)

# Parameters --------------------------------------------------------------

test_run <- TRUE # if TRUE, download only one file (to check URLs work)

vars <- c("pr", "tasmax", "tasmin")
periods <- tribble(
  ~scenario,    ~period,
  "historical", "1970_1989",
  "historical", "1990_2005",
  "rcp45",      "2066_2085",
  "rcp45",      "2086_2099",
  "rcp85",      "2066_2085",
  "rcp85",      "2086_2099"
)

thredds_base <- "https://www.reacchpna.org/thredds"
catalog_dir <- "NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE/macav2livneh"

# Folder names from the macav2livneh catalog (accessed 2026-10-20)
gcms <- c(
  "BNU-ESM", "CCSM4", "CNRM-CM5", "CSIRO-Mk3-6-0", "CanESM2",
  "GFDL-ESM2G", "GFDL-ESM2M", "HadGEM2-CC365", "HadGEM2-ES365",
  "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC-ESM",
  "MIROC-ESM-CHEM", "MIROC5", "MRI-CGCM3", "NorESM1-M",
  "bcc-csm1-1", "bcc-csm1-1-m", "inmcm4"
)

out_dir <- file.path(paths$large, "Data_raw/macaClimateProjections/Data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = 3600) # default 60 s is too short for ~100 MB files

# List files on the server ------------------------------------------------

#' List datasets in a THREDDS catalog folder
#'
#' @param url URL of the folder's catalog.xml.
#' @return Tibble with columns `file` (dataset name) and `url_path`
#'   (path used to build the download URL).
#' @examples
#' \dontrun{
#' url <- paste0(
#'   "https://www.reacchpna.org/thredds/catalog/",
#'   "NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE/macav2livneh/BNU-ESM/catalog.xml"
#' )
#' x <- list_thredds_files(url)
#' head(x)
#' }
list_thredds_files <- function(url) {
  doc <- xml2::read_xml(url)
  xml2::xml_ns_strip(doc)
  ds <- xml2::xml_find_all(doc, "//dataset[@urlPath]")
  tibble(
    file = xml2::xml_attr(ds, "name"),
    url_path = xml2::xml_attr(ds, "urlPath")
  )
}

files <- tibble(gcm_folder = gcms) |>
  mutate(listing = map(
    gcm_folder,
    \(g) list_thredds_files(
      file.path(thredds_base, "catalog", catalog_dir, g, "catalog.xml")
    ),
    .progress = TRUE
  )) |>
  unnest(listing)

# e.g. macav2livneh_pr_BNU-ESM_r1i1p1_rcp85_2066_2085_CONUS_monthly.nc
pattern <- "^macav2livneh_[^_]+_[^_]+_[^_]+_[^_]+_\\d{4}_\\d{4}_CONUS_[a-z]+\\.nc$"

parsed <- files |>
  filter(str_detect(file, pattern)) |>
  separate_wider_regex(
    file,
    patterns = c(
      "macav2livneh_", var = "[^_]+", "_", gcm = "[^_]+", "_",
      member = "[^_]+", "_", scenario = "[^_]+", "_",
      period = "\\d{4}_\\d{4}", "_CONUS_", freq = "[a-z]+", "\\.nc"
    ),
    cols_remove = FALSE
  )

# Select files to download ------------------------------------------------

to_get <- parsed |>
  filter(freq == "monthly", var %in% vars) |>
  semi_join(periods, by = c("scenario", "period")) |>
  mutate(
    url = file.path(thredds_base, "fileServer", url_path),
    dest = file.path(out_dir, file)
  )

# Expect one file per GCM x variable x scenario-period
n_expected <- length(gcms) * length(vars) * nrow(periods)
if (nrow(to_get) != n_expected) {
  warning("Found ", nrow(to_get), " files; expected ", n_expected)
}

if (test_run) {
  to_get <- slice(to_get, 1)
  message("test_run = TRUE: downloading only ", to_get$file)
}

# Download ----------------------------------------------------------------

#' Download one file, writing to a temporary name first so an interrupted
#' download never leaves a partial file under the final name.
#'
#' @param url Download URL.
#' @param dest Destination file path.
#' @return TRUE if the file exists after the attempt, otherwise FALSE.
#' @examples
#' \dontrun{
#' url <- paste0(
#'   "https://www.reacchpna.org/thredds/fileServer/",
#'   "NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE/macav2livneh/BNU-ESM/",
#'   "macav2livneh_tasmin_BNU-ESM_r1i1p1_rcp85_2086_2099_CONUS_monthly.nc"
#' )
#' dest <- file.path(tempdir(), basename(url))
#' download_one(url, dest)
#' terra::nlyr(terra::rast(dest)) # 168 (14 years x 12 months)
#' }
download_one <- function(url, dest) {
  if (file.exists(dest)) return(TRUE)
  tmp <- paste0(dest, ".part")
  ok <- tryCatch(
    download.file(url, tmp, mode = "wb", quiet = TRUE) == 0,
    error = function(e) {
      message("  failed: ", basename(dest), " (", conditionMessage(e), ")")
      FALSE
    }
  )
  if (ok) file.rename(tmp, dest) else unlink(tmp)
  file.exists(dest)
}

message("Downloading ", sum(!file.exists(to_get$dest)), " of ",
        nrow(to_get), " files to ", out_dir)

to_get$downloaded <- map2_lgl(to_get$url, to_get$dest, download_one,
                              .progress = TRUE)

failed <- filter(to_get, !downloaded)
if (nrow(failed) > 0) {
  warning(nrow(failed), " file(s) failed; re-run the script to retry")
}

# Check number of monthly layers ------------------------------------------
# Catches truncated files that downloaded without an error.

check <- to_get |>
  filter(downloaded) |>
  mutate(
    n_expected = (as.integer(str_sub(period, 6, 9)) -
                    as.integer(str_sub(period, 1, 4)) + 1) * 12,
    n_layers = map_int(dest, \(f) terra::nlyr(terra::rast(f)))
  ) |>
  filter(n_layers != n_expected)

if (nrow(check) > 0) {
  warning(nrow(check), " file(s) have the wrong number of layers; see `check`")
} else {
  message("All downloaded files have the expected number of monthly layers")
}