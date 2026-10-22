# Download MACAv2-METDATA monthly climate projections (pr, tasmax, tasmin)
# for all 20 GCMs from the NKN THREDDS server. MACAv2-METDATA is trained on,
# and on the same ~4 km grid as, gridMET.
#
# Files come in 5-year chunks; every chunk overlapping one of the normal
# windows below is downloaded (historical 1971-2000; rcp45 and rcp85
# 2040-2069 and 2070-2099).
#
# Files are written to
#   paths$large/Data_raw/macav2metdata
# which is where 01_summarise_projected_climate_data.R reads them from.
# Re-running skips files that already exist, so an interrupted run can be
# restarted.
# Note: if on HPC run this script as a job using hpc/download_maca.sh

source("config/paths.R") # provides paths
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(xml2)

# Parameters --------------------------------------------------------------

test_run <- FALSE # if TRUE, download only one file (to check URLs work)

vars <- c("pr", "tasmax", "tasmin")
# Normal windows; files are downloaded if they overlap any of these.
windows <- tribble(
  ~scenario,    ~win_start, ~win_end,
  "historical", 1971,       2000,
  "rcp45",      2040,       2069,
  "rcp45",      2070,       2099,
  "rcp85",      2040,       2069,
  "rcp85",      2070,       2099
)

thredds_base <- "https://www.reacchpna.org/thredds"
catalog_dir <- "MACAV2"

# Folder names from the MACAV2 (METDATA) catalog (accessed 2026-10-22)
gcms <- c(
  "BNU-ESM", "CCSM4", "CNRM-CM5", "CSIRO-Mk3-6-0", "CanESM2",
  "GFDL-ESM2G", "GFDL-ESM2M", "HadGEM2-CC365", "HadGEM2-ES365",
  "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC-ESM",
  "MIROC-ESM-CHEM", "MIROC5", "MRI-CGCM3", "NorESM1-M",
  "bcc-csm1-1", "bcc-csm1-1-m", "inmcm4"
)

out_dir <- file.path(paths$large, "Data_raw/macav2metdata")
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
#'   "MACAV2/BNU-ESM/catalog.xml"
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

#' Call a function, retrying after a pause if it errors
#'
#' The server intermittently refuses connections, so network calls are
#' retried rather than failing the whole run.
#'
#' @param f Function to call.
#' @param ... Arguments passed to `f`.
#' @param tries Maximum number of attempts.
#' @param wait Seconds to wait between attempts.
#' @return The value of `f(...)`; errors if all attempts fail.
#' @examples
#' \dontrun{
#' url <- paste0(
#'   "https://www.reacchpna.org/thredds/catalog/",
#'   "MACAV2/BNU-ESM/catalog.xml"
#' )
#' x <- with_retry(list_thredds_files, url, tries = 3, wait = 5)
#' }
with_retry <- function(f, ..., tries = 5, wait = 300) {
  for (i in seq_len(tries)) {
    out <- tryCatch(f(...), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    message("  attempt ", i, " of ", tries, " failed: ", conditionMessage(out))
    if (i < tries) Sys.sleep(wait)
  }
  stop(out)
}

# The catalog listing takes ~5 min, so it is cached and reused on re-runs.
# Delete this file to force a fresh listing.
listing_file <- file.path(out_dir, "thredds_listing.rds")

if (file.exists(listing_file)) {
  files <- readRDS(listing_file)
} else {
  files <- tibble(gcm_folder = gcms) |>
    mutate(listing = map(
      gcm_folder,
      \(g) with_retry(
        list_thredds_files,
        file.path(thredds_base, "catalog", catalog_dir, g, "catalog.xml")
      ),
      .progress = TRUE
    )) |>
    unnest(listing)
  saveRDS(files, listing_file)
}

# e.g. macav2metdata_pr_BNU-ESM_r1i1p1_rcp85_2066_2070_CONUS_monthly.nc
pattern <- "^macav2metdata_[^_]+_[^_]+_[^_]+_[^_]+_\\d{4}_\\d{4}_CONUS_[a-z]+\\.nc$"

parsed <- files |>
  filter(str_detect(file, pattern)) |>
  separate_wider_regex(
    file,
    patterns = c(
      "macav2metdata_", var = "[^_]+", "_", gcm = "[^_]+", "_",
      member = "[^_]+", "_", scenario = "[^_]+", "_",
      period = "\\d{4}_\\d{4}", "_CONUS_", freq = "[a-z]+", "\\.nc"
    ),
    cols_remove = FALSE
  )

# Select files to download ------------------------------------------------

to_get <- parsed |>
  filter(freq == "monthly", var %in% vars) |>
  mutate(file_start = as.integer(str_sub(period, 1, 4)),
         file_end = as.integer(str_sub(period, 6, 9))) |>
  # keep files that overlap any normal window for their scenario
  inner_join(windows, by = "scenario", relationship = "many-to-many") |>
  filter(file_start <= win_end, file_end >= win_start) |>
  distinct(file, .keep_all = TRUE) |>
  select(-win_start, -win_end) |>
  mutate(
    url = file.path(thredds_base, "fileServer", url_path),
    dest = file.path(out_dir, file)
  )

# Every GCM should have the same number of files
n_by_gcm <- count(to_get, gcm)
if (n_distinct(n_by_gcm$n) != 1 || nrow(n_by_gcm) != length(gcms)) {
  warning("Number of files differs among GCMs; see `n_by_gcm`")
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
#' @param tries Maximum number of attempts.
#' @param wait Seconds to wait between attempts.
#' @return TRUE if the file exists after the attempts, otherwise FALSE.
#' @examples
#' \dontrun{
#' url <- paste0(
#'   "https://www.reacchpna.org/thredds/fileServer/",
#'   "MACAV2/BNU-ESM/",
#'   "macav2metdata_tasmin_BNU-ESM_r1i1p1_rcp85_2096_2099_CONUS_monthly.nc"
#' )
#' dest <- file.path(tempdir(), basename(url))
#' download_one(url, dest)
#' terra::nlyr(terra::rast(dest)) # 48 (4 years x 12 months)
#' }
download_one <- function(url, dest, tries = 5, wait = 300) {
  if (file.exists(dest)) return(TRUE)
  Sys.sleep(2) # be gentle with the server
  tmp <- paste0(dest, ".part")
  ok <- tryCatch(
    with_retry(
      \() stopifnot(download.file(url, tmp, mode = "wb", quiet = TRUE) == 0),
      tries = tries, wait = wait
    ),
    error = function(e) {
      message("  giving up on ", basename(dest))
      FALSE
    }
  )
  if (!isFALSE(ok)) file.rename(tmp, dest) else unlink(tmp)
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