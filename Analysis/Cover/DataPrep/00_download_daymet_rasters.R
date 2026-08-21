#///////////////////////////////////////////////////////////////////////////
# Download Daymet V4R1 monthly + annual climate-summary COGs (CONUS "na")
# via the Earthdata CMR granule API, for use in
# GettingDaymetDataForModelPredictions.R
#
# The V4R1 files live in ORNL's cloud archive, so we query CMR for the real
# download URLs rather than constructing THREDDS paths (which returned 404s).
#
#   Monthly summaries: concept C2532007210-ORNL_CLOUD  (dataset 2131)
#   Annual  summaries: concept <SET BELOW>             (dataset 2130)
#
# AUTH: requires a free NASA Earthdata login. Your ~/.netrc (machine
# urs.earthdata.nasa.gov ...) plus ~/.urs_cookies is used by curl. If cloud
# links need a bearer token instead, the test stage will reveal an auth error
# and we adjust then.
#
# WORKFLOW:
#   1. Run with `stage <- "inspect"` first. It queries CMR, prints a few real
#      URLs + filenames, and tries ONE download. Confirm URLs/filenames/auth.
#   2. Then set `stage <- "download"` for the full 1980-2023 pull.

# July 2, 2026
#///////////////////////////////////////////////////////////////////////////

library(stringr)
library(purrr)
library(dplyr)
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

# --- settings ---------------------------------------------------------------

stage <- "download" # "inspect"   # "inspect" (verify first) or "download" (full run)

years <- 1980:2023

# CMR concept IDs (cloud archive). Monthly is confirmed from the dataset page.
# Annual (2130) concept id is NOT yet confirmed -- the inspect stage will look
# it up by short_name if this is left as NA.
concept_monthly <- "C2532007210-ORNL_CLOUD"
concept_annual <- "C2531982907-ORNL_CLOUD"
short_annual    <- "Daymet_Annual_V4R1" # used only if concept_annual is NA

# Destination directories (must match the processing script's paths)
monthly_dir <- "./Data_raw/dayMet/rawMonthlyData/orders/70e0da02b9d2d6e8faa8c97d211f3546/Daymet_Monthly_V4R1/data"
annual_dir  <- "./Data_raw/dayMet/yearly"

# Variable tokens the PROCESSING script's regex expects in the saved filenames
monthly_vars <- c("prcp_monttl", "tmax_monavg", "tmin_monavg")
annual_vars  <- c("prcp_annttl", "tmax_annavg", "tmin_annavg")

cmr_base <- "https://cmr.earthdata.nasa.gov/search/granules.json"

# --- CMR query --------------------------------------------------------------

#' Fetch all granule download URLs for a CMR collection
#'
#' @param concept_id Character CMR collection concept id, or NA to use short_name.
#' @param short_name Character collection short_name (used only if concept_id is NA).
#' @return Character vector of data-download URLs (http/https, ending .tif).
get_cmr_urls <- function(concept_id = NA_character_, short_name = NA_character_) {
  q <- if (!is.na(concept_id)) {
    str_glue("{cmr_base}?concept_id={concept_id}&page_size=2000")
  } else {
    str_glue("{cmr_base}?short_name={short_name}&page_size=2000")
  }
  res <- jsonlite::fromJSON(q, simplifyVector = FALSE)
  entries <- res$feed$entry
  if (length(entries) == 0) return(character(0))
  
  # each entry has a $links list; keep data-download links to .tif
  urls <- entries %>%
    map(~ map_chr(.x$links, "href", .default = NA_character_)) %>%
    unlist(use.names = FALSE)
  
  urls <- urls[!is.na(urls)]
  urls[str_detect(urls, "\\.tif$") & str_starts(urls, "http")]
}

#' Pick the URL for one variable/year from a pool of CMR URLs
#'
#' Matches on the variable token and the 4-digit year appearing in the URL.
#' @return Single URL string, or NA if not found / ambiguous.
pick_url <- function(urls, var, year) {
  hit <- urls[str_detect(urls, fixed(var)) & str_detect(urls, as.character(year))]
  # prefer the North America ("na") file if region tokens are present
  na_hit <- hit[str_detect(hit, "_na_")]
  if (length(na_hit) >= 1) hit <- na_hit
  if (length(hit) == 1) hit else NA_character_
}

#' Download one file to the processing-script's expected name
#' @return Logical success.
download_one <- function(url, dest_dir, var, year) {
  if (is.na(url)) { message("no URL: ", var, " ", year); return(FALSE) }
  fname <- str_glue("daymet_v4_{var}_na_{year}.tif")  # name the processor expects
  dest  <- file.path(dest_dir, fname)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("skip (exists): ", fname); return(TRUE)
  }
  ok <- tryCatch({
    download.file(url, dest, method = "curl", quiet = TRUE,
                  extra = "-n -L -c ~/.urs_cookies -b ~/.urs_cookies --fail --silent --show-error")
    file.exists(dest) && file.info(dest)$size > 0
  }, error = function(e) {
    message("FAILED: ", fname, " -- ", conditionMessage(e))
    if (file.exists(dest)) file.remove(dest); FALSE
  })
  if (ok) message("ok: ", fname)
  ok
}

# --- run --------------------------------------------------------------------

dir.create(monthly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annual_dir,  recursive = TRUE, showWarnings = FALSE)

message("Querying CMR for monthly granule URLs...")
urls_monthly <- get_cmr_urls(concept_id = concept_monthly)
message("  found ", length(urls_monthly), " monthly .tif URLs")

message("Querying CMR for annual granule URLs...")
urls_annual <- if (!is.na(concept_annual)) {
  get_cmr_urls(concept_id = concept_annual)
} else {
  get_cmr_urls(short_name = short_annual)
}
message("  found ", length(urls_annual), " annual .tif URLs")

if (stage == "inspect") {
  
  message("\n--- sample monthly URLs ---")
  print(utils::head(urls_monthly, 5))
  message("\n--- sample annual URLs ---")
  print(utils::head(urls_annual, 5))
  
  message("\nAttempting ONE test download (monthly prcp 2020)...")
  test_url <- pick_url(urls_monthly, "prcp_monttl", 2020)
  message("matched URL: ", ifelse(is.na(test_url), "<none>", test_url))
  download_one(test_url, monthly_dir, "prcp_monttl", 2020)
  
  message("\nInspect the URLs and test result above. If the test file downloaded,")
  message("set stage <- \"download\" and rerun for the full 1980-2023 pull.")
  
} else if (stage == "download") {
  
  plan <- bind_rows(
    expand.grid(var = monthly_vars, year = years, stringsAsFactors = FALSE) %>%
      mutate(set = "monthly"),
    expand.grid(var = annual_vars, year = years, stringsAsFactors = FALSE) %>%
      mutate(set = "annual")
  )
  
  results <- pmap_lgl(plan, function(var, year, set) {
    if (set == "monthly") {
      download_one(pick_url(urls_monthly, var, year), monthly_dir, var, year)
    } else {
      download_one(pick_url(urls_annual, var, year), annual_dir, var, year)
    }
  })
  
  plan$ok <- results
  n_fail <- sum(!plan$ok)
  message("\nDone. ", sum(plan$ok), " succeeded, ", n_fail, " failed.")
  if (n_fail > 0) {
    message("Failed (first 20):")
    print(utils::head(plan[!plan$ok, ], 20))
  }
}

length(list.files(monthly_dir, pattern = "_na_.*\\.tif$"))  # expect ~132
length(list.files(annual_dir,  pattern = "_na_.*\\.tif$"))  # expect ~132