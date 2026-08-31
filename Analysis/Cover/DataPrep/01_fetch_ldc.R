# 01_fetch_ldc.R
# Download Landscape Data Commons header and line-point intercept data.
#
# Resumable: each batch is cached to disk and skipped on rerun. Safe to
# kill and restart at any point. Expect ~10 hours for a full cold run
# (~106k plots, ~35M LPI rows).

library(tidyverse)
source("Functions/init.R")

# ---- settings -----------------------------------------------------------
ldc_dir    <- file.path(paths$large, "Data_raw/LandscapeDataCommonsDat")
lpi_dir    <- file.path(ldc_dir, "lpi_batches")
batch_size <- 200      # 300 works but sits closer to the server timeout
min_split  <- 5        # give up on a batch once it is this small
api_delay  <- 500      # ms between requests; do not go below ~250

dir.create(lpi_dir, showWarnings = FALSE, recursive = TRUE)

# ---- header -------------------------------------------------------------
hdr_path <- file.path(ldc_dir, "header_raw.rds")

if (file.exists(hdr_path)) {
  hdr_all <- readRDS(hdr_path)
  message("header: ", nrow(hdr_all), " rows (cached)")
} else {
  hdr_all <- trex::fetch_ldc(data_type = "header")
  stopifnot(is.data.frame(hdr_all), nrow(hdr_all) > 0)
  saveRDS(hdr_all, hdr_path, compress = "xz")
  message("header: ", nrow(hdr_all), " rows (fetched)")
}

# ---- LPI ----------------------------------------------------------------

#' Fetch LPI for a set of PrimaryKeys, halving the batch on failure
#'
#' The LDC API returns 504 on requests that take too long server-side, and
#' occasional plots fail individually. Splitting recursively handles both.
#'
#' @param keys Character vector of PrimaryKeys.
#' @param min_size Stop splitting below this size and abandon the batch.
#' @return A list with two elements:
#'   `data`      A data frame of the LPI records retrieved (0 rows if every
#'               queried plot legitimately has no LPI data).
#'   `abandoned` Character vector of PrimaryKeys that could not be retrieved
#'               due to repeated API failures. Empty when the query succeeded
#'               for all keys, even if some keys returned no rows. This lets
#'               the caller distinguish "plot has no LPI data" (cacheable)
#'               from "plot request failed" (retry).
fetch_lpi_keys <- function(keys, min_size = min_split) {
  x <- try(
    trex::fetch_ldc(
      data_type = "lpi",
      query_parameters = list("PrimaryKey" = list("=" = keys)),
      take = 10000,
      delay = api_delay,
      timeout = 1800
    ),
    silent = TRUE
  )
  
  if (!inherits(x, "try-error")) {
    return(list(data = x, abandoned = character(0)))
  }
  
  if (length(keys) <= min_size) {
    message("  giving up on ", length(keys), " keys")
    return(list(data = NULL, abandoned = keys))
  }
  
  message("  split: ", length(keys), " -> ", ceiling(length(keys) / 2))
  parts <- keys |>
    split(seq_along(keys) > length(keys) / 2) |>
    map(fetch_lpi_keys, min_size = min_size)
  
  list(
    data      = map(parts, "data") |> list_rbind(),
    abandoned = map(parts, "abandoned") |> unlist(use.names = FALSE)
  )
}

#' Fetch and cache one batch of plots
#'
#' Skips the batch if its cache file exists. Writes to a `.partial` file
#' first so an interrupted write leaves no truncated cache. Batches that
#' fail entirely are not cached and will be retried on the next run.
#'
#' @param keys Character vector of PrimaryKeys.
#' @param batch_id Batch number, used in the file name.
fetch_lpi_batch <- function(keys, batch_id) {
  out <- file.path(lpi_dir, sprintf("lpi_batch-%05d.rds", as.integer(batch_id)))
  if (file.exists(out)) return(invisible(NULL))
  
  res <- fetch_lpi_keys(keys)
  
  # Only cache once every requested plot was successfully queried. A key is
  # "abandoned" only when its request failed repeatedly; a key that returns
  # zero rows (no LPI data for that plot) is a valid, complete result and is
  # cached so it is not re-queried every run. Any abandoned keys mean the
  # batch is incomplete and is left uncached for retry next run.
  if (length(res$abandoned) > 0) {
    message(batch_id, ": INCOMPLETE, ", length(res$abandoned), "/", length(keys),
            " plots failed, not caching, will retry next run")
    return(invisible(NULL))
  }
  
  x <- res$data
  n_rows <- if (is.null(x)) 0L else nrow(x)
  
  tmp <- paste0(out, ".partial")
  saveRDS(x, tmp)
  file.rename(tmp, out)
  
  message(batch_id, ": ", n_rows, " rows, ",
          n_distinct(x$PrimaryKey), "/", length(keys), " plots with data")
  invisible(NULL)
}

all_keys    <- sort(unique(hdr_all$PrimaryKey))
key_batches <- split(all_keys, ceiling(seq_along(all_keys) / batch_size))

message("plots: ", length(all_keys), " in ", length(key_batches), " batches")
message("cached: ", length(list.files(lpi_dir, pattern = "^lpi_batch-\\d+\\.rds$")))

# clear any partial writes from a previous interrupted run
file.remove(list.files(lpi_dir, pattern = "\\.partial$", full.names = TRUE))

iwalk(key_batches, fetch_lpi_batch)

# ---- report -------------------------------------------------------------
done <- list.files(lpi_dir, pattern = "^lpi_batch-\\d+\\.rds$")

message("\ncomplete: ", length(done), "/", length(key_batches), " batches")

if (length(done) < length(key_batches)) {
  message("rerun this script to retry the missing batches")
} else {
  writeLines(format(Sys.time()), file.path(ldc_dir, "lpi_fetch_complete.txt"))
}