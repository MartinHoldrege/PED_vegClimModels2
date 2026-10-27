# Purpose:
# For each MACAv2-METDATA GCM, future scenario (RCP 4.5/8.5 x 2040-2069/
# 2070-2099) and climate variable, calculate the change from the same GCM's
# historical (1971-2000) normals, and summarise it across CONUS and within
# the 7 National Climate Assessment regions: area-weighted mean, 25th
# percentile, median and 75th percentile. Both absolute and percent change are summarised, so either
# can be plotted.
#
# Also prints the GCM colour palette (blue to red, ranked by CONUS mean change
# in MAT under RCP 4.5, 2070-2099) to paste into Functions/viz/fig_params.R
#
# Note: this is for visualize the climate futures, and doesn't not 
# need to be run for training the models
#
# Inputs:
#   Data_processed/WallToWallClimateData/maca/maca_<gcm>_<scenario>_<period>_CLIM.tif
#     (01_summarise_projected_climate_data.R)
#   NCA regions from state boundaries (Functions/nca_regions.R)
# Output:
#   Data_processed/WallToWallClimateData/maca/summaries/maca_delta_summaries.csv
#
# October 2026

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
library(terra)

# params ------------------------------------------------------------------

maca_dir <- file.path(paths$large, "Data_processed/WallToWallClimateData/maca")
out_dir <- file.path(maca_dir, "summaries")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
p_out <- file.path(out_dir, "maca_delta_summaries.csv")

hist_period <- "1971-2000"

# file list ---------------------------------------------------------------

files <- tibble(path = list.files(maca_dir, pattern = "^maca_.*_CLIM\\.tif$",
                                  full.names = TRUE)) |>
  mutate(file = basename(path)) |>
  tidyr::separate_wider_regex(
    file,
    patterns = c("maca_", gcm = ".+", "_",
                 scenario = "historical|rcp45|rcp85", "_",
                 period = "\\d{4}-\\d{4}", "_CLIM\\.tif")
  )

hist_files <- filter(files, scenario == "historical", period == hist_period)
fut_files <- filter(files, scenario != "historical")

# every future file needs its GCM's historical file
stopifnot(nrow(fut_files) > 0,
          all(fut_files$gcm %in% hist_files$gcm),
          !anyDuplicated(hist_files$gcm))

message(n_distinct(fut_files$gcm), " GCMs, ", nrow(fut_files),
        " future scenario-periods")

# grid, cell areas and regions --------------------------------------------
# All MACA outputs are on the same ~4 km lon/lat grid; cell areas differ with
# latitude, so statistics are weighted by cell area.

template <- rast(hist_files$path[1])[[1]]
cell_area <- values(cellSize(template, unit = "km"), mat = FALSE)

# NCA regions; cells with climate data but no region (coastal fringe) count
# towards CONUS only
region_names <- distinct(nca_region_lookup(), region_code, region)
region_code <- values(rasterize_nca_regions(template), mat = FALSE)

has_data <- !is.na(values(template, mat = FALSE))
message(round(100 * mean(is.na(region_code[has_data])), 2),
        "% of cells with data have no NCA region")

# helpers -----------------------------------------------------------------

#' Area-weighted quantiles
#'
#' Returns the smallest value whose cumulative weight share is >= each
#' probability.
#'
#' @param x Numeric vector.
#' @param w Non-negative weights, same length as x.
#' @param probs Probabilities.
#' @return Numeric vector, one value per probability.
#' @examples
#' wtd_quantile(1:4, w = c(1, 1, 1, 5), probs = 0.5) # 4
wtd_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]
  cw <- cumsum(w[o]) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

#' Area-weighted summary of one vector of changes
#'
#' @param x Numeric vector of changes (NAs dropped).
#' @param w Cell areas (km2), same length as x.
#' @return One-row tibble: n_cells, area_km2, mean, q25, median, q75.
#' @examples
#' x = c(1, 2, NA, 4)
#' w = c(1, 1, 1, 2)
#' summarise_wtd(x, w)
summarise_wtd <- function(x, w) {
  keep <- !is.na(x)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(NULL)
  q <- wtd_quantile(x, w, c(0.25, 0.5, 0.75))
  tibble(n_cells = length(x), area_km2 = sum(w),
         mean = sum(x * w) / sum(w),
         q25 = q[1], median = q[2], q75 = q[3])
}

#' Summarise one matrix of changes (cells x variables) by region
#'
#' @param m Matrix of changes, one column per variable.
#' @return Tibble with variable, region, region_code and the summary columns.
summarise_by_region <- function(m) {
  purrr::map_dfr(colnames(m), function(v) {
    x <- m[, v]
    conus <- summarise_wtd(x, cell_area) |>
      mutate(region = "CONUS", region_code = NA_integer_)
    regs <- purrr::map_dfr(sort(unique(na.omit(region_code))), function(rc) {
      i <- which(region_code == rc)
      summarise_wtd(x[i], cell_area[i]) |>
        mutate(region = region_names$region[region_names$region_code == rc],
               region_code = rc)
    })
    bind_rows(conus, regs) |>
      mutate(variable = v, .before = 1)
  })
}

# calculate ---------------------------------------------------------------

out <- purrr::map_dfr(unique(fut_files$gcm), function(g) {
  message("  ", g)
  r_hist <- rast(hist_files$path[hist_files$gcm == g])
  stopifnot(compareGeom(r_hist, template))
  hist <- values(r_hist)
  colnames(hist) <- climate_name_lookup(names(r_hist))

  rows <- filter(fut_files, gcm == g)
  purrr::map_dfr(seq_len(nrow(rows)), function(i) {
    r_fut <- rast(rows$path[i])
    stopifnot(compareGeom(r_fut, template),
              identical(names(r_fut), names(r_hist)))
    fut <- values(r_fut)
    colnames(fut) <- colnames(hist)

    d_abs <- fut - hist
    # percent change is undefined where the historical value is 0
    d_pct <- ifelse(hist == 0, NA, (fut - hist) / abs(hist) * 100)

    bind_rows(
      summarise_by_region(d_abs) |> mutate(delta_type = "abs"),
      summarise_by_region(d_pct) |> mutate(delta_type = "pct")
    ) |>
      mutate(gcm = g, scenario = rows$scenario[i], period = rows$period[i],
             .before = 1)
  })
})

readr::write_csv(out, p_out)
message("Wrote ", p_out)

# GCM colour palette ------------------------------------------------------
# Rank GCMs by CONUS mean change in MAT, RCP 4.5 2070-2099; coolest = blue.
# Paste the printed vector into Functions/constants.R as `gcm_colors`.

gcm_rank <- out |>
  filter(variable == "MAT", delta_type == "abs", region == "CONUS",
         scenario == "rcp45", period == "2070-2099") |>
  arrange(median)
stopifnot(nrow(gcm_rank) == n_distinct(out$gcm))

gcm_colors <- setNames(hcl.colors(nrow(gcm_rank), "RdYlBu", rev = TRUE),
                       gcm_rank$gcm)
print(select(gcm_rank, gcm, median))
dput(gcm_colors)
