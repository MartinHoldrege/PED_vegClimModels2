# Process eVIIRS phenology metrics to the Daymet 1 km grid.
#
# Pipeline, per metric:
#   1. find the 4 yearly 375 m GeoTIFFs
#   2. decode NoData/water codes -> NA, clamp to valid range (raw units)
#   3. project each year to the Daymet 1 km grid (circular for SOST/EOST),
#      carrying a per-year valid-fraction layer
#   4. collapse across years (mean AND median; circular mean for timing)
#   5. rescale to real units and write outputs
#
# Diagnostics are written alongside so the mean-vs-median and seam questions
# are answered from the data, not assumed. Run order is grouped so partial
# results are inspectable.

library(terra)

source('Functions/init.R')
source_functions()                      # provides daymet_grid, crs_daymet

# Paths -----------------------------------------------------------------------
raw_dir <- file.path(paths$large, "Data_raw", "phenology")
out_dir <- file.path(paths$large, "Data_processed", "phenology")
diag_dir <- out_dir
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

info_tbl <- phenology_metric_info()
files <- find_phenology_files(raw_dir)

# ---------------------------------------------------------------------------
# STEP 0: Inventory each raw file (grid, extent, crs, value range) ----------
# The two metadata files implied the yearly grids may differ; verify against
# the actual rasters before trusting any cross-year operation.
# ---------------------------------------------------------------------------
inventory <- files |>
  mutate(meta = map(path, \(p) {
    r <- terra::rast(p)
    tibble(
      nrow = nrow(r), ncol = ncol(r),
      xmin = terra::xmin(r), xmax = terra::xmax(r),
      ymin = terra::ymin(r), ymax = terra::ymax(r),
      crs_name = terra::crs(r, describe = TRUE)$name,
      vmin = terra::minmax(r)[1], vmax = terra::minmax(r)[2]
    )
  })) |>
  unnest(meta)

write_csv(inventory, file.path(diag_dir, "00_raw_inventory.csv"))
print(inventory, n = Inf)

# Flag cross-year grid mismatches per metric (extent/dim differences mean the
# yearly rasters are NOT pixel-aligned and must each be projected independently,
# which this pipeline does -- but we want to know).
grid_check <- inventory |>
  group_by(metric) |>
  summarise(
    n_distinct_dims   = n_distinct(paste(nrow, ncol)),
    n_distinct_extent = n_distinct(paste(xmin, xmax, ymin, ymax)),
    .groups = "drop"
  )
print(grid_check)


# ---------------------------------------------------------------------------
# STEP 1-4: Process one metric end to end -----------------------------------
# Wrapped in a function so it can be run per metric (memory-friendly) and
# re-run individually. Each metric writes: the collapsed multi-layer diagnostic
# raster, the final chosen layer, and a mean-vs-median summary row.
# ---------------------------------------------------------------------------

process_one_metric <- function(metric_code) {
  info <- filter(info_tbl, metric == metric_code)
  mf <- filter(files, metric == metric_code) |> arrange(year)
  message("Processing ", metric_code, " (", nrow(mf), " years)")

  # Decode + project each year; keep metric layer and valid_frac separately.
  per_year <- map(seq_len(nrow(mf)), \(i) {
    r <- terra::rast(mf$path[i])
    r <- decode_phenology(r, info)
    project_phenology_year(r, info, daymet_grid)
  })

  metric_stk <- terra::rast(map(per_year, \(x) x[[info$metric]]))
  names(metric_stk) <- paste0(metric_code, "_", mf$year)
  vf_stk <- terra::rast(map(per_year, \(x) x[["valid_frac"]]))
  names(vf_stk) <- paste0("validfrac_", mf$year)

  # Collapse across years (raw units): mean, median, diffs, (circular).
  collapsed <- collapse_years(metric_stk, info)

  # Per-year valid fraction summary + an overall "n valid years" layer.
  n_valid_years <- terra::app(metric_stk, \(x) sum(!is.na(x)))
  names(n_valid_years) <- "n_valid_years"
  mean_valid_frac <- terra::app(vf_stk, mean, na.rm = TRUE)
  names(mean_valid_frac) <- "mean_valid_frac"

  # Rescale the central-tendency layers to real units for the final product.
  final_mean   <- rescale_phenology(collapsed[["mean"]],   info)
  final_median <- rescale_phenology(collapsed[["median"]], info)
  names(final_mean) <- paste0(metric_code, "_mean")
  names(final_median) <- paste0(metric_code, "_median")
  final_layers <- c(final_mean, final_median)

  if (isTRUE(info$is_circular)) {
    final_circ <- rescale_phenology(collapsed[["circ_mean"]], info)
    names(final_circ) <- paste0(metric_code, "_circ_mean")
    final_layers <- c(final_layers, final_circ)
  }

  # Write: diagnostic stack (raw units) + final layers (real units) + support.
  terra::writeRaster(
    c(collapsed, n_valid_years, mean_valid_frac),
    file.path(out_dir, paste0(metric_code, "_diagnostics_1km.tif")),
    overwrite = TRUE
  )
  terra::writeRaster(
    final_layers,
    file.path(out_dir, paste0(metric_code, "_1km.tif")),
    overwrite = TRUE
  )

  # Numeric summary of the mean-vs-median divergence (the key decision input).
  ad <- collapsed[["mean_median_absdiff"]]
  q <- terra::global(ad, fun = \(x) quantile(x, c(.5, .9, .99, 1), na.rm = TRUE))
  summ <- tibble(
    metric = metric_code,
    type = info$type,
    md_absdiff_p50 = q[[1]], md_absdiff_p90 = q[[2]],
    md_absdiff_p99 = q[[3]], md_absdiff_max = q[[4]]
  )
  if (info$is_circular) {
    cl <- terra::global(collapsed[["circ_linear_absdiff"]],
                        fun = \(x) quantile(x, c(.9, .99, 1), na.rm = TRUE))
    summ <- summ |>
      mutate(circ_linear_p90 = cl[[1]], circ_linear_p99 = cl[[2]],
             circ_linear_max = cl[[3]])
  }
  summ
}

# Run all metrics, collect the mean-vs-median / circular  <- .
summaries <- map_dfr(info_tbl$metric, process_one_metric)
write_csv(summaries, file.path(diag_dir, "01_collapse_summary.csv"))
print(summaries)

# ---------------------------------------------------------------------------
# DECISION GUIDE (read the summary, then decide) ----------------------------
# - md_absdiff_*: how far mean and median diverge across years, in raw units.
#   For NDVI metrics (raw 101-200) a p99 of a few units is negligible -> use
#   mean. For timing (days) a p99 of a few days is negligible; tens of days in
#   a meaningful number of cells signals outlier/seam contamination -> prefer
#   median, or mask high-divergence cells.
# - circ_linear_*: for SOST/EOST, how far the linear cross-year mean sits from
#   the circular cross-year mean (days). Large only where the seam bites. If
#   small CONUS-wide, the simple linear mean is safe and you can ignore the
#   circular machinery for the cross-year step (spatial step already handled).
#
# If both diagnostics are small -> use the *_mean layers (your stated preference).
# ---------------------------------------------------------------------------
