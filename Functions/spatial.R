# various functions for working with spatial data



crs_daymet <-  "PROJCRS[\"unnamed\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"unknown\",\n            ELLIPSOID[\"Spheroid\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1,\n                    ID[\"EPSG\",9001]]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Lambert Conic Conformal (2SP)\",\n        METHOD[\"Lambert Conic Conformal (2SP)\",\n            ID[\"EPSG\",9802]],\n        PARAMETER[\"Latitude of false origin\",42.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",25,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",60,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]"


# same extent as used for exports in GEE
daymet_grid <- terra::rast(
  xmin = -1950750, xmax = 2428250,
  ymin = -1786500, ymax = 946500,
  resolution = 1000,        # 1 km
  crs = crs_daymet
)

#' Align rasters to a common extent
#'
#' Checks that all rasters share the same CRS, computes the intersection
#' of their extents, and crops each to that common extent.
#'
#' @param ... SpatRaster objects to align.
#' @param rast_list Alternatively, a named list of SpatRaster objects.
#'
#' @return A named list of cropped SpatRaster objects.
#' @export
align_raster_extents <- function(..., rast_list = NULL) {
  if (is.null(rast_list)) {
    rast_list <- list(...)
    # try to capture names from call
    nms <- as.character(match.call(expand.dots = FALSE)$...)
    if (length(nms) == length(rast_list)) names(rast_list) <- nms
  }
  stopifnot(length(rast_list) >= 2)
  
  # check CRS
  # check CRS
  ref <- rast_list[[1]]
  crs_ok <- purrr::map_lgl(rast_list[-1], \(r) terra::same.crs(r, ref))
  if (!all(crs_ok)) {
    bad <- names(rast_list[-1])[!crs_ok]
    stop("CRS mismatch for: ", paste(bad, collapse = ", "),
         "\n  relative to: ", names(rast_list)[1])
  }
  
  # harmonize CRS strings — use the most common representation
  crs_strings <- purrr::map_chr(rast_list, terra::crs)
  crs_freq <- sort(table(crs_strings), decreasing = TRUE)
  ref_crs <- names(crs_freq)[1]
  
  for (i in seq_along(rast_list)) {
    terra::crs(rast_list[[i]]) <- ref_crs
  }

  # intersect extents
  common_ext <- purrr::reduce(purrr::map(rast_list, terra::ext), terra::intersect)
  
  if (is.null(common_ext)) {
    stop("Rasters have no overlapping extent.")
  }
  
  # cat("Common extent:", as.vector(common_ext), "\n")
  
  # crop each
  purrr::map(rast_list, \(r) terra::crop(r, common_ext))
}

align_raster <- function(r, snap = read_mask()) {
  r <- terra::crop(r, snap)

  # Stop if the CRSs are not actually equivalent
  stopifnot(terra::same.crs(r, snap))
  
  # Equivalent CRS, different WKT text: copy snap's (no resampling)
  terra::crs(r) <- terra::crs(snap)
  
  stopifnot(terra::compareGeom(r, snap, stopOnError = FALSE,
                               rowcol = TRUE, ext = TRUE))
  r
}

# for just a single layer, for determining limits for plotting
raster_quantile <- function(r, prob) {
  stopifnot(nlyr(r) == 1)
  terra::global(r, fun = function(x) quantile(x, probs = prob, na.rm = TRUE))[[1]]
}


#' Assign each row a snap-grid cell
#'
#' @param df Data frame with longitude and latitude columns.
#' @param lon,lat Column names, as strings.
#' @param crs EPSG string for the source coordinates.
#' @return df with a `cell` column; NA where the point falls off the grid.
assign_cell <- function(df, lon, lat, crs, snap = read_mask()) {
  v <- terra::vect(df, geom = c(lon, lat), crs = crs) |>
    terra::project(terra::crs(snap))
  stopifnot('cell_id' %in% names(snap))
  df$cell <- terra::extract(snap, v)$cell_id
  df
}

#' Predict a fitted model onto a raster
#'
#' Generic; methods exist for cwexp biomass fits (`predict_raster.cwexp_fit`)
#' and cover classification models (`predict_raster.cover_classification`).
#' Generics are defined in other scripts under Functions/models/
#' 
#' @param fit Fitted model object.
#' @param rast `SpatRaster` of predictors.
#' @param ... Passed to methods.
#' @return `SpatRaster`.
predict_raster <- function(fit, rast, ...) {
  UseMethod("predict_raster")
}

#' RAP cover at given cell-years
#'
#' @param r SpatRaster with one band per year, named 'year_YYYY'.
#' @param cell_years Data frame with `cell` and `year` columns.
#' @return Numeric vector of cover, one value per row of `cell_years`.
#' @examples
#' r <- rap$tree
#' cell_years <- tibble(cell = c(4190, 4191), year = c(2015L, 2020L))
#' rap_at(r = r, cell_years = cell_years)
rap_at <- function(r, cell_years) {
  # one row per cell-year, one column per year
  vals <- as.matrix(r[cell_years$cell])
  col <- match(paste0("year_", cell_years$year), names(r))
  stopifnot(!anyNA(col))
  vals[cbind(seq_len(nrow(cell_years)), col)]
}

#' Change of each GCM's prediction from the historical prediction
#'
#' CONUS (whole-raster) means, area-weighted (cells of a lon/lat grid differ
#' in area), over the cells where both the GCM and the historical predictions
#' have a value. GCMs are ranked by delta and the lowest, median (rank
#' ceiling(n / 2)) and highest are labelled. As in 03_predict_raster.R.
#'
#' @param preds SpatRaster, one layer per GCM, named by GCM.
#' @param hist SpatRaster, one layer: the historical prediction (e.g.
#'   gridMET), on the same grid.
#' @return Tibble, one row per GCM: gcm, mean_gcm, delta (mean of GCM -
#'   historical), rank (1 = lowest delta) and role ("lowest", "median",
#'   "highest", else NA).
#' @examples
#' hist <- terra::rast(nrows = 2, ncols = 2, vals = c(10, 20, NA, 40))
#' preds <- c(hist + 1, hist - 2, hist + 5)
#' names(preds) <- c("a", "b", "c")
#' gcm_deltas(preds = preds, hist = hist)
gcm_deltas <- function(preds, hist) {
  stopifnot(terra::nlyr(hist) == 1, terra::compareGeom(preds, hist),
            !anyDuplicated(names(preds)))
  area <- terra::cellSize(hist, unit = "km")
  area_mean <- function(r) {
    a <- terra::mask(area, r)
    terra::global(r * a, "sum", na.rm = TRUE)[[1]] /
      terra::global(a, "sum", na.rm = TRUE)[[1]]
  }
  
  purrr::map(names(preds), \(g) {
    delta_r <- preds[[g]] - hist
    tibble::tibble(gcm = g,
                   mean_gcm = area_mean(terra::mask(preds[[g]], delta_r)),
                   delta = area_mean(delta_r))
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(rank = rank(delta, ties.method = "first"),
                  role = dplyr::case_when(
                    rank == 1 ~ "lowest",
                    rank == ceiling(dplyr::n() / 2) ~ "median",
                    rank == dplyr::n() ~ "highest"))
}
