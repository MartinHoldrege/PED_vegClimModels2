# various functions for working with spatial data



crs_daymet <-  "PROJCRS[\"unnamed\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"unknown\",\n            ELLIPSOID[\"Spheroid\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1,\n                    ID[\"EPSG\",9001]]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Lambert Conic Conformal (2SP)\",\n        METHOD[\"Lambert Conic Conformal (2SP)\",\n            ID[\"EPSG\",9802]],\n        PARAMETER[\"Latitude of false origin\",42.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",25,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",60,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]"


# same extent as used for exports in GEE
daymet_grid <- terra::rast(
  xmin = -1950750, xmax = 2428250,
  ymin = -1785500, ymax = 945500,
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
  
  cat("Common extent:", as.vector(common_ext), "\n")
  
  # crop each
  purrr::map(rast_list, \(r) terra::crop(r, common_ext))
}
