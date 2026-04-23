# Mapping functions for CONUS raster visualization

# at the top of map_plots.R
.map_cache <- new.env(parent = emptyenv())

#' Plot a SpatRaster as a ggplot map with US state boundaries
#'
#' Creates a ggplot-based map from a SpatRaster with consistent CONUS
#' extent and state boundary basemap. Supports continuous and categorical
#' data, multi-layer faceting, and flexible color scales.
#'
#' @param rast A SpatRaster. Single-layer for one map, multi-layer for
#'   faceted maps.
#' @param colorscale A ggplot2 scale object (e.g., `scale_fill_viridis_c()`,
#'   `scale_fill_distiller()`). If NULL, uses `scale_fill_viridis_c()`.
#' @param title Character; plot title.
#' @param legend_name Character; legend title. Default uses layer name.
#' @param basemap_color Character; color for state boundaries.
#' @param basemap_size Numeric; line width for state boundaries.
#' @param na_color Character; color for NA values.
#' @param maxcell Integer; maximum number of cells to plot (for speed).
#'   Default 1e6. Set higher for publication quality.
#' @param states_sf Optional sf object of state boundaries. If NULL,
#'   loaded from `spData::us_states` and transformed to raster CRS.
#'
#' @return A ggplot object.
#' @examples
#' # r <- terra::rast(system.file("ex/elev.tif", package = "terra"))
#' # plot_map_conus(r, title = "Elevation")
#' @export
plot_map_conus <- function(rast,
                           colorscale = NULL,
                           title = NULL,
                           legend_name = NULL,
                           basemap_color = "white",
                           basemap_size = 0.2,
                           na_color = "grey90",
                           maxcell = 1e6,
                           states_sf = NULL) {
  
  stopifnot(
    inherits(rast, "SpatRaster"),
    requireNamespace("tidyterra", quietly = TRUE)
  )
  
  # get state boundaries
  if (is.null(states_sf)) {
    states_sf <- get_or_cache_states(crs = terra::crs(rast))
  }
  
  # default color scale
  if (is.null(colorscale)) {
    colorscale <- ggplot2::scale_fill_viridis_c(
      name = legend_name,
      na.value = na_color,
      option = "viridis"
    )
  }
  
  # build plot
  g <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = rast, maxcell = maxcell) +
    colorscale +
    ggplot2::geom_sf(data = states_sf,
                     fill = NA,
                     color = basemap_color,
                     linewidth = basemap_size) +
    ggplot2::coord_sf(
      xlim = terra::ext(rast)[1:2],
      ylim = terra::ext(rast)[3:4],
      crs = terra::crs(rast),
      expand = FALSE
    ) +
    ggplot2::labs(title = title) +
    map_theme()
  
  # if multi-layer, facet
  if (terra::nlyr(rast) > 1) {
    g <- g + ggplot2::facet_wrap(~ lyr, ncol = 2)
  }
  
  g
}


#' Minimal clean theme for maps
#'
#' @return A ggplot2 theme object.
#' @export
map_theme <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2),
      legend.position = "right",
      strip.text = ggplot2::element_text(size = 10, face = "bold")
    )
}


#' Get CONUS state boundaries as sf, excluding AK/HI
#'
#' @param crs CRS to transform to (character or terra CRS object).
#'
#' @return An sf object of CONUS state boundaries.
#' @export
get_conus_states <- function(crs) {
  stopifnot(requireNamespace("spData", quietly = TRUE))
  
  states <- spData::us_states
  # spData::us_states already excludes AK, HI, and territories
  sf::st_transform(states, crs = crs)
}



#' Get or cache CONUS state boundaries
#'
#' @param crs CRS to transform to.
#' @return An sf object.
get_or_cache_states <- function(crs) {
  crs_string <- terra::crs(crs, proj = TRUE)
  
  if (!is.null(.map_cache$states) && 
      identical(sf::st_crs(.map_cache$states)$proj4string, crs_string)) {
    return(.map_cache$states)
  }
  
  .map_cache$states <- get_conus_states(crs = crs)
  .map_cache$states
}

# --- Convenience color scales ---

#' Viridis continuous color scale for biomass maps
#' @param name Legend title.
#' @param limits Numeric vector of length 2.
#' @param option Viridis palette option (default "viridis").
#' @param ... Passed to `scale_fill_viridis_c`.
#' @export
colorscale_biomass <- function(name = "Biomass\n(Mg/ha)",
                               limits = NULL,
                               option = "viridis",
                               ...) {
  ggplot2::scale_fill_viridis_c(
    name = name,
    limits = limits,
    na.value = "grey90",
    option = option,
    oob = scales::oob_squish,
    ...
  )
}


#' Diverging color scale for residual/difference maps
#' @param name Legend title.
#' @param limits Symmetric limits (e.g., c(-50, 50)).
#' @param midpoint Midpoint of the diverging scale.
#' @param ... Passed to `scale_fill_gradient2`.
#' @export
colorscale_diverging <- function(name = "Difference",
                                 limits = NULL,
                                 midpoint = 0,
                                 ...) {
  ggplot2::scale_fill_gradient2(
    name = name,
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = midpoint,
    limits = limits,
    na.value = "grey90",
    oob = scales::oob_squish,
    ...
  )
}


#' Continuous color scale for cover (0-1 or 0-100%)
#' @param name Legend title.
#' @param limits Numeric limits.
#' @param ... Passed to `scale_fill_viridis_c`.
#' @export
colorscale_cover <- function(name = "Cover",
                             limits = c(0, 1),
                             ...) {
  ggplot2::scale_fill_viridis_c(
    name = name,
    limits = limits,
    na.value = "grey90",
    option = "mako",
    direction = -1,
    oob = scales::oob_squish,
    ...
  )
}

# example usage
if(FALSE) {
  # quick test of plot_map_conus using load_conus_rasters output
  source("Functions/init.R")
  source_functions()
  
  obj <- load_conus_rasters()

  # single layer: herbaceous biomass
  plot_map_conus(obj$herb_bio, 
                 colorscale = colorscale_biomass(limits = c(0, 5)),
                 title = "RAP Herbaceous AGB")
  
  # single layer: one climate variable
  plot_map_conus(obj$climate[["MAT"]], 
                 title = "Mean Annual Temperature")
  
  # multi-layer: cover (auto-faceted)
  plot_map_conus(obj$cover, 
                 colorscale = colorscale_cover(),
                 title = "RAP Cover")
}