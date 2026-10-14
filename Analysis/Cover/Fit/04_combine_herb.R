# 04_combine_herb.R
#
# Herbaceous cover by functional type (C3 grass, C4 grass, forb), for current
# and future climate. Combines four models:
#   herbaceous cover (cover/herbaceous)
#   C3, C4 and forb shares of herbaceous cover (proportion/c3, c4, forb)
# The combination (which version of each model) is
# cover_specs$combined$herbaceous[[--vmc]].
#
# The three shares are fit separately, so they don't sum to 100%; they are
# rescaled to do so, and cover of each type is its rescaled share of the
# predicted herbaceous cover (split_cover()). The share_sum layer keeps their
# sum before rescaling, to show how much rescaling did.
#
# Tree cover by leaf type is not written here: tree_leaf_cover() computes it
# on the fly from the combined tree cover and the needleleaf share
# (cover_specs$combined$tree_all).
#
# Inputs (Data_processed/CoverData/Predictions/, from 03_predict_raster.R):
#   <cover_type>_<cover_model>_<vc>-<vmc>_<scenario>.tif for the four models
#
# Outputs (paths from combined_path()):
#   Data_processed/CoverData/Predictions/combined_herbaceous_<vmc>_<scenario>.tif
#     layers: share_c3, share_c4, share_forb (% of herbaceous cover, sum to
#     100), cover_c3, cover_c4, cover_forb (% cover), share_sum (%)
#   Data_processed/CoverData/Fit/combined_herbaceous_<vmc>_components.csv
#     the component models
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

# spec --------------------------------------------------------------------

stopifnot(opt$cover_type == "combined", opt$cover_model == "herbaceous")
spec <- cover_specs$combined$herbaceous[[opt$vmc]]
if (is.null(spec)) stop("no spec for cover_specs$combined$herbaceous$", opt$vmc)

comp <- spec$components  # one row per component model
print(comp)

components_file <- combined_path("herbaceous", opt$vmc, "components")
dir.create(dirname(components_file), recursive = TRUE, showWarnings = FALSE)
write_csv(comp, components_file)

scenarios <- c("current", "BNU-ESM", "IPSL-CM5A-MR")
types <- c("c3", "c4", "forb")

# prediction of one component for one scenario: % cover, or % of
# herbaceous cover for the shares
read_component <- function(component, scenario) {
  read_component_pred(comp[comp$component == component, ], scenario)
}

# combine and write, one scenario at a time --------------------------------

walk(scenarios, \(scenario) {
  shares <- map(types, \(type) read_component(type, scenario)) |>
    terra::rast()
  names(shares) <- types
  herb <- read_component("herbaceous", scenario)

  out <- split_cover(shares, herb)
  terra::writeRaster(out, combined_path("herbaceous", opt$vmc, "prediction",
                                        scenario),
                     overwrite = TRUE, datatype = "FLT4S")
  message("Wrote ", scenario)
})
