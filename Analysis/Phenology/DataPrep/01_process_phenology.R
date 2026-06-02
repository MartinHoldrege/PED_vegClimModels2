# Process eVIIRS phenology metrics to the Daymet 1 km grid.
#
# Per metric: decode NoData -> project each year to Daymet 1 km -> collapse
# across years (circular mean for timing metrics, arithmetic mean otherwise)
# -> rescale -> write a single 1 km layer.
#
# Final step: enforce SOST <= MAXT <= EOST on the averaged layers, using the
# circular wrap (SOST-365, EOST+365) where that brings them into order and
# within range, else NA. Counts of NA-ed cells are written out.

library(terra)

source('Functions/init.R')
source_functions()                       # daymet_grid, crs_daymet

raw_dir <- file.path(paths$large, "Data_raw", "phenology")
out_dir <- file.path(paths$large, "Data_processed", "phenology")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

info_tbl <- phenology_metric_info()
files    <- find_phenology_files(raw_dir)


# Process one metric: decode, project each year, collapse, rescale ------------
process_one_metric <- function(metric_code) {
  info <- filter(info_tbl, metric == metric_code)
  mf   <- filter(files, metric == metric_code) |> arrange(year)
  message("Processing ", metric_code, " (", nrow(mf), " years)")
  
  per_year <- map(mf$path, \(p) {
    terra::rast(p) |>
      decode_phenology(info) |>
      project_phenology_year(info, daymet_grid)
  })
  
  stk <- terra::rast(map(per_year, \(x) x[[info$metric]]))
  collapsed <- collapse_years(stk, info)        # circular or arithmetic mean
  rescale_phenology(collapsed, info)            # single layer, real units
}

# Build all metric layers into one named stack.
metric_layers <- map(info_tbl$metric[1], process_one_metric)
pheno <- terra::rast(metric_layers)
names(pheno) <- info_tbl$metric


# Enforce ordering SOST <= MAXT <= EOST on the averaged layers ----------------
# MAXT is the anchor (bounded 1-365, never wraps). SOST after MAXT is assumed
# to have wrapped from the prior year: try SOST-365, keep if >= -150 else NA.
# EOST before MAXT is assumed to wrap into next year: try EOST+365, keep if
# <= 450 else NA. Where MAXT is NA, SOST and EOST cannot be ordered -> NA.
enforce_ordering <- function(sost, maxt, eost) {
  # SOST: pull back a year if it sits after the peak.
  sost_fix <- terra::ifel(sost > maxt, sost - 365, sost)
  sost_fix <- terra::ifel(sost_fix < -150, NA, sost_fix)
  
  # EOST: push forward a year if it sits before the peak.
  eost_fix <- terra::ifel(eost < maxt, eost + 365, eost)
  eost_fix <- terra::ifel(eost_fix > 450, NA, eost_fix)
  
  # No anchor -> cannot order.
  sost_fix <- terra::ifel(is.na(maxt), NA, sost_fix)
  eost_fix <- terra::ifel(is.na(maxt), NA, eost_fix)
  
  list(SOST = sost_fix, EOST = eost_fix)
}

fixed <- enforce_ordering(pheno[["SOST"]], pheno[["MAXT"]], pheno[["EOST"]])

# Count cells newly set to NA by the ordering step (were valid before).
n_sost_na <- terra::global(
  (!is.na(pheno[["SOST"]])) & is.na(fixed$SOST), "sum", na.rm = TRUE
)[[1]]
n_eost_na <- terra::global(
  (!is.na(pheno[["EOST"]])) & is.na(fixed$EOST), "sum", na.rm = TRUE
)[[1]]

write_csv(
  tibble(metric = c("SOST", "EOST"), n_set_na = c(n_sost_na, n_eost_na)),
  file.path(out_dir, "ordering_na_counts.csv")
)

# Replace SOST/EOST with the ordered versions.
pheno[["SOST"]] <- fixed$SOST
pheno[["EOST"]] <- fixed$EOST


# Write final layers ----------------------------------------------------------
terra::writeRaster(pheno, file.path(out_dir, "diagnostics",
                                    "phenology_eVIIRS_1000m.tif"),
                   overwrite = TRUE)

