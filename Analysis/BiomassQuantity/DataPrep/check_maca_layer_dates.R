# //////////////////////////////////////////////////////////////////////////
# check_maca_layer_dates.R
#
# Quick standalone check to confirm that terra reads MACA layer dates the way
# 01_get_projected-climate-data.R expects, i.e. that str_sub(names, 1, 4) yields
# a 4-digit year. Run this once per model (or just once) before the full job.
#
# Assumes the repo is initialised: paths$large0 and (for the read test) the
# MACA files exist. Does NOT need daymet_grid.
# //////////////////////////////////////////////////////////////////////////

source('Functions/init.R')
library(terra)
library(stringr)

# ---- match these to the main script ----
climate_model <- "BNU-ESM"     # or "IPSL-CM5A-MR"
rcp           <- "rcp85"
chunk         <- "2066_2085"   # just test one chunk
maca_var      <- "tasmin"


source_dir <- file.path(paths$large0, 'Data_raw/macaClimateProjections/Data')
f <- file.path(source_dir,
               paste0("macav2livneh_", maca_var, "_", climate_model,
                      "_r1i1p1_", rcp, "_", chunk, "_CONUS_monthly.nc"))

stopifnot(file.exists(f))
r <- terra::rast(f)

# 1. How many layers, and what does terra::time() return?
cat("n layers:", terra::nlyr(r), "\n")
cat("time() class:", class(terra::time(r)), "\n")
cat("first few time() values:\n"); print(utils::head(terra::time(r)))

# 2. Apply the same naming the main script uses, then extract the year.
names(r) <- terra::time(r)
cat("\nfirst few layer names after names(r) <- time(r):\n")
print(utils::head(names(r)))

yr <- str_sub(names(r), 1, 4)
cat("\nfirst few extracted year strings:\n"); print(utils::head(yr))

# 3. The checks that matter:
#    - every extracted year must be 4 digits
#    - as.integer() must not introduce NAs
ok_format <- all(str_detect(yr, "^\\d{4}$"))
yr_int    <- suppressWarnings(as.integer(yr))
ok_int    <- !any(is.na(yr_int))

cat("\nall years 4-digit:", ok_format, "\n")
cat("as.integer() produced no NAs:", ok_int, "\n")
cat("year range seen:", paste(range(yr_int), collapse = " - "), "\n")
cat("layers per year (should all be 12):\n"); print(table(yr_int))

if (ok_format && ok_int) {
  cat("\nPASS: year extraction works as the main script assumes.\n")
} else {
  cat("\nFAIL: layer naming is not what the main script expects.\n",
      "Inspect names(r) above and adjust the year-extraction step\n",
      "(layer_years <- ...) in 07_GetForecastedClimateData.R accordingly.\n")
}

