

# path to where large files are stored (update paths as needed), that aren't in working directory
# note the folder structure there should be identical to the local
# folder structure, except just store big files there so don't 
# take up local storage
paths <- list()
if (dir.exists("/caldera/hovenweep")) {
  # on the HPC, large files live inside the repo itself
  paths$large <- getwd()
  paths$large0 <- getwd()
} else {
  if (dir.exists("D:/USGS")) {
    .base <- "D:/USGS"
  } else if (dir.exists("E:/USGS")) {
    .base <- "E:/USGS"
  } else {
    stop("large-file directory not found")
  }
  paths$large <- file.path(.base, "large_files/PED_vegClimModels2")
  paths$large0 <- file.path(.base, "large_files/PED_vegClimModels")
}