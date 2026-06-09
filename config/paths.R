

# path to where large files are stored (update paths as needed), that aren't in working directory
# note the folder structure there should be identical to the local
# folder structure, except just store big files there so don't 
# take up local storage
paths <- list()
if(dir.exists("E:/USGS")) {
  .base <- "E:/USGS"
} else if(dir.exists("D:/USGS")) {
  .base <- "D:/USGS"
} else {
  warning('directory not found')
}
paths$large <- file.path(.base, "large_files/PED_vegClimModels2")

# used in a few cases to avoid copying raw data
# unnecessarily, this is the directory created by Alice:
paths$large0 <- file.path(.base, "large_files/PED_vegClimModels")