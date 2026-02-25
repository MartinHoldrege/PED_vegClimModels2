

# path to where large files are stored (update paths as needed), that aren't in working directory
# note the folder structure there should be identical to the local
# folder structure, except just store big files there so don't 
# take up local storage
paths <- list()
if(dir.exists("E:/USGS")) {
  paths$large <- "E:/USGS/large_files/PED_vegClimModels"
} else  if(dir.exists("D:/USGS")) {
  paths$large <- "D:/USGS/large_files/PED_vegClimModels"
} else {
  warning('directory not found')
}
