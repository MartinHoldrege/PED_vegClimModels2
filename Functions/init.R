# initialization script that should be sourced
# runs code that is likely needed in most scripts



# load main packages ------------------------------------------------------

packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "readr",
  "stringr",
  "forcats",
  "tibble",
  "lubridate"
)

invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}))

options(readr.show_col_types = FALSE)

# source key scripts ------------------------------------------------------

source('Functions/params.R')
source('config/paths.R')
source('Functions/constants.R')


# optional function to use to source all functions
source_functions <- function(path = "./Functions") {
  files <- sort(list.files(path, pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
  files <- stringr::str_subset(files, "init.R$|params.R$", negate = TRUE)
  invisible(lapply(files, source))
}


