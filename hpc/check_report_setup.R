# check_report_setup.R
#
# Checks that the cover evaluation reports (Analysis/Cover/Evaluate/
# 01_*_diagnostics.Rmd, rendered by main_cover.R) can be made on the HPC:
#   1. R packages the reports use; missing CRAN packages are installed to the
#      user library if install_missing = TRUE
#   2. pandoc (rmarkdown::render() needs it)
#   3. a PNG device that works without a display
#   4. a test render of a small Rmd (pandoc and PNG figures together)
#   5. the repo's functions load and plot_map_conus() draws a map
# Stops at the first failure, saying what to fix.
#
# Run from the repository root:
#   module load cray-R-spatial/4.6
#   Rscript hpc/check_report_setup.R
# Installing packages needs internet, which compute nodes may not have, so run
# it on a login node the first time.
#
# If pandoc is missing (step 2), check `module avail pandoc`. Otherwise, on a
# login node:
#   mkdir -p ~/software && cd ~/software
#   wget https://github.com/jgm/pandoc/releases/download/3.5/pandoc-3.5-linux-amd64.tar.gz
#   tar xzf pandoc-3.5-linux-amd64.tar.gz
# and before running R (and in hpc/run_cover.sh):
#   export RSTUDIO_PANDOC=$HOME/software/pandoc-3.5/bin
#
# November 2026

install_missing <- TRUE

# 1. packages ---------------------------------------------------------------
# used directly in the Rmds and main_cover.R; the repo's function files may
# need others, which step 5 catches

pkgs <- c("rmarkdown", "knitr", "callr", "patchwork",
          "PresenceAbsence", "glmnet", "ranger", "mgcv", "grpreg",
          "grplasso", "sf", "terra", "scales", "viridisLite", "viridis", 'yardstick')

missing <- setdiff(pkgs, rownames(installed.packages()))

if (length(missing) > 0) {
  message("missing packages: ", paste(missing, collapse = ", "))
  if (!install_missing) stop("install them, or set install_missing = TRUE")

  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")

  still_missing <- setdiff(missing,
                           rownames(installed.packages(noCache = TRUE)))
  if (length(still_missing) > 0) {
    stop("failed to install: ", paste(still_missing, collapse = ", "))
  }
}
message("1. packages ok")

# 2. pandoc -----------------------------------------------------------------

if (!rmarkdown::pandoc_available()) {
  stop("pandoc not found: see the header of this script")
}
message("2. pandoc ", rmarkdown::pandoc_version(), " at ",
        rmarkdown::find_pandoc()$dir)

# 3. PNG device -------------------------------------------------------------
# knitr draws figures with png(), whose type comes from
# getOption("bitmapType"); without a display only "cairo" works

#' Whether png() can write a file with a given device type
#'
#' @param type png() device type, e.g. "cairo".
#' @return TRUE if a PNG was written.
#' @examples
#' type <- "cairo"
#' png_works(type = type)
png_works <- function(type) {
  f <- tempfile(fileext = ".png")
  tryCatch({
    png(f, type = type)
    plot(1:10)
    dev.off()
    file.exists(f)
  }, error = function(e) FALSE)
}

if (!png_works(getOption("bitmapType"))) {
  if (!(capabilities("cairo") && png_works("cairo"))) {
    stop("no working PNG device: ask HPC support about cairo, or install ragg")
  }
  # carry on with cairo, so the steps below test the setup the fix gives
  options(bitmapType = "cairo")
  message("PNG works only with cairo: add options(bitmapType = \"cairo\") ",
          "to the top of main_cover.R")
}
message("3. PNG device ok (", getOption("bitmapType"), ")")

# 4. test render ------------------------------------------------------------

test_dir <- tempfile("render_test_")
dir.create(test_dir)
test_rmd <- file.path(test_dir, "test.Rmd")
writeLines(c("---", "title: test", "output: html_document", "---", "",
             "```{r}",
             "library(ggplot2)",
             "ggplot(mtcars, aes(wt, mpg)) + geom_point()",
             "```"),
           test_rmd)

out <- rmarkdown::render(test_rmd, quiet = TRUE)
stopifnot(file.exists(out))
message("4. test render ok: ", out)

# 5. repo functions and a map -----------------------------------------------
# needs the snap mask (hpc/sync_cover_to_hpc.sh, sync_predict_data); also
# shows whether the map functions read any other file or need internet

source("Functions/init.R")
source_functions()

g <- plot_map_conus(read_mask(),
                    colorscale = ggplot2::scale_fill_viridis_c(),
                    title = "check", maxcell = 1e4)
f <- tempfile(fileext = ".png")
ggplot2::ggsave(f, g, width = 6, height = 4, dpi = 72)
stopifnot(file.exists(f))
message("5. repo functions and map ok")

message("\nAll checks passed. Next, render one real report (its fit and ",
        "predictions must exist):\n",
        "  Rscript Analysis/Cover/main_cover.R --run_index=<i> ",
        "--run_fit=FALSE --run_predict=FALSE --run_evaluate=TRUE")