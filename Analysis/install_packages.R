# install_repo_packages.R -----------------------------------------------------
# Scan .R and .Rmd files in a repository for package references, then try to
# install whatever is missing. Deliberately uses base R only, since it has to
# run before anything else is installed.
#
# Usage: open a fresh R session at the repo root and source this file.
#        Leave dry_run = TRUE for a first look, then set it to FALSE.

## Settings -------------------------------------------------------------------

repo_dir <- "."          # repo root to scan
dry_run <- FALSE         # TRUE = report only; FALSE = actually install
repos <- "https://cloud.r-project.org"
exclude_dirs <- c(".git", "renv", "packrat", "node_modules", ".Rproj.user")

# Packages that aren't on CRAN; name = package, value = "user/repo".
github_pkgs <- c(
  rSOILWAT2 = "DrylandEcology/rSOILWAT2"
)

## Functions ------------------------------------------------------------------

#' Find R and Rmd files in a directory tree
#'
#' @param dir Directory to search recursively.
#' @param exclude Directory names to skip.
#' @return Character vector of file paths.
find_r_files <- function(dir, exclude = character()) {
  files <- list.files(
    dir,
    pattern = "\\.(R|r|Rmd|rmd|qmd|Rnw)$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(exclude) > 0) {
    drop <- paste0("(^|/)(", paste(exclude, collapse = "|"), ")(/|$)")
    files <- files[!grepl(drop, files)]
  }
  files
}

#' Extract package names referenced in one file
#'
#' Matches library()/require()/requireNamespace()/loadNamespace() calls and
#' pkg:: / pkg::: references. Whole-line comments are ignored; inline comments
#' are not, so a few spurious names may appear.
#'
#' @param file Path to an R or Rmd file.
#' @return Character vector of package names (may contain duplicates).
extract_packages <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  txt <- paste(lines, collapse = "\n")
  
  patterns <- c(
    "(?:library|require|requireNamespace|loadNamespace)\\s*\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)",
    "([A-Za-z][A-Za-z0-9._]*):::?"
  )
  
  unlist(lapply(patterns, function(p) {
    hits <- regmatches(txt, gregexpr(p, txt, perl = TRUE))[[1]]
    if (length(hits) == 0) return(character())
    sub(paste0("^.*?", p, ".*$"), "\\1", hits, perl = TRUE)
  }))
}

#' Is a package installed?
#'
#' Uses system.file() rather than installed.packages(), which caches its result
#' and so can miss packages installed earlier in the same session.
#'
#' @param pkg Package name.
#' @return TRUE or FALSE.
is_installed <- function(pkg) {
  nzchar(system.file(package = pkg))
}

#' Install one package, reporting success without stopping on failure
#'
#' Skips packages that are already present, so the script can be re-run after an
#' interrupted or partial install.
#'
#' @param pkg Package name.
#' @param github Named vector mapping package names to GitHub "user/repo".
#' @return "already installed", "installed", "failed", or "error: <message>".
install_one <- function(pkg, github = github_pkgs) {
  if (is_installed(pkg)) return("already installed")
  
  status <- tryCatch({
    if (pkg %in% names(github)) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = repos)
      }
      remotes::install_github(github[[pkg]], upgrade = "never")
    } else {
      install.packages(pkg, repos = repos)
    }
    NA_character_
  }, error = function(e) paste("error:", conditionMessage(e)))
  
  # install.packages() often warns rather than errors, so verify directly.
  if (requireNamespace(pkg, quietly = TRUE)) {
    "installed"
  } else if (!is.na(status)) {
    status
  } else {
    "failed"
  }
}

## Scan -----------------------------------------------------------------------

files <- find_r_files(repo_dir, exclude_dirs)
message("Scanning ", length(files), " files.")

found <- sort(unique(unlist(lapply(files, extract_packages))))

base_pkgs <- rownames(installed.packages(priority = "base"))
found <- setdiff(found, c(base_pkgs, "R"))

missing <- found[!vapply(found, is_installed, logical(1))]

on_cran <- rownames(available.packages(repos = repos))
unknown <- setdiff(missing, c(on_cran, names(github_pkgs)))

message("Referenced: ", length(found),
        " | already installed: ", length(found) - length(missing),
        " | missing: ", length(missing))
if (length(unknown) > 0) {
  message("Not on CRAN and not in github_pkgs (check these by hand):\n  ",
          paste(unknown, collapse = ", "))
}

## Install --------------------------------------------------------------------

to_install <- setdiff(missing, unknown)

if (dry_run) {
  message("\nDry run. Would attempt to install:\n  ",
          paste(to_install, collapse = ", "))
} else {
  results <- vapply(to_install, function(p) {
    message("\n--- ", p, " ---")
    install_one(p)
  }, character(1))
  
  out <- data.frame(package = to_install, status = unname(results),
                    stringsAsFactors = FALSE)
  ok <- out$status %in% c("installed", "already installed")
  print(out[order(!ok, out$package), ], row.names = FALSE)
  message("\nInstalled: ", sum(ok), " | not installed: ", sum(!ok))
  # write.csv(out, "package_install_results.csv", row.names = FALSE)
  message("\nResults written to package_install_results.csv")
}
