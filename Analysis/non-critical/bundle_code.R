# Bundle a repo's code files into a single text file, so that can have ai read
# for context if needed
#
# Reads directly from the original repo (no intermediate copy needed).
# Prints a file tree at the top, then concatenates each file with a header
# showing its relative path. The original folder is never modified.

library(fs)
library(stringr)
source('Functions/init.R')
# ---- settings ---------------------------------------------------------------

# Folder to bundle (the original repo).
src <- paths$large0

# Output file.
out <- file.path(paths$large, 'Data_processed', 'PED_code_bundle.txt')

# Extensions to include (case-insensitive).
keep_ext <- c("R", "Rmd", "qmd")

# Folders to exclude, as paths RELATIVE to src. May be nested.
# Matches the folder and everything under it. Use forward slashes.
# e.g. "Analysis/viz" excludes src/Analysis/viz/ and all its contents.
exclude_dirs <- c(
  ".git",
  ".Rproj.user",
  "renv",
  ".quarto",
  "Analysis/BiomassPhenology",
  "Analysis/BiomassQuantity",
  "Analysis/VegComposition/Visualizations"
)

# Optional: split output if it exceeds this many MB. Inf = always one file.
max_mb <- Inf

# ---- gather files -----------------------------------------------------------

src <- path_norm(src)

files <- dir_ls(src, recurse = TRUE, type = "file")
files <- files[str_to_lower(path_ext(files)) %in% str_to_lower(keep_ext)]

rel <- path_rel(files, start = src)

#' Is a relative path inside any excluded directory?
#' @param rel_paths Character vector of paths relative to src.
#' @param ex Character vector of excluded dirs (relative, forward slashes).
#' @return Logical vector, TRUE if the path should be excluded.
is_excluded <- function(rel_paths, ex) {
  if (length(ex) == 0) return(rep(FALSE, length(rel_paths)))
  ex <- str_replace(ex, "/+$", "")  # drop any trailing slash
  out <- rep(FALSE, length(rel_paths))
  for (e in ex) {
    # exclude the dir itself or anything beneath it
    out <- out | rel_paths == e | str_starts(rel_paths, str_c(e, "/"))
  }
  out
}

keep <- !is_excluded(rel, exclude_dirs)
files <- files[keep]
rel   <- rel[keep]

if (length(files) == 0) stop("No matching files under: ", src)

# stable, grouped-by-folder order
ord   <- order(rel)
files <- files[ord]
rel   <- rel[ord]

# ---- build a simple indented file tree --------------------------------------

#' Make an indented text tree from relative paths
#' @param rel_paths Character vector of relative file paths.
#' @return Character vector, one line per directory/file.
make_tree <- function(rel_paths) {
  parts <- str_split(rel_paths, "/")
  seen_dirs <- character(0)
  lines <- character(0)
  for (p in parts) {
    if (length(p) > 1) {
      for (d in seq_len(length(p) - 1)) {
        dir_key <- str_c(p[seq_len(d)], collapse = "/")
        if (!dir_key %in% seen_dirs) {
          seen_dirs <- c(seen_dirs, dir_key)
          lines <- c(lines, str_c(strrep("  ", d - 1), p[d], "/"))
        }
      }
    }
    lines <- c(lines, str_c(strrep("  ", length(p) - 1), p[length(p)]))
  }
  lines
}

tree_lines <- make_tree(rel)

# ---- assemble the bundle ----------------------------------------------------

header <- c(
  "================================================================",
  str_c("PROJECT CODE BUNDLE: ", path_file(src)),
  str_c("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  str_c("Files included: ", length(files), "  |  Extensions: ",
        str_c(keep_ext, collapse = ", ")),
  if (length(exclude_dirs)) str_c("Excluded folders: ",
                                  str_c(exclude_dirs, collapse = ", ")),
  "================================================================",
  "",
  "FILE TREE",
  "----------------------------------------------------------------",
  tree_lines,
  "",
  "================================================================",
  "FILE CONTENTS",
  "================================================================",
  ""
)

#' Build the text block for one file
file_block <- function(f, r) {
  body <- tryCatch(
    readLines(f, warn = FALSE, encoding = "UTF-8"),
    error = function(e) str_c("[could not read file: ", conditionMessage(e), "]")
  )
  c(
    "",
    "================================================================",
    str_c("FILE: ", r),
    "================================================================",
    body,
    ""
  )
}

blocks <- mapply(file_block, files, rel, SIMPLIFY = FALSE)

# ---- write out (single file, or split by size) ------------------------------

write_lines_to <- function(lines, file) {
  writeLines(lines, file, useBytes = TRUE)
  message("Wrote ", file, " (",
          round(file_info(file)$size / 1024^2, 2), " MB)")
}

if (is.infinite(max_mb)) {
  write_lines_to(c(header, unlist(blocks)), out)
} else {
  limit_bytes <- max_mb * 1024^2
  part <- 1L
  buf <- header
  buf_bytes <- sum(nchar(buf, type = "bytes")) + length(buf)
  flush <- function() {
    f <- str_replace(out, "\\.txt$", str_c("_part", part, ".txt"))
    write_lines_to(buf, f)
  }
  for (b in blocks) {
    b_bytes <- sum(nchar(b, type = "bytes")) + length(b)
    if (buf_bytes + b_bytes > limit_bytes && length(buf) > length(header)) {
      flush(); part <- part + 1L
      buf <- header; buf_bytes <- sum(nchar(buf, type = "bytes")) + length(buf)
    }
    buf <- c(buf, b); buf_bytes <- buf_bytes + b_bytes
  }
  flush()
}
