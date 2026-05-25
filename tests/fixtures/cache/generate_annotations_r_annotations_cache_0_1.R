#!/usr/bin/env Rscript
# Generate tests/fixtures/cache/annotations_r_annotations_cache_0_1.rds.
#
# Run from a v0.3.2 checkout at the repository root:
#   Rscript tests/fixtures/cache/generate_annotations_r_annotations_cache_0_1.R
#
# This uses the official v0.3.2 R implementation, whose annotation writer
# produces annotations_cache/0.1.

args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1]])
cache_dir <- dirname(normalizePath(file_arg))
repo_root <- normalizePath(file.path(cache_dir, "..", "..", ".."), mustWork = TRUE)

r_files <- file.path(
  repo_root,
  "R-package",
  "R",
  c("utils.R", "verbosity.R", "hash.R", "variant_match.R", "bfile_utils.R", "reference.R", "annotations.R")
)
for (path in r_files) {
  source(path)
}

reference <- load_reference(file.path("tests", "fixtures", "reference", "sharded", "@.bim"))
annotations <- load_annotations(
  file.path("tests", "fixtures", "annotations", c("anno1.bed", "anno2.bed")),
  reference
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
save_annotations_cache(
  annotations,
  file.path(cache_dir, "annotations_r_annotations_cache_0_1.rds")
)
