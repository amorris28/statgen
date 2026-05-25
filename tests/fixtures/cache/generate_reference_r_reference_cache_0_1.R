#!/usr/bin/env Rscript
# Generate tests/fixtures/cache/reference_r_reference_cache_0_1.rds.
# Run from the repository root:
#   Rscript tests/fixtures/cache/generate_reference_r_reference_cache_0_1.R

statgen_cache_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
source(file.path(statgen_cache_dir, "cache_fixture_helpers.R"))
statgen_generate_r_cache_fixture("reference")
