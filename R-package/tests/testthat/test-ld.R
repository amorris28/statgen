copy_ld_fixture <- function() {
  src <- system.file("extdata", "ld", "python", package = "statgen", mustWork = TRUE)
  dst <- tempfile("ld_npz_")
  dir.create(dst)
  ok <- file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE)
  stopifnot(all(ok))
  dst
}

test_that("prepare_ld_npz_for_r enables LD loading and operations", {
  old_verbosity <- get_verbosity()
  on.exit(set_verbosity(old_verbosity), add = TRUE)
  set_verbosity("quiet")

  ld_root <- copy_ld_fixture()

  expect_error(load_ld(ld_root), "missing r_reference_cache")

  manifest <- prepare_ld_npz_for_r(ld_root, extract_npz = TRUE)
  expect_identical(manifest$runtime_format, "python_npz_csc32")
  expect_identical(manifest$r_reference_cache, "reference_cache.rds")
  expect_match(manifest$r_reference_cache_md5, "^[0-9a-f]{32}$")
  expect_true(file.exists(file.path(ld_root, manifest$r_reference_cache)))
  cache_dir <- file.path(ld_root, paste0(manifest$shards[[1]]$file, ".d"))
  expect_true(dir.exists(cache_dir))
  expect_true(file.exists(file.path(cache_dir, "data.npy")))
  marker <- jsonlite::fromJSON(file.path(cache_dir, ".statgen-extracted-npz.json"))
  expect_identical(marker$source_file, manifest$shards[[1]]$file)
  expect_identical(marker$source_file_md5, manifest$shards[[1]]$file_md5)

  manifest2 <- prepare_ld_npz_for_r(ld_root)
  expect_identical(manifest2$r_reference_cache, manifest$r_reference_cache)
  expect_match(manifest2$r_reference_cache_md5, "^[0-9a-f]{32}$")

  report <- validate_ld_distribution(ld_root, check_payload_structure = TRUE)
  expect_true(report$ok)

  ref <- load_ld_reference(ld_root)
  ld <- load_ld(ld_root)
  expect_s3_class(ld, "LDPanel")
  expect_equal(num_snp(ld), 8L)
  expect_equal(num_snp(ref), 8L)
  expect_equal(vapply(shards(ld), function(s) s$label, character(1)), c("1", "X"))
  expect_equal(default_chrX_sex(ld), "female")

  expect_equal(
    round(a1freq(ld), 6),
    c(0.30, 0.40, 0.20, 0.35, 0.15, 0.28, 0.32, 0.22)
  )
  expect_equal(
    round(a1freq(ld, chrX_sex = "combined"), 6),
    c(0.30, 0.40, 0.20, 0.35, 0.15, 0.25, 0.30, 0.20)
  )

  expect_equal(
    round(multiply_r2(ld, seq_len(num_snp(ld))), 6),
    c(2.89, 3.81, 3.09, 6.95, 6.96, 8.9575, 10.515, 8.8575)
  )
  expect_equal(
    round(fast_prune(c(9, 9, 7, 1, 2, 6, 5, 4), ld, 0.35), 6),
    c(9, NaN, 7, NaN, 2, 6, NaN, 4)
  )

  ld_no_r <- load_ld(ld_root, retain_ld_r = FALSE)
  expect_true(is.null(ld_no_r$shard_groups[[1]][[1]]$ld_r))
  expect_equal(round(multiply_r2(ld_no_r, seq_len(num_snp(ld_no_r)))[[1]], 6), 2.89)
})
