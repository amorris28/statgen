test_that("load_reference reads bundled BIM fixture", {
  path <- system.file("extdata", "reference_chr1.bim", package = "statgen", mustWork = TRUE)
  ref <- load_reference(path)

  expect_s3_class(ref, "ReferencePanel")
  expect_equal(num_snp(ref), 5L)
  expect_equal(vapply(shards(ref), function(s) s$label, character(1)), "1")
  expect_equal(bp(ref), c(100L, 200L, 300L, 400L, 500L))
  expect_equal(snp(ref)[1:2], c("rs1001", "rs1002"))
  expect_equal(a1(ref)[1], "A")
  expect_equal(a2(ref)[1], "G")
  expect_true(inherits(a1_hash64(ref), "integer64"))
  expect_true(inherits(a2_hash64(ref), "integer64"))
})

test_that("reference cache round-trips and subsets", {
  path <- system.file("extdata", "reference_chr1.bim", package = "statgen", mustWork = TRUE)
  ref <- load_reference(path)
  cache <- tempfile(fileext = ".rds")

  save_reference_cache(ref, cache)
  loaded <- load_reference_cache(cache)

  expect_equal(num_snp(loaded), num_snp(ref))
  expect_equal(chr(loaded), chr(ref))
  expect_equal(snp(loaded), snp(ref))
  expect_equal(bp(loaded), bp(ref))
  expect_equal(a1(loaded), a1(ref))
  expect_equal(a2(loaded), a2(ref))
  expect_equal(as.character(a1_hash64(loaded)), as.character(a1_hash64(ref)))
  expect_true(validate_checksums(loaded))

  subset <- load_reference_cache(cache, shards = "1")
  expect_equal(num_snp(subset), 5L)
})
