.phase2_reference_template <- function() {
  ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen", mustWork = TRUE)
  file.path(dirname(ref_path), "reference_chr@.bim")
}

.phase2_capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

test_that("load_sumstats aligns TSV rows and cache round-trips", {
  sum_path <- system.file("extdata", "traits.tsv.gz", package = "statgen", mustWork = TRUE)
  ref <- load_reference(.phase2_reference_template())

  expect_warning(s <- load_sumstats(sum_path, ref), "zvec has missing values")
  expect_s3_class(s, "Sumstats")
  expect_equal(num_snp(s), 8L)
  expect_equal(is_present(s), c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(zvec(s), c(2.5, NaN, 1.8, -1.2, NaN, 3.0, 0.5, NaN))
  expect_equal(nvec(s), c(1000, 1000, 1000, 1000, NaN, 500, 500, NaN))
  expect_true(is.infinite(logpvec(s)[[3]]))
  expect_true(is.null(beta_vec(s)))

  cache <- tempfile(fileext = ".rds")
  save_sumstats_cache(s, cache)
  expect_warning(loaded <- load_sumstats_cache(cache), "zvec has missing values")
  expect_equal(logpvec(loaded), logpvec(s))
  expect_equal(nvec(loaded), nvec(s))
  expect_equal(is_present(loaded), is_present(s))

  x_only <- select_shards(loaded, "X")
  expect_equal(num_snp(x_only), 3L)
  expect_equal(is_present(x_only), c(TRUE, TRUE, FALSE))
})

test_that("load_sumstats accepts documented column aliases", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tN",
    "1\t100\trs1001\tA\tG\t0\t2\t100"
  ), path)

  captured <- .phase2_capture_warnings(load_sumstats(path, ref))
  s <- captured$value
  expect_true(is.infinite(logpvec(s)[[1]]))
  expect_equal(zvec(s)[[1]], 2)
  expect_equal(nvec(s)[[1]], 100)
  expect_true(is_present(s)[[1]])
})

test_that("create_sumstats preserves optional fields and save_cache convenience", {
  ref <- load_reference(.phase2_reference_template())
  n <- num_snp(ref)
  s <- create_sumstats(
    ref,
    p = rep(0.5, n),
    z = seq_len(n),
    n = rep(100, n),
    beta = rep(0.1, n),
    se = rep(0.2, n),
    eaf = rep(0.3, n),
    info = rep(0.9, n)
  )

  expect_equal(nvec(s), rep(100, n))
  expect_equal(beta_vec(s), rep(0.1, n))
  expect_equal(se_vec(s), rep(0.2, n))
  expect_equal(eaf_vec(s), rep(0.3, n))
  expect_equal(info_vec(s), rep(0.9, n))

  cache <- tempfile(fileext = ".rds")
  save_cache(s, cache)
  loaded <- load_sumstats_cache(cache, shards = "X")
  expect_equal(num_snp(loaded), 3L)
  expect_equal(nvec(loaded), rep(100, 3L))
})

test_that("load_sumstats warns for absent zvec and nvec", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chr\tbp\ta1\ta2\tp",
    "1\t100\tA\tG\t0.1"
  ), path)

  captured <- .phase2_capture_warnings(load_sumstats(path, ref))
  expect_true(any(grepl("zvec is absent", captured$warnings)))
  expect_true(any(grepl("nvec is absent", captured$warnings)))
  expect_null(zvec(captured$value))
  expect_null(nvec(captured$value))

  cache <- tempfile(fileext = ".rds")
  save_sumstats_cache(captured$value, cache)
  loaded <- .phase2_capture_warnings(load_sumstats_cache(cache))
  expect_true(any(grepl("zvec is absent", loaded$warnings)))
  expect_true(any(grepl("nvec is absent", loaded$warnings)))
})

test_that("load_annotations paints BED intervals and cache round-trips", {
  bed1 <- system.file("extdata", "anno1.bed", package = "statgen", mustWork = TRUE)
  bed2 <- system.file("extdata", "anno2.bed", package = "statgen", mustWork = TRUE)
  ref <- load_reference(.phase2_reference_template())
  ann <- load_annotations(c(bed1, bed2), ref)

  expect_s3_class(ann, "AnnotationPanel")
  expect_equal(num_snp(ann), 8L)
  expect_equal(num_annot(ann), 2L)
  expect_equal(annonames(ann), c("anno1", "anno2"))
  expect_s4_class(annomat(ann), "lgCMatrix")
  ann_dense <- as.matrix(annomat(ann))
  expect_equal(colnames(annomat(ann)), c("anno1", "anno2"))
  expect_equal(as.vector(ann_dense[, 1L]), c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(as.vector(ann_dense[, 2L]), c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))

  cache <- tempfile(fileext = ".rds")
  save_cache(ann, cache)
  loaded <- load_annotations_cache(cache, shards = "X")
  expect_equal(annonames(loaded), c("anno1", "anno2"))
  expect_equal(colnames(annomat(loaded)), c("anno1", "anno2"))
  expect_equal(as.matrix(annomat(loaded)), as.matrix(annomat(ann))[6:8, ])
})

test_that("load_annotations accepts scalar path and leading comments", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  bed <- tempfile(fileext = ".bed")
  writeLines(c("# leading comment", "", "1\t99\t200"), bed)

  ann <- load_annotations(bed, ref)
  expect_equal(num_annot(ann), 1L)
  expect_equal(as.vector(as.matrix(annomat(ann))), c(TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("annotation factories, selection, and union work", {
  ref <- load_reference(.phase2_reference_template())
  n <- num_snp(ref)
  a <- create_annotation(ref, rep(c(1, 0), length.out = n), "a")
  b <- create_annotations(
    ref,
    annotation_matrix = matrix(rep(c(0, 1), length.out = n), ncol = 1L),
    annotation_names = "b"
  )
  u <- union_annotations(a, b)

  expect_equal(annonames(u), c("a", "b"))
  expect_equal(colnames(annomat(u)), c("a", "b"))
  selected <- select_annotations(u, names = c("b", "a"))
  expect_equal(annonames(selected), c("b", "a"))
  expect_equal(colnames(annomat(selected)), c("b", "a"))
  expect_error(select_annotations(u, "missing"), "unknown annotation")
})

test_that("annotations cache validation rejects bad schema and zero shards", {
  ref <- load_reference(.phase2_reference_template(), shards = "1")
  ann <- create_annotation(ref, rep(TRUE, num_snp(ref)), "a")
  cache <- tempfile(fileext = ".rds")
  save_annotations_cache(ann, cache)
  payload <- readRDS(cache)

  bad_schema <- tempfile(fileext = ".rds")
  payload$metadata$schema <- "annotations_cache/bad"
  saveRDS(payload, bad_schema)
  expect_error(load_annotations_cache(bad_schema), "Unsupported annotations cache schema")

  zero_shards <- tempfile(fileext = ".rds")
  saveRDS(
    list(
      metadata = list(
        schema = "annotations_cache/0.1",
        n_shards = 0L,
        shard_labels = character(),
        shard_checksums = character(),
        shard_start0 = numeric(),
        shard_stop0 = numeric()
      ),
      annomat = Matrix::Matrix(matrix(logical(), nrow = 0L, ncol = 1L), sparse = TRUE),
      annonames = "a"
    ),
    zero_shards
  )
  expect_error(load_annotations_cache(zero_shards), "n_shards must be at least 1")
})
