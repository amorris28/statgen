.read_ld_npz_shard <- function(path, check_payload_structure = FALSE, retain_ld_r = TRUE, expected_file_md5 = NULL) {
  .require_ld_file(path)
  payload <- .read_npz_ld_payload(path, expected_file_md5 = expected_file_md5)
  meta <- payload$metadata
  .validate_ld_shard_metadata(meta, path, expected_format = .ld_npz_format)
  .validate_csc_payload_dimensions(path, payload, meta)
  if (isTRUE(check_payload_structure)) {
    .validate_csc_payload_structure(path, payload$data, payload$indices, payload$indptr, meta$num_snp, meta$nnz)
  }
  ld_r <- .csc_payload_to_dgCMatrix(payload$data, payload$indices, payload$indptr, meta$num_snp)
  if (isTRUE(check_payload_structure)) {
    .validate_sparse_matrix_values(path, ld_r)
  }
  shard <- .new_ld_shard(
    chr = meta$chr,
    sex = meta$sex,
    num_snp = meta$num_snp,
    ld_r = ld_r,
    a1freq = payload$a1freq,
    reference_checksum = meta$reference_checksum,
    num_monomorphic_snps = meta$num_monomorphic_snps,
    retain_ld_r = retain_ld_r
  )
  if (shard$num_monomorphic_snps > 0L) {
    warning(
      sprintf(
        "%s: LD shard contains %d monomorphic SNPs; LD involving those SNPs is undefined and represented by omitted off-diagonal entries",
        path, shard$num_monomorphic_snps
      ),
      call. = FALSE
    )
  }
  list(shard = shard, metadata = meta)
}

.csc_payload_to_dgCMatrix <- function(data, indices, indptr, num_snp) {
  .new_dgCMatrix_from_slots(
    i = as.integer(indices),
    p = as.integer(indptr),
    x = as.numeric(data),
    dim = as.integer(c(num_snp, num_snp))
  )
}

.new_dgCMatrix_from_slots <- function(i, p, x, dim, dimnames = list(NULL, NULL)) {
  loadNamespace("Matrix")
  methods::new(
    "dgCMatrix",
    i = as.integer(i),
    p = as.integer(p),
    Dim = as.integer(dim),
    Dimnames = dimnames,
    x = as.numeric(x),
    factors = list()
  )
}

.validate_csc_payload_dimensions <- function(path, payload, meta) {
  .validate_numeric_1d(payload$data, path, "data")
  .validate_integer_1d(payload$indices, path, "indices")
  .validate_integer_1d(payload$indptr, path, "indptr")
  .validate_numeric_1d(payload$a1freq, path, "a1freq")
  shape <- as.integer(payload$shape)
  num_snp <- as.integer(meta$num_snp)
  nnz <- as.integer(meta$nnz)
  if (length(shape) != 2L || any(shape != num_snp)) {
    stop(sprintf("%s: shape must equal [num_snp, num_snp]", path), call. = FALSE)
  }
  if (length(payload$indptr) != num_snp + 1L) {
    stop(sprintf("%s: indptr length must be num_snp + 1", path), call. = FALSE)
  }
  if (length(payload$data) != length(payload$indices) || length(payload$data) != nnz) {
    stop(sprintf("%s: data, indices, and metadata nnz length mismatch", path), call. = FALSE)
  }
  if (as.integer(payload$indptr[[length(payload$indptr)]]) != nnz) {
    stop(sprintf("%s: indptr[-1] must equal metadata nnz", path), call. = FALSE)
  }
  if (length(payload$a1freq) != num_snp) {
    stop(sprintf("%s: a1freq length must equal num_snp", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_csc_payload_structure <- function(path, data, indices, indptr, num_snp, nnz) {
  if (!identical(as.integer(indptr[[1L]]), 0L)) {
    stop(sprintf("%s: indptr[0] must be 0", path), call. = FALSE)
  }
  if (length(indptr) > 1L && any(indptr[-1L] < indptr[-length(indptr)])) {
    stop(sprintf("%s: indptr must be monotonic nondecreasing", path), call. = FALSE)
  }
  if (length(indices) && (min(indices) < 0L || max(indices) >= as.integer(num_snp))) {
    stop(sprintf("%s: sparse row indices out of bounds", path), call. = FALSE)
  }
  if (any(!is.finite(data))) {
    stop(sprintf("%s: sparse data must be finite", path), call. = FALSE)
  }
  if (as.integer(nnz) < as.integer(num_snp)) {
    stop(sprintf("%s: sparse matrix must include explicit unit diagonal", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_sparse_matrix_values <- function(path, mat) {
  diag <- Matrix::diag(mat)
  if (length(diag) != nrow(mat) || any(abs(diag - 1) > 1e-7)) {
    stop(sprintf("%s: ld_r diagonal must be explicit unit", path), call. = FALSE)
  }
  diff <- mat - Matrix::t(mat)
  if (length(diff@x) && max(abs(diff@x)) > 1e-6) {
    stop(sprintf("%s: ld_r must be symmetric", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_numeric_1d <- function(value, path, name) {
  if (!is.numeric(value) || !is.null(dim(value))) {
    stop(sprintf("%s: %s must be a one-dimensional numeric vector", path, name), call. = FALSE)
  }
  invisible(NULL)
}

.validate_integer_1d <- function(value, path, name) {
  if (!is.integer(value) || !is.null(dim(value))) {
    stop(sprintf("%s: %s must be a one-dimensional integer vector", path, name), call. = FALSE)
  }
  invisible(NULL)
}
