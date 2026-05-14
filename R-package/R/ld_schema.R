.ld_manifest_schema <- "1.0"
.ld_runtime_format <- "python_npz_csc32"
.ld_npz_format <- "statgen_ld_npz_csc32"
.ld_valid_chrx_sex <- c("female", "male", "combined")
.ld_required_arrays <- c("data", "indices", "indptr", "shape", "a1freq", "metadata")

.read_ld_manifest <- function(path, require_r_reference = FALSE) {
  .require_ld_file(path)
  manifest <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!identical(manifest$object_type, "ld_panel_manifest")) {
    stop(sprintf("%s: object_type must be 'ld_panel_manifest'", path), call. = FALSE)
  }
  if (!identical(manifest$schema_version, .ld_manifest_schema)) {
    stop(sprintf("%s: unsupported LD manifest schema_version %s", path, sQuote(as.character(manifest$schema_version))), call. = FALSE)
  }
  if (!identical(manifest$runtime_format, .ld_runtime_format)) {
    stop(sprintf("%s: unsupported LD runtime_format %s", path, sQuote(as.character(manifest$runtime_format))), call. = FALSE)
  }
  .validate_plain_relative_filename(manifest$reference_cache, sprintf("%s: reference_cache", path))
  .validate_md5_hex(manifest$reference_cache_md5, sprintf("%s: reference_cache_md5", path))
  has_r_reference <- !is.null(manifest$r_reference_cache)
  has_r_md5 <- !is.null(manifest$r_reference_cache_md5)
  if (xor(has_r_reference, has_r_md5)) {
    stop(sprintf("%s: r_reference_cache and r_reference_cache_md5 must be present together", path), call. = FALSE)
  }
  if (isTRUE(require_r_reference) && !has_r_reference) {
    stop(sprintf("%s: missing r_reference_cache/r_reference_cache_md5 for R LD loading", path), call. = FALSE)
  }
  if (has_r_reference) {
    .validate_plain_relative_filename(manifest$r_reference_cache, sprintf("%s: r_reference_cache", path))
    .validate_md5_hex(manifest$r_reference_cache_md5, sprintf("%s: r_reference_cache_md5", path))
  }
  if (!is.list(manifest$shards) || !length(manifest$shards)) {
    stop(sprintf("%s: shards must be a non-empty list", path), call. = FALSE)
  }
  seen <- character()
  reference_identity_by_chr <- list()
  for (i in seq_along(manifest$shards)) {
    entry <- manifest$shards[[i]]
    .validate_ld_manifest_entry(entry, sprintf("%s:shards[%d]", path, i - 1L))
    key <- paste(entry$chr, if (is.null(entry$sex)) "<null>" else entry$sex, sep = "\r")
    if (key %in% seen) {
      stop(sprintf("%s: duplicate manifest entry for chr %s sex %s", path, entry$chr, sQuote(as.character(entry$sex))), call. = FALSE)
    }
    seen <- c(seen, key)
    identity <- list(
      reference_bim = entry$reference_bim,
      num_snp = as.integer(entry$num_snp),
      reference_checksum = entry$reference_checksum
    )
    prior <- reference_identity_by_chr[[entry$chr]]
    if (is.null(prior)) {
      reference_identity_by_chr[[entry$chr]] <- identity
    } else if (!identical(prior, identity)) {
      stop(
        sprintf(
          "%s: manifest entries for chr %s must share reference_bim, num_snp, and reference_checksum",
          path, entry$chr
        ),
        call. = FALSE
      )
    }
  }
  manifest
}

.validate_ld_manifest_entry <- function(entry, where) {
  required <- c("chr", "sex", "file", "file_md5", "num_snp", "nnz", "reference_checksum", "reference_bim")
  missing <- setdiff(required, names(entry))
  if (length(missing)) {
    stop(sprintf("%s: missing required fields: %s", where, paste(missing, collapse = ", ")), call. = FALSE)
  }
  .validate_chr_sex(entry$chr, entry$sex, where)
  .validate_relative_path(entry$file, sprintf("%s: file", where))
  .validate_md5_hex(entry$file_md5, sprintf("%s: file_md5", where))
  .validate_positive_int(entry$num_snp, sprintf("%s: num_snp", where))
  .validate_nonnegative_int(entry$nnz, sprintf("%s: nnz", where))
  if (!is.character(entry$reference_checksum) || length(entry$reference_checksum) != 1L || identical(entry$reference_checksum, "")) {
    stop(sprintf("%s: reference_checksum must be a non-empty string", where), call. = FALSE)
  }
  .validate_plain_relative_filename(entry$reference_bim, sprintf("%s: reference_bim", where))
  invisible(NULL)
}

.validate_ld_shard_metadata <- function(meta, path, expected_format) {
  required <- c(
    "object_type", "schema_version", "format", "chr", "sex", "num_snp", "nnz",
    "matrix", "diagonal", "value", "reference_checksum", "reference_bim",
    "num_monomorphic_snps"
  )
  missing <- setdiff(required, names(meta))
  if (length(missing)) {
    stop(sprintf("%s: metadata missing required fields: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!identical(meta$object_type, "ld_shard")) {
    stop(sprintf("%s: metadata object_type must be 'ld_shard'", path), call. = FALSE)
  }
  if (!identical(meta$schema_version, .ld_manifest_schema)) {
    stop(sprintf("%s: unsupported metadata schema_version %s", path, sQuote(as.character(meta$schema_version))), call. = FALSE)
  }
  if (!identical(meta$format, expected_format)) {
    stop(sprintf("%s: metadata format must be %s", path, sQuote(expected_format)), call. = FALSE)
  }
  .validate_chr_sex(meta$chr, meta$sex, sprintf("%s: metadata", path))
  .validate_positive_int(meta$num_snp, sprintf("%s: metadata num_snp", path))
  .validate_nonnegative_int(meta$nnz, sprintf("%s: metadata nnz", path))
  .validate_nonnegative_int(meta$num_monomorphic_snps, sprintf("%s: metadata num_monomorphic_snps", path))
  if (as.integer(meta$num_monomorphic_snps) > as.integer(meta$num_snp)) {
    stop(sprintf("%s: metadata num_monomorphic_snps must not exceed num_snp", path), call. = FALSE)
  }
  if (as.numeric(meta$nnz) >= 2^31) {
    stop(sprintf("%s: CSC32 metadata nnz must be < 2^31", path), call. = FALSE)
  }
  if (!identical(meta$matrix, "symmetric") || !identical(meta$diagonal, "explicit_unit") || !identical(meta$value, "r")) {
    stop(sprintf("%s: metadata matrix/diagonal/value fields are invalid", path), call. = FALSE)
  }
  if (!is.character(meta$reference_checksum) || length(meta$reference_checksum) != 1L || identical(meta$reference_checksum, "")) {
    stop(sprintf("%s: metadata reference_checksum must be a non-empty string", path), call. = FALSE)
  }
  .validate_plain_relative_filename(meta$reference_bim, sprintf("%s: metadata reference_bim", path))
  if (!is.null(meta$sparse_layout) && !identical(meta$sparse_layout, "csc")) {
    stop(sprintf("%s: metadata sparse_layout must be 'csc'", path), call. = FALSE)
  }
  if (!is.null(meta$index_base) && !identical(as.integer(meta$index_base), 0L)) {
    stop(sprintf("%s: metadata index_base must be 0", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_manifest_entry_agreement <- function(entry, meta, path) {
  for (key in c("chr", "sex", "num_snp", "nnz", "reference_checksum", "reference_bim")) {
    if (!identical(entry[[key]], meta[[key]])) {
      stop(sprintf("%s: manifest/per-file metadata mismatch for %s", path, key), call. = FALSE)
    }
  }
  invisible(NULL)
}

.validate_chr_sex <- function(chr_label, sex, where) {
  if (!is.character(chr_label) || length(chr_label) != 1L || !(chr_label %in% .canonical_chr_order)) {
    stop(sprintf("%s: chr must be one of 1-22 or X", where), call. = FALSE)
  }
  if (identical(chr_label, "X")) {
    .validate_chrx_sex(sex, sprintf("%s: sex", where))
  } else if (!is.null(sex)) {
    stop(sprintf("%s: autosomal LD shards must have sex null", where), call. = FALSE)
  }
  invisible(NULL)
}

.validate_chrx_sex <- function(sex, where) {
  if (!is.character(sex) || length(sex) != 1L || is.na(sex) || !(sex %in% .ld_valid_chrx_sex)) {
    stop(sprintf("%s: chrX sex must be one of female, male, combined", where), call. = FALSE)
  }
  sex
}

.validate_positive_int <- function(value, where) {
  value_num <- suppressWarnings(as.numeric(value))
  if (length(value_num) != 1L || is.na(value_num) || !is.finite(value_num) || floor(value_num) != value_num || value_num <= 0) {
    stop(sprintf("%s must be a positive integer", where), call. = FALSE)
  }
  invisible(NULL)
}

.validate_nonnegative_int <- function(value, where) {
  value_num <- suppressWarnings(as.numeric(value))
  if (length(value_num) != 1L || is.na(value_num) || !is.finite(value_num) || floor(value_num) != value_num || value_num < 0) {
    stop(sprintf("%s must be a non-negative integer", where), call. = FALSE)
  }
  invisible(NULL)
}

.validate_plain_relative_filename <- function(value, where) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || identical(value, "")) {
    stop(sprintf("%s must be a non-empty filename", where), call. = FALSE)
  }
  if (grepl("[/\\\\]", value) || grepl("(^|[.])[.]", value) || value %in% c(".", "..")) {
    stop(sprintf("%s must be a plain relative filename", where), call. = FALSE)
  }
  value
}

.validate_relative_path <- function(value, where) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || identical(value, "")) {
    stop(sprintf("%s must be a non-empty relative path", where), call. = FALSE)
  }
  if (grepl("^(/|[A-Za-z]:)", value) || any(strsplit(value, "[/\\\\]")[[1L]] == "..")) {
    stop(sprintf("%s must be a relative path inside the LD distribution", where), call. = FALSE)
  }
  value
}

.validate_md5_hex <- function(value, where) {
  if (!is.character(value) || length(value) != 1L || !grepl("^[0-9a-f]{32}$", value)) {
    stop(sprintf("%s must be a lowercase MD5 hex string", where), call. = FALSE)
  }
  invisible(NULL)
}

.validate_ld_root <- function(path, where) {
  root <- .validate_path_scalar(path, "path")
  if (!dir.exists(root)) {
    stop(sprintf("%s: path must identify a panel root directory", where), call. = FALSE)
  }
  root
}

.require_ld_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop(sprintf("LD file not found: %s", path), call. = FALSE)
  }
  invisible(NULL)
}

.md5_file <- function(path) {
  digest::digest(path, algo = "md5", file = TRUE)
}

.write_ld_manifest <- function(manifest, path) {
  text <- jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null", pretty = TRUE)
  tmp <- tempfile(".ld_manifest_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(text, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) {
    stop(sprintf("Failed to replace LD manifest: %s", path), call. = FALSE)
  }
  invisible(NULL)
}
