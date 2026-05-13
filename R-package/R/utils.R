.canonical_chr_order <- c(as.character(seq_len(22L)), "X")
.chr_rank <- stats::setNames(seq_along(.canonical_chr_order), .canonical_chr_order)
.ignored_chr <- c("Y", "MT")

is_present <- function(x, ...) UseMethod("is_present")

.validate_requested_shards <- function(requested, available_labels, where) {
  available_labels <- as.character(available_labels)
  if (is.null(requested)) {
    if (!length(available_labels)) {
      stop(sprintf("%s: no supported shards are present", where), call. = FALSE)
    }
    return(available_labels)
  }
  if (!is.character(requested) || !length(requested) || anyNA(requested)) {
    stop(sprintf("%s: shards must be a non-empty character vector of unique canonical contig labels", where), call. = FALSE)
  }
  seen <- character()
  prev_rank <- 0L
  for (label in requested) {
    if (!(label %in% names(.chr_rank))) {
      stop(sprintf("%s: unsupported shard label %s; expected canonical labels 1-22 or X", where, sQuote(label)), call. = FALSE)
    }
    if (label %in% seen) {
      stop(sprintf("%s: duplicate shard label %s in shards", where, sQuote(label)), call. = FALSE)
    }
    rank <- .chr_rank[[label]]
    if (rank <= prev_rank) {
      stop(sprintf("%s: shards must be in canonical subsequence order", where), call. = FALSE)
    }
    if (!(label %in% available_labels)) {
      stop(sprintf("%s: requested shard %s is not present", where, sQuote(label)), call. = FALSE)
    }
    seen <- c(seen, label)
    prev_rank <- rank
  }
  requested
}

.validate_char_scalar <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be a character scalar", name), call. = FALSE)
  }
  value
}

.validate_path_scalar <- function(path, name) {
  .validate_char_scalar(path, name)
}

.validate_equal_lengths <- function(values, message) {
  lengths <- vapply(values, length, integer(1))
  if (length(unique(lengths)) != 1L) {
    stop(message, call. = FALSE)
  }
}

.validate_cache_payload_metadata <- function(payload, schema, cache_name, extra_fields = character()) {
  if (!is.list(payload) || is.null(payload$metadata)) {
    stop(sprintf("Invalid %s cache: expected an RDS list with metadata", cache_name), call. = FALSE)
  }
  meta <- payload$metadata
  if (!identical(meta$schema, schema)) {
    stop(sprintf("Unsupported %s cache schema: %s", cache_name, sQuote(as.character(meta$schema))), call. = FALSE)
  }
  n_shards <- suppressWarnings(as.integer(meta$n_shards))
  if (length(n_shards) != 1L || is.na(n_shards) || n_shards < 1L) {
    stop(sprintf("Invalid %s cache: n_shards must be at least 1", cache_name), call. = FALSE)
  }
  for (field in c("shard_labels", "shard_checksums", "shard_start0", "shard_stop0", extra_fields)) {
    if (is.null(meta[[field]]) || length(meta[[field]]) != n_shards) {
      stop(sprintf("Invalid %s cache: metadata.%s length mismatch", cache_name, field), call. = FALSE)
    }
  }
  total <- meta$shard_stop0[[n_shards]]
  if (!identical(as.integer(meta$shard_start0[[1]]), 0L) || any(meta$shard_stop0 < meta$shard_start0)) {
    stop(sprintf("Invalid %s cache: shard offsets are invalid", cache_name), call. = FALSE)
  }
  if (n_shards > 1L && any(meta$shard_start0[-1L] != meta$shard_stop0[-n_shards])) {
    stop(sprintf("Invalid %s cache: shard offsets are not contiguous", cache_name), call. = FALSE)
  }
  list(meta = meta, n_shards = n_shards, total = total)
}
