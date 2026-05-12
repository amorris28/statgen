.canonical_chr_order <- c(as.character(seq_len(22L)), "X")
.chr_rank <- stats::setNames(seq_along(.canonical_chr_order), .canonical_chr_order)
.ignored_chr <- c("Y", "MT")

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

.validate_path_scalar <- function(path, name) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop(sprintf("%s must be a character scalar", name), call. = FALSE)
  }
  path
}

.validate_equal_lengths <- function(values, message) {
  lengths <- vapply(values, length, integer(1))
  if (length(unique(lengths)) != 1L) {
    stop(message, call. = FALSE)
  }
}
