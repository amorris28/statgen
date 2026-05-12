.check_unique_variant_keys <- function(bp, a1_hash64, a2_hash64, where, object_name = "variant") {
  .validate_equal_lengths(
    list(bp = bp, a1_hash64 = a1_hash64, a2_hash64 = a2_hash64),
    "numeric variant key vector lengths mismatch"
  )
  if (length(bp) < 2L) {
    return(invisible(NULL))
  }
  ord <- order(bp, a1_hash64, a2_hash64)
  duplicate <- bp[ord][-1L] == bp[ord][-length(ord)] &
    a1_hash64[ord][-1L] == a1_hash64[ord][-length(ord)] &
    a2_hash64[ord][-1L] == a2_hash64[ord][-length(ord)]
  if (any(duplicate)) {
    stop(sprintf("%s: duplicate %s matching key", where, object_name), call. = FALSE)
  }
  invisible(NULL)
}

.match_variant_keys <- function(ref_bp, ref_a1_hash64, ref_a2_hash64,
                                src_bp, src_a1_hash64, src_a2_hash64,
                                shard_label, source_name) {
  ref_a1_hash64 <- bit64::as.integer64(ref_a1_hash64)
  ref_a2_hash64 <- bit64::as.integer64(ref_a2_hash64)
  src_a1_hash64 <- bit64::as.integer64(src_a1_hash64)
  src_a2_hash64 <- bit64::as.integer64(src_a2_hash64)

  .check_unique_variant_keys(
    ref_bp, ref_a1_hash64, ref_a2_hash64,
    sprintf("reference shard %s", shard_label),
    "reference"
  )
  .check_unique_variant_keys(
    src_bp, src_a1_hash64, src_a2_hash64,
    sprintf("source shard %s", shard_label),
    source_name
  )

  match_src <- rep.int(NA_integer_, length(ref_bp))
  if (!length(ref_bp) || !length(src_bp)) {
    return(match_src)
  }

  pairs <- .variant_key_pair_indices(
    ref_bp, ref_a1_hash64, ref_a2_hash64,
    src_bp, src_a1_hash64, src_a2_hash64
  )
  if (length(pairs$ref_index)) {
    match_src[pairs$ref_index] <- pairs$src_index
  }
  match_src
}

.count_variant_key_intersections <- function(ref_bp, ref_a1_hash64, ref_a2_hash64,
                                             src_bp, src_a1_hash64, src_a2_hash64) {
  ref_a1_hash64 <- bit64::as.integer64(ref_a1_hash64)
  ref_a2_hash64 <- bit64::as.integer64(ref_a2_hash64)
  src_a1_hash64 <- bit64::as.integer64(src_a1_hash64)
  src_a2_hash64 <- bit64::as.integer64(src_a2_hash64)

  .check_unique_variant_keys(
    ref_bp, ref_a1_hash64, ref_a2_hash64,
    "reference",
    "reference"
  )
  .check_unique_variant_keys(
    src_bp, src_a1_hash64, src_a2_hash64,
    "source",
    "variant"
  )

  if (!length(ref_bp) || !length(src_bp)) {
    return(0L)
  }
  length(.variant_key_pair_indices(
    ref_bp, ref_a1_hash64, ref_a2_hash64,
    src_bp, src_a1_hash64, src_a2_hash64
  )$ref_index)
}

.variant_key_pair_indices <- function(ref_bp, ref_a1_hash64, ref_a2_hash64,
                                      src_bp, src_a1_hash64, src_a2_hash64) {
  n_ref <- length(ref_bp)
  n_src <- length(src_bp)
  if (!n_ref || !n_src) {
    return(list(ref_index = integer(), src_index = integer()))
  }

  bp_all <- c(as.integer(ref_bp), as.integer(src_bp))
  a1_all <- c(bit64::as.integer64(ref_a1_hash64), bit64::as.integer64(src_a1_hash64))
  a2_all <- c(bit64::as.integer64(ref_a2_hash64), bit64::as.integer64(src_a2_hash64))
  is_source <- c(rep.int(FALSE, n_ref), rep.int(TRUE, n_src))
  local_index <- c(seq_len(n_ref), seq_len(n_src))

  # Sort one combined table and find adjacent equal tuple keys. This keeps the
  # allele hashes as bit64 integer64 values throughout and avoids string keys.
  ord <- order(bp_all, a1_all, a2_all)
  if (length(ord) < 2L) {
    return(list(ref_index = integer(), src_index = integer()))
  }

  bp_sorted <- bp_all[ord]
  a1_sorted <- a1_all[ord]
  a2_sorted <- a2_all[ord]
  source_sorted <- is_source[ord]
  index_sorted <- local_index[ord]

  left <- seq_len(length(ord) - 1L)
  right <- left + 1L
  same_key <- bp_sorted[left] == bp_sorted[right] &
    a1_sorted[left] == a1_sorted[right] &
    a2_sorted[left] == a2_sorted[right]
  opposite_side <- source_sorted[left] != source_sorted[right]
  pair_left <- left[same_key & opposite_side]
  if (!length(pair_left)) {
    return(list(ref_index = integer(), src_index = integer()))
  }
  pair_right <- pair_left + 1L
  left_is_source <- source_sorted[pair_left]
  list(
    ref_index = ifelse(left_is_source, index_sorted[pair_right], index_sorted[pair_left]),
    src_index = ifelse(left_is_source, index_sorted[pair_left], index_sorted[pair_right])
  )
}
