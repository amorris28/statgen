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
