.parse_bim <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  df <- tryCatch(
    utils::read.table(
      path,
      header = FALSE,
      sep = "",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      colClasses = rep("character", 6L),
      fill = FALSE
    ),
    error = function(e) {
      stop(sprintf("%s: expected 6 whitespace-delimited columns", path), call. = FALSE)
    }
  )
  if (ncol(df) != 6L) {
    stop(sprintf("%s: expected 6 whitespace-delimited columns, got %d", path, ncol(df)), call. = FALSE)
  }
  names(df) <- c("chr", "snp", "cm", "bp", "a1", "a2")
  bad_chr <- is.na(df$chr) | df$chr == ""
  if (any(bad_chr)) {
    stop(sprintf("%s:%d: chr must be non-empty", path, which(bad_chr)[[1]]), call. = FALSE)
  }
  chr_style <- startsWith(tolower(df$chr), "chr")
  if (any(chr_style)) {
    stop(sprintf("%s:%d: chr-style labels (e.g., chr1/chrX) are not allowed", path, which(chr_style)[[1]]), call. = FALSE)
  }
  known_chr <- df$chr %in% c(.canonical_chr_order, .ignored_chr)
  if (!all(known_chr)) {
    idx <- which(!known_chr)[[1]]
    stop(sprintf("%s:%d: unsupported chr label %s; expected 1-22, X (Y/MT are ignored)", path, idx, sQuote(df$chr[[idx]])), call. = FALSE)
  }
  bad_allele <- is.na(df$a1) | is.na(df$a2) | df$a1 == "" | df$a2 == ""
  if (any(bad_allele)) {
    stop(sprintf("%s:%d: a1 and a2 must be non-empty", path, which(bad_allele)[[1]]), call. = FALSE)
  }
  bad_a1 <- !grepl("^[ACGT]+$", df$a1)
  if (any(bad_a1)) {
    idx <- which(bad_a1)[[1]]
    stop(sprintf("%s:%d: a1 must be uppercase DNA bases (A/C/G/T): %s", path, idx, sQuote(df$a1[[idx]])), call. = FALSE)
  }
  bad_a2 <- !grepl("^[ACGT]+$", df$a2)
  if (any(bad_a2)) {
    idx <- which(bad_a2)[[1]]
    stop(sprintf("%s:%d: a2 must be uppercase DNA bases (A/C/G/T): %s", path, idx, sQuote(df$a2[[idx]])), call. = FALSE)
  }
  cm <- suppressWarnings(as.numeric(df$cm))
  if (any(is.na(cm) | !is.finite(cm))) {
    idx <- which(is.na(cm) | !is.finite(cm))[[1]]
    stop(sprintf("%s:%d: cm is not a number: %s", path, idx, sQuote(df$cm[[idx]])), call. = FALSE)
  }
  bp_num <- suppressWarnings(as.numeric(df$bp))
  bad_bp <- is.na(bp_num) | !is.finite(bp_num) | floor(bp_num) != bp_num
  if (any(bad_bp)) {
    idx <- which(bad_bp)[[1]]
    stop(sprintf("%s:%d: bp is not an integer: %s", path, idx, sQuote(df$bp[[idx]])), call. = FALSE)
  }
  df$bp <- as.integer(bp_num)
  df$a1_hash64 <- .allele_hash64(df$a1)
  df$a2_hash64 <- .allele_hash64(df$a2)
  .validate_bim_order_and_duplicates(df[df$chr %in% .canonical_chr_order, , drop = FALSE], path)
  df
}

.validate_bim_order_and_duplicates <- function(df, path) {
  if (!nrow(df)) {
    return(invisible(NULL))
  }
  rank <- unname(.chr_rank[df$chr])
  bad_order <- c(FALSE, rank[-1L] < rank[-length(rank)] | (rank[-1L] == rank[-length(rank)] & df$bp[-1L] < df$bp[-nrow(df)]))
  if (any(bad_order)) {
    stop(sprintf("%s:%d: rows must be sorted by (chr_rank, bp) in canonical contig order", path, which(bad_order)[[1]]), call. = FALSE)
  }

  for (label in unique(df$chr)) {
    keep <- df$chr == label
    .check_unique_variant_keys(
      df$bp[keep],
      df$a1_hash64[keep],
      df$a2_hash64[keep],
      sprintf("%s: shard %s", path, label),
      "variant"
    )
  }
  invisible(NULL)
}
