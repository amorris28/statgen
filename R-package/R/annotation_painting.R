.annotation_merge_intervals <- function(starts, ends) {
  if (!length(starts)) {
    return(matrix(integer(), ncol = 2L, dimnames = list(NULL, c("start", "end"))))
  }
  ord <- order(starts, ends)
  starts <- starts[ord]
  ends <- ends[ord]

  running_max_end <- cummax(ends)
  new_group <- c(TRUE, starts[-1L] > running_max_end[-length(running_max_end)])
  group_id <- cumsum(new_group)

  cbind(
    start = starts[new_group],
    end = as.integer(unname(tapply(ends, group_id, max)))
  )
}

.annotation_binary_intervals_by_chr <- function(chr_values, starts, ends) {
  keep <- !(chr_values %in% .ignored_chr)
  if (!any(keep)) {
    return(list())
  }
  split_idx <- split(which(keep), chr_values[keep])
  lapply(split_idx, function(idx) .annotation_merge_intervals(starts[idx], ends[idx]))
}

.annotation_validate_numeric_nonoverlap <- function(chr_values, starts, ends, path, row_base0) {
  keep <- !(chr_values %in% .ignored_chr)
  if (!any(keep)) {
    return(invisible(NULL))
  }
  split_idx <- split(which(keep), chr_values[keep])
  for (idx in split_idx) {
    if (length(idx) <= 1L) {
      next
    }
    ord <- order(starts[idx], ends[idx])
    sorted_idx <- idx[ord]
    sorted_starts <- starts[sorted_idx]
    sorted_ends <- ends[sorted_idx]
    prev_max <- cummax(sorted_ends)[-length(sorted_ends)]
    overlap <- sorted_starts[-1L] < prev_max
    if (any(overlap)) {
      row <- sorted_idx[which(overlap)[[1]] + 1L]
      stop(sprintf(
        "%s: row %d: numeric annotation intervals overlap on chromosome %s",
        path, row_base0 + row, chr_values[[row]]
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

.annotation_paint_mask <- function(bp, intervals) {
  pos0 <- as.integer(bp) - 1L
  idx <- findInterval(pos0, intervals[, "start"])
  valid <- idx > 0L
  out <- rep.int(0, length(bp))
  out[valid] <- as.integer(pos0[valid] < intervals[idx[valid], "end"])
  out
}

.annotation_paint_binary_column <- function(intervals_by_chr, reference) {
  out <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(num_snp(reference), 1L)
  )
  offsets <- shard_offsets(reference)
  for (i in seq_along(shards(reference))) {
    ref_shard <- shards(reference)[[i]]
    intervals <- intervals_by_chr[[ref_shard$label]]
    if (is.null(intervals) || !nrow(intervals)) {
      next
    }
    start <- offsets$start0[[i]] + 1L
    stop <- offsets$stop0[[i]]
    mask <- .annotation_paint_mask(bp(ref_shard), intervals)
    nz <- which(mask != 0)
    if (length(nz)) {
      out[start + nz - 1L, 1L] <- 1
    }
  }
  methods::as(out, "dgCMatrix")
}

.annotation_paint_numeric_sparse <- function(bp, intervals, values) {
  n <- length(bp)
  k <- ncol(values)
  if (!length(intervals)) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, k)))
  }
  pos0 <- as.integer(bp) - 1L
  idx <- findInterval(pos0, intervals[, "start"])
  valid <- idx > 0L
  if (!any(valid)) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, k)))
  }
  valid_rows <- which(valid)
  inside <- pos0[valid] < intervals[idx[valid], "end"]
  if (!any(inside)) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, k)))
  }
  hit_rows <- valid_rows[inside]
  hit_values <- values[idx[valid][inside], , drop = FALSE]
  nz <- which(hit_values != 0, arr.ind = TRUE)
  if (!nrow(nz)) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, k)))
  }
  Matrix::sparseMatrix(
    i = hit_rows[nz[, "row"]],
    j = nz[, "col"],
    x = hit_values[nz],
    dims = c(n, k)
  )
}

.annotation_paint_numeric <- function(chr_values, starts, ends, values, reference) {
  k <- ncol(values)
  blocks <- lapply(shards(reference), function(ref_shard) {
    Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(num_snp(ref_shard), k))
  })
  keep <- !(chr_values %in% .ignored_chr)
  if (any(keep)) {
    split_idx <- split(which(keep), chr_values[keep])
    labels <- vapply(shards(reference), function(s) s$label, character(1))
    for (label in names(split_idx)) {
      shard_idx <- match(label, labels)
      if (is.na(shard_idx)) {
        next
      }
      idx <- split_idx[[label]]
      ord <- order(starts[idx], ends[idx])
      sorted_idx <- idx[ord]
      intervals <- cbind(start = starts[sorted_idx], end = ends[sorted_idx])
      blocks[[shard_idx]] <- .annotation_paint_numeric_sparse(
        bp(shards(reference)[[shard_idx]]),
        intervals,
        values[sorted_idx, , drop = FALSE]
      )
    }
  }
  methods::as(do.call(rbind, blocks), "dgCMatrix")
}
