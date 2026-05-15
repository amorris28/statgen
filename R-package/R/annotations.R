.annotations_cache_schema <- "annotations_cache/0.1"

annomat <- function(x, ...) UseMethod("annomat")
annonames <- function(x, ...) UseMethod("annonames")
num_annot <- function(x, ...) UseMethod("num_annot")
select_annotations <- function(x, names, ...) UseMethod("select_annotations")
union_annotations <- function(x, other, mode = "by_name", ...) UseMethod("union_annotations")

load_annotations <- function(bed_paths, reference) {
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  paths <- .coerce_bed_paths(bed_paths)
  for (path in paths) {
    if (!file.exists(path)) {
      stop(sprintf("BED file not found: %s", path), call. = FALSE)
    }
  }
  names <- .annotation_names_from_paths(paths)
  if (anyDuplicated(names)) {
    stop("duplicate annotation names derived from BED basenames", call. = FALSE)
  }

  columns <- lapply(paths, function(path) {
    Matrix::Matrix(.paint_bed_to_reference(.parse_bed(path), reference), ncol = 1L, sparse = TRUE)
  })
  create_annotations(reference, annotation_matrix = do.call(cbind, columns), annotation_names = names)
}

create_annotation <- function(reference, annovec, annoname) {
  name <- as.character(annoname)
  if (length(name) != 1L || is.na(name) || identical(name, "")) {
    stop("annoname must be a non-empty character scalar", call. = FALSE)
  }
  vec <- .coerce_annovec(annovec, num_snp(reference))
  create_annotations(reference, annotation_matrix = Matrix::Matrix(vec, ncol = 1L, sparse = TRUE), annotation_names = name)
}

create_annotations <- function(reference, annotation_matrix, annotation_names) {
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  names <- .coerce_annonames(annotation_names)
  mat <- .as_lgC_binary_matrix(annotation_matrix, "annotation_matrix")
  expected <- c(num_snp(reference), length(names))
  if (!identical(as.integer(dim(mat)), as.integer(expected))) {
    stop(sprintf(
      "annotation_matrix shape mismatch: expected (%d, %d), got (%d, %d)",
      expected[[1]], expected[[2]], dim(mat)[[1]], dim(mat)[[2]]
    ), call. = FALSE)
  }

  out <- vector("list", length(shards(reference)))
  offsets <- shard_offsets(reference)
  for (i in seq_along(shards(reference))) {
    ref_shard <- shards(reference)[[i]]
    start <- offsets$start0[[i]] + 1L
    stop <- offsets$stop0[[i]]
    out[[i]] <- .new_annotation_shard(
      ref_shard$label,
      ref_shard$checksum,
      mat[start:stop, , drop = FALSE]
    )
  }
  .new_annotation_panel(out, names)
}

save_annotations_cache <- function(panel, path) {
  if (!inherits(panel, "AnnotationPanel")) {
    stop("panel must be an AnnotationPanel", call. = FALSE)
  }
  path <- .validate_path_scalar(path, "path")
  offsets <- shard_offsets(panel)
  payload <- list(
    metadata = list(
      schema = .annotations_cache_schema,
      n_shards = length(panel$shards),
      shard_labels = vapply(panel$shards, function(s) s$label, character(1)),
      shard_checksums = vapply(panel$shards, function(s) s$reference_checksum, character(1)),
      shard_start0 = offsets$start0,
      shard_stop0 = offsets$stop0
    ),
    annomat = annomat(panel),
    annonames = annonames(panel)
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, path)
  invisible(NULL)
}

load_annotations_cache <- function(path, shards = NULL) {
  path <- .validate_path_scalar(path, "path")
  payload <- readRDS(path)
  cache <- .validate_annotations_cache_payload(payload)
  meta <- payload$metadata
  selected <- .validate_requested_shards(shards, meta$shard_labels, "load_annotations_cache")
  label_to_index <- stats::setNames(seq_along(meta$shard_labels), meta$shard_labels)
  names <- .coerce_annonames(payload$annonames)

  out <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    label <- selected[[i]]
    idx <- label_to_index[[label]]
    start <- meta$shard_start0[[idx]] + 1L
    stop <- meta$shard_stop0[[idx]]
    out[[i]] <- .new_annotation_shard(
      label,
      meta$shard_checksums[[idx]],
      cache$annomat[start:stop, , drop = FALSE]
    )
  }
  .new_annotation_panel(out, names)
}

num_snp.AnnotationShard <- function(x, ...) x$num_snp
num_snp.AnnotationPanel <- function(x, ...) x$num_snp
num_annot.AnnotationShard <- function(x, ...) x$num_annot
num_annot.AnnotationPanel <- function(x, ...) x$num_annot
shards.AnnotationPanel <- function(x, ...) x$shards
shard_offsets.AnnotationPanel <- function(x, ...) x$shard_offsets
annomat.AnnotationShard <- function(x, ...) x$annomat
annomat.AnnotationPanel <- function(x, ...) {
  mats <- lapply(x$shards, annomat)
  mat <- do.call(rbind, mats)
  .set_annotation_colnames(mat, x$annonames)
}
annonames.AnnotationPanel <- function(x, ...) x$annonames

select_shards.AnnotationPanel <- function(x, shards, ...) {
  available <- vapply(x$shards, function(s) s$label, character(1))
  selected <- .validate_requested_shards(shards, available, "AnnotationPanel.select_shards")
  by_label <- stats::setNames(x$shards, available)
  .new_annotation_panel(unname(by_label[selected]), x$annonames)
}

select_annotations.AnnotationPanel <- function(x, names, ...) {
  requested <- .coerce_annonames(names)
  idx <- match(requested, x$annonames)
  missing <- requested[is.na(idx)]
  if (length(missing)) {
    stop(sprintf("unknown annotation name(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  out <- lapply(x$shards, function(s) {
    .new_annotation_shard(s$label, s$reference_checksum, s$annomat[, idx, drop = FALSE])
  })
  .new_annotation_panel(out, requested)
}

union_annotations.AnnotationPanel <- function(x, other, mode = "by_name", ...) {
  if (!identical(mode, "by_name")) {
    stop("union_annotations supports only mode = 'by_name'", call. = FALSE)
  }
  if (!inherits(other, "AnnotationPanel")) {
    stop("other must be an AnnotationPanel", call. = FALSE)
  }
  overlap <- intersect(x$annonames, other$annonames)
  if (length(overlap)) {
    stop(sprintf("annotation name collision(s): %s", paste(overlap, collapse = ", ")), call. = FALSE)
  }
  .check_annotation_reference_alignment(x, other, "union_annotations")
  out <- vector("list", length(x$shards))
  for (i in seq_along(x$shards)) {
    a <- x$shards[[i]]
    b <- other$shards[[i]]
    out[[i]] <- .new_annotation_shard(
      a$label,
      a$reference_checksum,
      do.call(cbind, list(a$annomat, b$annomat))
    )
  }
  .new_annotation_panel(out, c(x$annonames, other$annonames))
}

.check_annotation_reference_alignment <- function(x, other, context) {
  if (length(x$shards) != length(other$shards)) {
    stop(sprintf(
      "%s requires compatible reference alignment: shard count mismatch",
      context
    ), call. = FALSE)
  }
  for (i in seq_along(x$shards)) {
    a <- x$shards[[i]]
    b <- other$shards[[i]]
    if (!identical(a$label, b$label)) {
      stop(sprintf(
        "%s requires compatible reference alignment: shard label mismatch: %s vs %s",
        context, sQuote(a$label), sQuote(b$label)
      ), call. = FALSE)
    }
    if (!identical(as.integer(a$num_snp), as.integer(b$num_snp))) {
      stop(sprintf(
        "%s requires compatible reference alignment: shard %s row count mismatch",
        context, sQuote(a$label)
      ), call. = FALSE)
    }
    if (!identical(a$reference_checksum, b$reference_checksum)) {
      stop(sprintf(
        "%s requires compatible reference alignment: shard %s reference_checksum mismatch",
        context, sQuote(a$label)
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

save_cache.AnnotationPanel <- function(x, path, ...) save_annotations_cache(x, path)

print.AnnotationPanel <- function(x, ...) {
  labels <- vapply(x$shards, function(s) s$label, character(1))
  cat(sprintf(
    "<AnnotationPanel: %d SNPs x %d annotations across %d shard(s): %s>\n",
    x$num_snp, x$num_annot, length(labels), paste(labels, collapse = ", ")
  ))
  invisible(x)
}

.new_annotation_shard <- function(label, reference_checksum, mat) {
  mat <- .as_lgC_binary_matrix(mat)
  structure(
    list(
      label = as.character(label),
      reference_checksum = as.character(reference_checksum),
      num_snp = as.integer(dim(mat)[[1]]),
      num_annot = as.integer(dim(mat)[[2]]),
      annomat = mat
    ),
    class = "AnnotationShard"
  )
}

.new_annotation_panel <- function(shard_list, names) {
  if (!length(shard_list)) {
    stop("AnnotationPanel requires at least one shard", call. = FALSE)
  }
  names <- .coerce_annonames(names)
  pos <- 0L
  offsets <- data.frame(
    shard_label = character(length(shard_list)),
    start0 = integer(length(shard_list)),
    stop0 = integer(length(shard_list)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(shard_list)) {
    if (!identical(as.integer(shard_list[[i]]$num_annot), as.integer(length(names)))) {
      stop("all annotation shards must share the same annotation columns", call. = FALSE)
    }
    offsets$shard_label[[i]] <- shard_list[[i]]$label
    offsets$start0[[i]] <- pos
    pos <- pos + shard_list[[i]]$num_snp
    offsets$stop0[[i]] <- pos
  }
  structure(
    list(
      shards = shard_list,
      annonames = names,
      num_snp = pos,
      num_annot = length(names),
      shard_offsets = offsets
    ),
    class = "AnnotationPanel"
  )
}

.set_annotation_colnames <- function(mat, names) {
  colnames(mat) <- names
  mat
}

.coerce_bed_paths <- function(bed_paths) {
  if (!is.character(bed_paths) || !length(bed_paths) || anyNA(bed_paths) || any(bed_paths == "")) {
    stop("bed_paths must be a non-empty character vector of BED files", call. = FALSE)
  }
  if (is.character(bed_paths) && length(bed_paths) == 1L) {
    return(bed_paths)
  }
  as.character(bed_paths)
}

.annotation_names_from_paths <- function(paths) {
  tools::file_path_sans_ext(basename(paths))
}

.coerce_annonames <- function(names) {
  if (!is.character(names) || !length(names) || anyNA(names) || any(names == "")) {
    stop("annonames must be a non-empty character vector of unique strings", call. = FALSE)
  }
  if (anyDuplicated(names)) {
    stop("annonames must be unique", call. = FALSE)
  }
  as.character(names)
}

.coerce_annovec <- function(annovec, n) {
  vec <- as.vector(annovec)
  if (length(vec) != n) {
    stop(sprintf("annovec length mismatch: expected %d, got %d", n, length(vec)), call. = FALSE)
  }
  if (any(is.na(vec)) || any(!(vec %in% c(FALSE, TRUE, 0, 1)))) {
    stop("annovec must be binary (0/1)", call. = FALSE)
  }
  as.logical(vec)
}

.as_lgC_binary_matrix <- function(mat, name = "annomat") {
  if (inherits(mat, "lgCMatrix")) {
    return(mat)
  }
  if (inherits(mat, "Matrix")) {
    dims <- dim(mat)
    entries <- Matrix::summary(mat)
    if (nrow(entries)) {
      nonzero <- entries$x != 0
      bad <- !(entries$x[nonzero] %in% c(TRUE, 1))
      if (any(bad)) {
        stop(sprintf("%s contains non-binary values", name), call. = FALSE)
      }
      return(Matrix::sparseMatrix(
        i = entries$i[nonzero],
        j = entries$j[nonzero],
        x = rep.int(TRUE, sum(nonzero)),
        dims = dims
      ))
    }
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = logical(), dims = dims))
  }
  arr <- as.matrix(mat)
  if (length(dim(arr)) != 2L) {
    stop(sprintf("%s must be a 2D matrix", name), call. = FALSE)
  }
  if (any(is.na(arr)) || any(!(arr %in% c(FALSE, TRUE, 0, 1)))) {
    stop(sprintf("%s contains non-binary values", name), call. = FALSE)
  }
  Matrix::Matrix(arr != 0, sparse = TRUE)
}

.parse_bed <- function(path) {
  # First pass: scan line-by-line to validate structure and count header lines to
  # skip. Avoids seek(), which R documents as unreliable on Windows text connections.
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  skip_count <- 0L
  first_data <- NULL
  found_data <- FALSE

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) {
      if (!found_data) {
        stop(sprintf("%s: BED file is empty", path), call. = FALSE)
      }
      break
    }
    if (!found_data) {
      if (identical(line, "") || startsWith(line, "#")) {
        skip_count <- skip_count + 1L
      } else {
        first_data <- line
        found_data <- TRUE
      }
    } else {
      if (identical(line, "") || startsWith(line, "#")) {
        stop(sprintf("%s: blank or comment line after BED data row", path), call. = FALSE)
      }
    }
  }
  close(con)
  on.exit(NULL)

  if (length(strsplit(first_data, "\t", fixed = TRUE)[[1L]]) < 3L) {
    stop(sprintf("%s: BED must have at least 3 tab-separated columns", path), call. = FALSE)
  }

  # Second pass: read.table opens the file fresh, no seek() needed.
  df <- tryCatch(
    utils::read.table(
      file = path,
      skip = skip_count,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      colClasses = "character",
      fill = FALSE
    ),
    error = function(e) {
      stop(sprintf("%s: malformed BED", path), call. = FALSE)
    }
  )
  if (ncol(df) < 3L) {
    stop(sprintf("%s: BED must have at least 3 tab-separated columns", path), call. = FALSE)
  }
  chr <- df[[1]]
  bad_chr <- is.na(chr) | chr == ""
  if (any(bad_chr)) {
    stop(sprintf("%s: row %d: chromosome label must be non-empty", path, which(bad_chr)[[1]]), call. = FALSE)
  }
  starts <- suppressWarnings(as.numeric(df[[2]]))
  ends <- suppressWarnings(as.numeric(df[[3]]))
  bad_start <- is.na(starts) | !is.finite(starts) | floor(starts) != starts | starts < 0
  if (any(bad_start)) {
    stop(sprintf("%s: row %d: BED start must be a non-negative integer", path, which(bad_start)[[1]]), call. = FALSE)
  }
  bad_end <- is.na(ends) | !is.finite(ends) | floor(ends) != ends | ends < 0
  if (any(bad_end)) {
    stop(sprintf("%s: row %d: BED end must be a non-negative integer", path, which(bad_end)[[1]]), call. = FALSE)
  }
  if (any(ends < starts)) {
    stop(sprintf("%s: row %d: BED interval end must be >= start", path, which(ends < starts)[[1]]), call. = FALSE)
  }
  split_idx <- split(seq_along(chr), chr)
  lapply(split_idx, function(idx) .merge_intervals(as.integer(starts[idx]), as.integer(ends[idx])))
}

.merge_intervals <- function(starts, ends) {
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

.paint_bed_to_reference <- function(intervals_by_chr, reference) {
  out <- rep.int(FALSE, num_snp(reference))
  offsets <- shard_offsets(reference)
  for (i in seq_along(shards(reference))) {
    ref_shard <- shards(reference)[[i]]
    intervals <- intervals_by_chr[[ref_shard$label]]
    if (is.null(intervals) || !nrow(intervals)) {
      next
    }
    start <- offsets$start0[[i]] + 1L
    stop <- offsets$stop0[[i]]
    out[start:stop] <- .paint_mask(bp(ref_shard), intervals)
  }
  out
}

.paint_mask <- function(bp, intervals) {
  pos0 <- as.integer(bp) - 1L
  idx <- findInterval(pos0, intervals[, "start"])
  valid <- idx > 0L
  out <- rep.int(FALSE, length(bp))
  out[valid] <- pos0[valid] < intervals[idx[valid], "end"]
  out
}

.validate_annotations_cache_payload <- function(payload) {
  cache <- .validate_cache_payload_metadata(payload, .annotations_cache_schema, "annotations")
  total <- cache$total
  if (is.null(payload$annomat) || is.null(payload$annonames)) {
    stop("Invalid annotations cache: missing annomat or annonames", call. = FALSE)
  }
  annomat <- .as_lgC_binary_matrix(payload$annomat)
  if (!identical(as.integer(dim(annomat)[[1]]), as.integer(total))) {
    stop("Invalid annotations cache: annomat row count mismatch", call. = FALSE)
  }
  if (!identical(as.integer(dim(annomat)[[2]]), as.integer(length(payload$annonames)))) {
    stop("Invalid annotations cache: annomat column count mismatch", call. = FALSE)
  }
  list(annomat = annomat)
}
