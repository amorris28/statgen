.annotations_cache_schema <- "annotations_cache/0.2"
.annotations_old_cache_schema <- "annotations_cache/0.1"

annomat <- function(x, ...) UseMethod("annomat")
annonames <- function(x, ...) UseMethod("annonames")
num_annot <- function(x, ...) UseMethod("num_annot")
is_binary <- function(x, ...) UseMethod("is_binary")
annotation_metadata <- function(x, ...) UseMethod("annotation_metadata")
select_annotations <- function(x, names, ...) UseMethod("select_annotations")
union_annotations <- function(x, other, mode = "by_name", ...) UseMethod("union_annotations")

load_annotations <- function(bed_paths, reference,
                             annotation_metadata = NULL,
                             annotation_metadata_paths = NULL) {
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  if (!is.null(annotation_metadata) && !is.null(annotation_metadata_paths)) {
    stop("load_annotations accepts at most one of annotation_metadata and annotation_metadata_paths", call. = FALSE)
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

  parsed <- lapply(paths, .parse_binary_bed)
  columns <- lapply(parsed, function(x) {
    .annotation_paint_binary_column(x$intervals, reference)
  })
  source_counts <- vapply(parsed, function(x) x$num_source_intervals, integer(1))
  metadata <- .load_annotations_metadata(paths, annotation_metadata, annotation_metadata_paths, source_counts)
  create_annotations(
    reference,
    annotation_matrix = do.call(cbind, columns),
    annotation_names = names,
    is_binary = rep.int(TRUE, length(names)),
    annotation_metadata = metadata
  )
}

load_annotation <- function(path, reference, has_header = FALSE, value_columns = NULL,
                            group_column = NULL,
                            annotation_names = NULL, annotation_metadata = NULL,
                            annotation_metadata_path = NULL) {
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  if (!is.null(annotation_metadata) && !is.null(annotation_metadata_path)) {
    stop("load_annotation accepts at most one of annotation_metadata and annotation_metadata_path", call. = FALSE)
  }
  if (!is.null(value_columns) && !is.null(group_column)) {
    stop("group_column and value_columns are mutually exclusive", call. = FALSE)
  }
  if (!is.null(group_column) && !is.null(annotation_names)) {
    stop("annotation_names is invalid with group_column", call. = FALSE)
  }
  if (!is.null(group_column) && !is.null(annotation_metadata)) {
    stop("annotation_metadata is invalid with group_column", call. = FALSE)
  }
  if (!is.null(group_column) && !is.null(annotation_metadata_path)) {
    stop("annotation_metadata_path is invalid with group_column", call. = FALSE)
  }
  path <- .validate_path_scalar(path, "path")
  if (!file.exists(path)) {
    stop(sprintf("annotation file not found: %s", path), call. = FALSE)
  }

  table <- .read_annotation_table(path, isTRUE(has_header), value_columns, group_column, infer_single_value = TRUE)
  df <- table$df
  n_cols <- ncol(df)
  value_columns1 <- table$value_columns
  value_column_metadata <- table$value_column_metadata
  group_column1 <- table$group_column
  group_column_metadata <- if (is.null(table$group_column_metadata)) NULL else table$group_column_metadata[[1L]]
  interval <- .validate_annotation_intervals(df, path, table$row_base0)
  chr_values <- interval$chr
  starts <- interval$starts
  ends <- interval$ends

  if (!is.null(group_column1)) {
    group_values <- as.character(df[[group_column1]])
    empty <- is.na(group_values) | group_values == ""
    if (any(empty)) {
      stop(sprintf("%s: row %d: group_column values must be non-empty", path, table$row_base0 + which(empty)[[1]]), call. = FALSE)
    }
    names <- unique(group_values)
    columns <- vector("list", length(names))
    metadata <- character(length(names))
    # PERF: loop retained because grouped binary interval union is defined independently per group.
    for (i in seq_along(names)) {
      rows <- group_values == names[[i]]
      columns[[i]] <- .annotation_paint_binary_column(
        .annotation_binary_intervals_by_chr(chr_values[rows], starts[rows], ends[rows]),
        reference
      )
      metadata[[i]] <- .annotation_generated_metadata(
        path,
        isTRUE(has_header),
        sum(rows),
        group_column = group_column_metadata,
        group_value = names[[i]]
      )
    }
    mat <- do.call(cbind, columns)
    binary <- rep.int(TRUE, length(names))
  } else if (is.null(value_columns1) && n_cols == 3L) {
    names <- if (is.null(annotation_names)) {
      .default_annotation_names(path, n_cols, table$header_fields, value_columns1)
    } else {
      .coerce_annonames(annotation_names, "annotation_names")
    }
    if (length(names) != 1L) {
      stop("annotation_names length mismatch: expected 1", call. = FALSE)
    }
    mat <- .annotation_paint_binary_column(
      .annotation_binary_intervals_by_chr(chr_values, starts, ends),
      reference
    )
    binary <- TRUE
    metadata <- .load_annotation_metadata(
      path, annotation_metadata, annotation_metadata_path, 1L,
      NULL, NULL, n_cols, binary = TRUE,
      has_header = isTRUE(has_header), num_source_intervals = nrow(df)
    )
  } else {
    if (is.null(value_columns1)) {
      if (n_cols == 4L) {
        value_columns1 <- 4L
        value_column_metadata <- if (is.null(table$header_fields)) 4L else table$header_fields[[4L]]
      } else {
        stop(sprintf("%s: input with five or more columns requires explicit value_columns or group_column", path), call. = FALSE)
      }
    }
    names <- if (is.null(annotation_names)) {
      .default_annotation_names(path, n_cols, table$header_fields, value_columns1)
    } else {
      .coerce_annonames(annotation_names, "annotation_names")
    }
    if (length(names) != length(value_columns1)) {
      stop(sprintf(
        "annotation_names length mismatch: expected %d, got %d",
        length(value_columns1), length(names)
      ), call. = FALSE)
    }

    .annotation_validate_numeric_nonoverlap(chr_values, starts, ends, path, table$row_base0)
    values <- .annotation_numeric_values(df, value_columns1, path, table$row_base0)
    mat <- .annotation_paint_numeric(chr_values, starts, ends, values, reference)
    binary <- rep.int(FALSE, length(value_columns1))
    metadata <- .load_annotation_metadata(
      path, annotation_metadata, annotation_metadata_path, length(value_columns1),
      value_columns1, value_column_metadata, n_cols, binary = FALSE,
      has_header = isTRUE(has_header), num_source_intervals = nrow(df)
    )
  }

  create_annotations(
    reference,
    annotation_matrix = mat,
    annotation_names = names,
    is_binary = binary,
    annotation_metadata = metadata
  )
}

create_annotation <- function(reference, annovec, annotation_name,
                              is_binary = NULL, annotation_metadata = NULL) {
  name <- as.character(annotation_name)
  if (length(name) != 1L || is.na(name) || identical(name, "")) {
    stop("annotation_name must be a non-empty character scalar", call. = FALSE)
  }
  vec <- .coerce_annovec(annovec, num_snp(reference))
  metadata <- if (is.null(annotation_metadata)) "" else as.character(annotation_metadata)
  if (length(metadata) != 1L || is.na(metadata)) {
    stop("annotation_metadata must be a character scalar", call. = FALSE)
  }
  create_annotations(
    reference,
    annotation_matrix = Matrix::Matrix(vec, ncol = 1L, sparse = TRUE),
    annotation_names = name,
    is_binary = is_binary,
    annotation_metadata = metadata
  )
}

create_annotations <- function(reference, annotation_matrix, annotation_names,
                               is_binary = NULL, annotation_metadata = NULL) {
  if (!inherits(reference, "ReferencePanel")) {
    stop("reference must be a ReferencePanel", call. = FALSE)
  }
  names <- .coerce_annonames(annotation_names)
  mat <- .as_dgC_numeric_matrix(annotation_matrix, "annotation_matrix")
  expected <- c(num_snp(reference), length(names))
  if (!identical(as.integer(dim(mat)), as.integer(expected))) {
    stop(sprintf(
      "annotation_matrix shape mismatch: expected (%d, %d), got (%d, %d)",
      expected[[1]], expected[[2]], dim(mat)[[1]], dim(mat)[[2]]
    ), call. = FALSE)
  }

  binary <- .coerce_or_infer_is_binary(is_binary, mat)
  .validate_declared_binary(mat, binary, "annotation_matrix")
  metadata <- if (is.null(annotation_metadata)) {
    rep.int("", length(names))
  } else {
    .coerce_annotation_metadata(annotation_metadata, length(names), "annotation_metadata")
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
  .new_annotation_panel(out, names, binary, metadata)
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
    annonames = annonames(panel),
    is_binary = is_binary(panel),
    annotation_metadata = annotation_metadata(panel)
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
      cache$annomat[start:stop, , drop = FALSE],
      validate_values = FALSE
    )
  }
  .new_annotation_panel(out, names, cache$is_binary, cache$annotation_metadata)
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
is_binary.AnnotationPanel <- function(x, ...) x$is_binary
annotation_metadata.AnnotationPanel <- function(x, ...) x$annotation_metadata

select_shards.AnnotationPanel <- function(x, shards, ...) {
  available <- vapply(x$shards, function(s) s$label, character(1))
  selected <- .validate_requested_shards(shards, available, "AnnotationPanel.select_shards")
  by_label <- stats::setNames(x$shards, available)
  .new_annotation_panel(unname(by_label[selected]), x$annonames, x$is_binary, x$annotation_metadata)
}

select_annotations.AnnotationPanel <- function(x, names, ...) {
  requested <- .coerce_annonames(names, "names")
  idx <- match(requested, x$annonames)
  missing <- requested[is.na(idx)]
  if (length(missing)) {
    stop(sprintf("unknown annotation name(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  out <- lapply(x$shards, function(s) {
    .new_annotation_shard(s$label, s$reference_checksum, s$annomat[, idx, drop = FALSE])
  })
  .new_annotation_panel(out, requested, x$is_binary[idx], x$annotation_metadata[idx])
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
  .new_annotation_panel(
    out,
    c(x$annonames, other$annonames),
    c(x$is_binary, other$is_binary),
    c(x$annotation_metadata, other$annotation_metadata)
  )
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
    if (is.null(b$reference_checksum) || !identical(a$reference_checksum, b$reference_checksum)) {
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

.new_annotation_shard <- function(label, reference_checksum, mat, validate_values = TRUE) {
  mat <- .as_dgC_numeric_matrix(mat, validate_values = validate_values)
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

.new_annotation_panel <- function(shard_list, names, is_binary, annotation_metadata) {
  if (!length(shard_list)) {
    stop("AnnotationPanel requires at least one shard", call. = FALSE)
  }
  names <- .coerce_annonames(names)
  binary <- .coerce_is_binary_vector(is_binary, length(names))
  metadata <- .coerce_annotation_metadata(annotation_metadata, length(names), "annotation_metadata")
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
      # Store fields directly here; callers should use accessors outside constructors.
      is_binary = binary,
      annotation_metadata = metadata,
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
  as.character(bed_paths)
}

.annotation_names_from_paths <- function(paths) {
  tools::file_path_sans_ext(basename(paths))
}

.coerce_annonames <- function(names, name = "annonames") {
  if (!is.character(names) || !length(names) || anyNA(names) || any(names == "")) {
    stop(sprintf("%s must be a non-empty character vector of unique strings", name), call. = FALSE)
  }
  if (anyDuplicated(names)) {
    stop(sprintf("%s must be unique", name), call. = FALSE)
  }
  as.character(names)
}

.coerce_annotation_metadata <- function(value, expected_len, name) {
  if (!is.character(value)) {
    stop(sprintf("%s must be a character vector", name), call. = FALSE)
  }
  if (length(value) != expected_len) {
    stop(sprintf("%s length mismatch: expected %d, got %d", name, expected_len, length(value)), call. = FALSE)
  }
  if (anyNA(value)) {
    stop(sprintf("%s must not contain missing values", name), call. = FALSE)
  }
  as.character(value)
}

.coerce_annovec <- function(annovec, n) {
  vec <- as.numeric(as.vector(annovec))
  if (length(vec) != n) {
    stop(sprintf("annovec length mismatch: expected %d, got %d", n, length(vec)), call. = FALSE)
  }
  if (any(!is.finite(vec))) {
    stop(sprintf("annovec[%d] must be finite numeric", which(!is.finite(vec))[[1]]), call. = FALSE)
  }
  vec
}

.as_dgC_numeric_matrix <- function(mat, name = "annomat", validate_values = TRUE, drop_zeros = TRUE) {
  if (inherits(mat, "dgCMatrix")) {
    out <- mat
  } else if (inherits(mat, "Matrix")) {
    out <- Matrix::Matrix(mat, sparse = TRUE)
    if (drop_zeros) {
      out <- Matrix::drop0(out)
    }
    out <- methods::as(out, "dgCMatrix")
  } else {
    arr <- as.matrix(mat)
    if (length(dim(arr)) != 2L) {
      stop(sprintf("%s must be a 2D matrix", name), call. = FALSE)
    }
    if (validate_values && any(!is.finite(arr))) {
      idx <- which(!is.finite(arr), arr.ind = TRUE)[1, , drop = TRUE]
      stop(sprintf("%s contains non-finite value at row %d column %d", name, idx[[1]], idx[[2]]), call. = FALSE)
    }
    out <- Matrix::Matrix(arr, sparse = TRUE)
    if (drop_zeros) {
      out <- Matrix::drop0(out)
    }
    out <- methods::as(out, "dgCMatrix")
  }
  if (length(dim(out)) != 2L) {
    stop(sprintf("%s must be a 2D matrix", name), call. = FALSE)
  }
  if (validate_values) {
    entries <- Matrix::summary(out)
    if (nrow(entries) && any(!is.finite(entries$x))) {
      bad <- which(!is.finite(entries$x))[[1]]
      stop(sprintf("%s contains non-finite value: %s", name, entries$x[[bad]]), call. = FALSE)
    }
  }
  if (drop_zeros) {
    out <- Matrix::drop0(out)
  }
  methods::as(out, "dgCMatrix")
}

.coerce_or_infer_is_binary <- function(value, mat) {
  k <- dim(mat)[[2]]
  if (is.null(value)) {
    out <- rep.int(TRUE, k)
    entries <- Matrix::summary(mat)
    if (nrow(entries)) {
      out[unique(entries$j[entries$x != 1])] <- FALSE
    }
    return(out)
  }
  .coerce_is_binary_vector(value, k)
}

.coerce_is_binary_vector <- function(value, expected_len) {
  if (is.null(value)) {
    stop("is_binary must be supplied", call. = FALSE)
  }
  if (is.character(value)) {
    stop("is_binary must be a logical vector", call. = FALSE)
  }
  if (length(value) != expected_len) {
    stop(sprintf("is_binary length mismatch: expected %d, got %d", expected_len, length(value)), call. = FALSE)
  }
  if (anyNA(value)) {
    stop("is_binary must not contain missing values", call. = FALSE)
  }
  bad <- !(value %in% c(FALSE, TRUE, 0, 1))
  if (any(bad)) {
    stop(sprintf("is_binary[%d] must be logical", which(bad)[[1]]), call. = FALSE)
  }
  as.logical(value)
}

.validate_declared_binary <- function(mat, binary, name) {
  if (!any(binary)) {
    return(invisible(NULL))
  }
  entries <- Matrix::summary(mat[, binary, drop = FALSE])
  if (nrow(entries) && any(entries$x != 1)) {
    bad <- which(entries$x != 1)[[1]]
    stop(sprintf("%s declared binary must be binary (0/1); found non-binary value %s", name, entries$x[[bad]]), call. = FALSE)
  }
  invisible(NULL)
}

.parse_binary_bed <- function(path) {
  table <- .read_annotation_table(path, has_header = FALSE, value_columns = NULL, group_column = NULL, infer_single_value = FALSE)
  interval <- .validate_annotation_intervals(table$df, path, table$row_base0)
  list(
    intervals = .annotation_binary_intervals_by_chr(interval$chr, interval$starts, interval$ends),
    num_source_intervals = nrow(table$df)
  )
}

.read_annotation_table <- function(path, has_header, value_columns, group_column, infer_single_value) {
  probe <- .probe_annotation_table(path)
  n_cols <- probe$n_cols
  if (n_cols < 3L) {
    stop(sprintf("%s: BED must have at least 3 tab-separated columns", path), call. = FALSE)
  }
  header_fields <- NULL
  if (has_header) {
    header_fields <- strsplit(probe$first_line, "\t", fixed = TRUE)[[1L]]
    if (any(header_fields == "")) {
      stop(sprintf("%s: header names must be non-empty", path), call. = FALSE)
    }
    if (anyDuplicated(header_fields)) {
      stop(sprintf("%s: header names must be unique", path), call. = FALSE)
    }
  }
  value_normalized <- .normalize_column_selectors(value_columns, header_fields, n_cols, path, "value_columns")
  group_normalized <- .normalize_column_selectors(group_column, header_fields, n_cols, path, "group_column")
  value_columns1 <- value_normalized$columns
  value_column_metadata <- value_normalized$metadata
  group_column1 <- group_normalized$columns
  group_column_metadata <- group_normalized$metadata
  if (!is.null(value_columns1) && !is.null(group_column1)) {
    stop("group_column and value_columns are mutually exclusive", call. = FALSE)
  }
  if (!is.null(group_column1) && length(group_column1) != 1L) {
    stop("group_column must identify exactly one source column", call. = FALSE)
  }
  if (isTRUE(infer_single_value) && is.null(value_columns1) && is.null(group_column1) && n_cols == 4L) {
    value_columns1 <- 4L
    value_column_metadata <- if (is.null(header_fields)) 4L else header_fields[[4L]]
  }

  char_cols <- seq_len(n_cols)
  numeric_cols <- integer()
  if (!is.null(value_columns1)) {
    numeric_cols <- value_columns1
    char_cols <- setdiff(char_cols, numeric_cols)
  }
  col_names <- if (has_header) header_fields else paste0("V", seq_len(n_cols))
  fread_args <- list(
    input = path,
    sep = "\t",
    header = has_header,
    skip = probe$skip,
    fill = TRUE,
    blank.lines.skip = TRUE,
    colClasses = list(
      character = col_names[char_cols],
      numeric = col_names[numeric_cols]
    ),
    na.strings = c("NA", "NaN", "nan", "Inf", "inf", "-Inf", "-inf"),
    data.table = FALSE,
    showProgress = FALSE
  )
  df <- tryCatch(
    do.call(data.table::fread, fread_args),
    error = function(e) {
      stop(sprintf("%s: malformed annotation input", path), call. = FALSE)
    },
    warning = function(w) {
      stop(sprintf("%s: malformed annotation input", path), call. = FALSE)
    }
  )
  if (nrow(df) == 0L) {
    stop(sprintf("%s: BED file is empty", path), call. = FALSE)
  }
  if (ncol(df) != n_cols) {
    stop(sprintf("%s: malformed annotation input", path), call. = FALSE)
  }
  names(df) <- col_names
  .reject_late_annotation_comments(df, path)
  list(
    df = df,
    header_fields = header_fields,
    row_base0 = if (has_header) 1L else 0L,
    value_columns = value_columns1,
    value_column_metadata = value_column_metadata,
    group_column = group_column1,
    group_column_metadata = group_column_metadata
  )
}

.probe_annotation_table <- function(path) {
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) {
      stop(sprintf("%s: BED file is empty", path), call. = FALSE)
    }
    if (identical(line, "") || startsWith(line, "#")) {
      skip <- skip + 1L
      next
    }
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    return(list(skip = skip, first_line = line, n_cols = length(fields)))
  }
}

.reject_late_annotation_comments <- function(df, path) {
  first <- as.character(df[[1L]])
  late <- is.na(first) | first == "" | startsWith(first, "#")
  if (any(late)) {
    stop(sprintf("%s: blank or comment line after BED data row", path), call. = FALSE)
  }
}

.normalize_column_selectors <- function(selectors, header_fields, n_cols, path, name) {
  if (is.null(selectors)) {
    return(list(columns = NULL, metadata = NULL))
  }
  selector_list <- as.list(selectors)
  if (!length(selector_list)) {
    stop(sprintf("%s must be a non-empty vector", name), call. = FALSE)
  }
  out <- integer(length(selector_list))
  metadata <- vector("list", length(selector_list))
  for (i in seq_along(selector_list)) {
    selector <- selector_list[[i]]
    if (is.character(selector)) {
      if (is.null(header_fields)) {
        stop(sprintf("named %s are invalid when has_header = FALSE", name), call. = FALSE)
      }
      idx <- match(selector, header_fields)
      if (is.na(idx)) {
        stop(sprintf("%s: unknown %s name %s", path, name, sQuote(selector)), call. = FALSE)
      }
      col <- idx
      metadata[[i]] <- selector
    } else if (is.numeric(selector) && length(selector) == 1L && is.finite(selector) && selector == floor(selector)) {
      col <- as.integer(selector)
      metadata[[i]] <- col
    } else {
      stop(sprintf("%s must contain column names or integer indices", name), call. = FALSE)
    }
    if (col < 4L || col > n_cols) {
      stop(sprintf("%s: %s is out of range; selected columns must be physical columns 4 or later", path, name), call. = FALSE)
    }
    out[[i]] <- col
  }
  if (anyDuplicated(out)) {
    stop(sprintf("%s must not contain duplicates", name), call. = FALSE)
  }
  list(columns = out, metadata = metadata)
}

.validate_annotation_intervals <- function(df, path, row_base0) {
  chr_values <- as.character(df[[1L]])
  bad_chr <- is.na(chr_values) | chr_values == ""
  if (any(bad_chr)) {
    stop(sprintf("%s: row %d: chromosome label must be non-empty", path, row_base0 + which(bad_chr)[[1]]), call. = FALSE)
  }
  .validate_annotation_chr_labels(chr_values, path, row_base0)

  starts <- suppressWarnings(as.numeric(df[[2L]]))
  bad_start <- is.na(starts) | !is.finite(starts) | floor(starts) != starts | starts < 0
  if (any(bad_start)) {
    stop(sprintf("%s: row %d: BED start must be a non-negative integer", path, row_base0 + which(bad_start)[[1]]), call. = FALSE)
  }
  ends <- suppressWarnings(as.numeric(df[[3L]]))
  bad_end <- is.na(ends) | !is.finite(ends) | floor(ends) != ends | ends < 0
  if (any(bad_end)) {
    stop(sprintf("%s: row %d: BED end must be a non-negative integer", path, row_base0 + which(bad_end)[[1]]), call. = FALSE)
  }
  bad_len <- ends < starts
  if (any(bad_len)) {
    stop(sprintf("%s: row %d: BED interval end must be >= start", path, row_base0 + which(bad_len)[[1]]), call. = FALSE)
  }
  list(chr = chr_values, starts = as.integer(starts), ends = as.integer(ends))
}

.validate_annotation_chr_labels <- function(chr_values, path, row_base0) {
  chr_style <- startsWith(tolower(chr_values), "chr")
  if (any(chr_style)) {
    stop(sprintf("%s: row %d: chr-style labels (e.g., chr1/chrX) are not allowed", path, row_base0 + which(chr_style)[[1]]), call. = FALSE)
  }
  known <- chr_values %in% c(.canonical_chr_order, .ignored_chr)
  if (!all(known)) {
    idx <- which(!known)[[1]]
    stop(sprintf("%s: row %d: unsupported chr label %s; expected 1-22, X (Y/MT are ignored)", path, row_base0 + idx, sQuote(chr_values[[idx]])), call. = FALSE)
  }
}

.annotation_numeric_values <- function(df, value_columns1, path, row_base0) {
  values <- as.matrix(df[, value_columns1, drop = FALSE])
  storage.mode(values) <- "double"
  bad <- which(!is.finite(values), arr.ind = TRUE)
  if (nrow(bad)) {
    row <- bad[[1, "row"]]
    col <- value_columns1[[bad[[1, "col"]]]]
    stop(sprintf("%s: row %d: annotation value column %d must be finite numeric", path, row_base0 + row, col), call. = FALSE)
  }
  values
}

.default_annotation_names <- function(path, n_cols, header_fields, value_columns1) {
  if (is.null(value_columns1)) {
    if (n_cols == 3L) {
      return(.annotation_names_from_paths(path))
    }
    if (n_cols == 4L) {
      if (is.null(header_fields)) {
        return(.annotation_names_from_paths(path))
      }
      return(header_fields[[4L]])
    }
    stop(sprintf("%s: input with five or more columns requires explicit value_columns or group_column", path), call. = FALSE)
  }
  if (is.null(header_fields) && n_cols >= 5L) {
    stop(sprintf("%s: headerless input with five or more columns requires annotation_names", path), call. = FALSE)
  }
  if (is.null(header_fields)) {
    return(.annotation_names_from_paths(path))
  }
  header_fields[value_columns1]
}

.load_annotations_metadata <- function(paths, metadata, metadata_paths, source_counts) {
  if (!is.null(metadata)) {
    return(.coerce_annotation_metadata(metadata, length(paths), "annotation_metadata"))
  }
  if (!is.null(metadata_paths)) {
    if (!is.character(metadata_paths) || length(metadata_paths) != length(paths)) {
      stop(sprintf("annotation_metadata_paths length mismatch: expected %d, got %d", length(paths), length(metadata_paths)), call. = FALSE)
    }
    out <- character(length(paths))
    for (i in seq_along(paths)) {
      meta_path <- metadata_paths[[i]]
      if (is.na(meta_path) || identical(meta_path, "")) {
        out[[i]] <- .annotation_generated_metadata(paths[[i]], FALSE, source_counts[[i]])
      } else {
        out[[i]] <- .read_sidecar_exact(meta_path)
      }
    }
    return(out)
  }
  vapply(seq_along(paths), function(i) {
    .annotation_generated_metadata(paths[[i]], FALSE, source_counts[[i]])
  }, character(1))
}

.load_annotation_metadata <- function(path, metadata, metadata_path, expected_len,
                                      source_columns, source_column_metadata,
                                      n_cols, binary, has_header,
                                      num_source_intervals) {
  if (!is.null(metadata)) {
    return(.coerce_annotation_metadata(metadata, expected_len, "annotation_metadata"))
  }
  if (!is.null(metadata_path)) {
    metadata_path <- .validate_path_scalar(metadata_path, "annotation_metadata_path")
    if (binary) {
      return(.read_sidecar_exact(metadata_path))
    }
    lines <- .read_column_metadata_sidecar(metadata_path, n_cols)
    return(lines[source_columns])
  }
  if (binary) {
    return(.annotation_generated_metadata(path, has_header, num_source_intervals))
  }
  vapply(seq_len(expected_len), function(i) {
    .annotation_generated_metadata(
      path,
      has_header,
      num_source_intervals,
      value_column = source_column_metadata[[i]]
    )
  }, character(1))
}

.annotation_generated_metadata <- function(path, has_header, num_source_intervals,
                                           value_column = NULL,
                                           group_column = NULL,
                                           group_value = NULL) {
  payload <- list(
    source_file = as.character(path),
    source_file_has_header = if (isTRUE(has_header)) 1L else 0L,
    num_source_intervals = as.integer(num_source_intervals)
  )
  if (!is.null(value_column)) {
    payload$value_column <- value_column
  }
  if (!is.null(group_column)) {
    payload$group_column <- group_column
  }
  if (!is.null(group_value)) {
    payload$group_value <- group_value
  }
  jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null")
}

.read_sidecar_exact <- function(path) {
  path <- .validate_path_scalar(path, "annotation_metadata_path")
  size <- file.info(path)$size
  if (is.na(size)) {
    stop(sprintf("annotation metadata sidecar not found: %s", path), call. = FALSE)
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  readChar(con, nchars = size, useBytes = TRUE)
}

.read_column_metadata_sidecar <- function(path, n_cols) {
  text <- .read_sidecar_exact(path)
  lines <- strsplit(text, "\r\n|\n|\r", perl = TRUE)[[1L]]
  if (length(lines) == 1L && identical(lines, "")) {
    lines <- character()
  }
  if (length(lines) && identical(lines[[length(lines)]], "") && grepl("(\r\n|\n|\r)$", text, perl = TRUE)) {
    lines <- lines[-length(lines)]
  }
  if (length(lines) != n_cols) {
    stop(sprintf("%s: annotation metadata sidecar line count mismatch: expected %d, got %d", path, n_cols, length(lines)), call. = FALSE)
  }
  if (any(lines == "")) {
    stop(sprintf("%s: annotation metadata sidecar line %d is empty", path, which(lines == "")[[1]]), call. = FALSE)
  }
  lines
}

.validate_annotations_cache_payload <- function(payload) {
  if (!is.list(payload) || is.null(payload$metadata)) {
    stop("Invalid annotations cache: expected an RDS list with metadata", call. = FALSE)
  }
  schema <- as.character(payload$metadata$schema)
  if (identical(schema, .annotations_old_cache_schema)) {
    cache <- .validate_cache_payload_metadata(payload, .annotations_old_cache_schema, "annotations")
    names <- .coerce_annonames(payload$annonames)
    annomat <- .as_dgC_numeric_matrix(payload$annomat, validate_values = FALSE, drop_zeros = FALSE)
    if (!identical(as.integer(dim(annomat)[[1]]), as.integer(cache$total))) {
      stop("Invalid annotations cache: annomat row count mismatch", call. = FALSE)
    }
    if (!identical(as.integer(dim(annomat)[[2]]), as.integer(length(names)))) {
      stop("Invalid annotations cache: annomat column count mismatch", call. = FALSE)
    }
    return(list(
      annomat = annomat,
      is_binary = rep.int(TRUE, length(names)),
      annotation_metadata = rep.int("", length(names))
    ))
  }
  cache <- .validate_cache_payload_metadata(payload, .annotations_cache_schema, "annotations")
  total <- cache$total
  if (is.null(payload$annomat) || is.null(payload$annonames)) {
    stop("Invalid annotations cache: missing annomat or annonames", call. = FALSE)
  }
  names <- .coerce_annonames(payload$annonames)
  annomat <- .as_dgC_numeric_matrix(payload$annomat, validate_values = FALSE, drop_zeros = FALSE)
  if (!identical(as.integer(dim(annomat)[[1]]), as.integer(total))) {
    stop("Invalid annotations cache: annomat row count mismatch", call. = FALSE)
  }
  if (!identical(as.integer(dim(annomat)[[2]]), as.integer(length(names)))) {
    stop("Invalid annotations cache: annomat column count mismatch", call. = FALSE)
  }
  binary <- .coerce_is_binary_vector(payload$is_binary, length(names))
  metadata <- .coerce_annotation_metadata(payload$annotation_metadata, length(names), "annotation_metadata")
  list(annomat = annomat, is_binary = binary, annotation_metadata = metadata)
}
