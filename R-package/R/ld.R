reference <- function(x, ...) UseMethod("reference")
default_chrX_sex <- function(x, ...) UseMethod("default_chrX_sex")
a1freq <- function(x, chrX_sex = NULL, ...) UseMethod("a1freq")
multiply_r2 <- function(ld_panel, M, chrX_sex = NULL, ...) UseMethod("multiply_r2")

load_ld <- function(path, shards = NULL, default_chrX_sex = NULL, retain_ld_r = TRUE) {
  root <- .validate_ld_root(path, "load_ld")
  manifest <- .read_ld_manifest(file.path(root, "ld_manifest.json"), require_r_reference = TRUE)
  default_chrX_sex <- .validate_chrx_sex(
    if (is.null(default_chrX_sex)) "female" else default_chrX_sex,
    "default_chrX_sex"
  )
  reference_panel <- .load_ld_manifest_reference(root, manifest, shards)
  ref_shards <- shards(reference_panel)
  groups <- .load_ld_shard_groups(root, manifest, ref_shards, retain_ld_r = isTRUE(retain_ld_r))
  .new_ld_panel(groups, default_chrX_sex = default_chrX_sex, reference = reference_panel)
}

load_ld_reference <- function(path, shards = NULL) {
  root <- .validate_ld_root(path, "load_ld_reference")
  manifest <- .read_ld_manifest(file.path(root, "ld_manifest.json"), require_r_reference = TRUE)
  .load_ld_manifest_reference(root, manifest, shards)
}

validate_ld_distribution <- function(path, check_payload_structure = FALSE) {
  root <- .validate_ld_root(path, "validate_ld_distribution")
  manifest <- .read_ld_manifest(file.path(root, "ld_manifest.json"), require_r_reference = TRUE)
  python_reference_cache_path <- file.path(root, manifest$reference_cache)
  .require_ld_file(python_reference_cache_path)
  if (!identical(.md5_file(python_reference_cache_path), manifest$reference_cache_md5)) {
    stop(sprintf("%s: reference_cache_md5 does not match manifest", python_reference_cache_path), call. = FALSE)
  }
  reference_cache_path <- .ld_r_reference_cache_path(root, manifest)
  .require_ld_file(reference_cache_path)
  if (!identical(.md5_file(reference_cache_path), manifest$r_reference_cache_md5)) {
    stop(sprintf("%s: r_reference_cache_md5 does not match manifest", reference_cache_path), call. = FALSE)
  }
  reference_panel <- load_reference_cache(reference_cache_path)
  ref_by_label <- stats::setNames(shards(reference_panel), vapply(shards(reference_panel), function(s) s$label, character(1)))
  seen_reference_bim <- character()

  for (entry in manifest$shards) {
    shard_path <- file.path(root, entry$file)
    .require_ld_file(shard_path)
    if (!identical(.md5_file(shard_path), entry$file_md5)) {
      stop(sprintf("%s: file_md5 does not match manifest", shard_path), call. = FALSE)
    }
    loaded <- .read_ld_npz_shard(
      shard_path,
      check_payload_structure = isTRUE(check_payload_structure),
      retain_ld_r = TRUE,
      expected_file_md5 = entry$file_md5
    )
    .validate_manifest_entry_agreement(entry, loaded$metadata, shard_path)
    ref_shard <- ref_by_label[[entry$chr]]
    if (is.null(ref_shard)) {
      stop(sprintf("%s: reference cache missing shard %s", reference_cache_path, entry$chr), call. = FALSE)
    }
    .validate_ld_reference_entry(ref_shard, entry, reference_cache_path)

    reference_key <- entry$chr
    if (!(reference_key %in% seen_reference_bim)) {
      .validate_bundled_reference_bim(file.path(root, entry$reference_bim), entry)
      seen_reference_bim <- c(seen_reference_bim, reference_key)
    }
  }
  list(ok = TRUE)
}

prepare_ld_npz_for_r <- function(npz_root, extract_npz = FALSE) {
  npz_root <- .validate_ld_root(npz_root, "prepare_ld_npz_for_r")
  if (!is.logical(extract_npz) || length(extract_npz) != 1L || is.na(extract_npz)) {
    stop("extract_npz must be a scalar logical value", call. = FALSE)
  }
  manifest_path <- file.path(npz_root, "ld_manifest.json")
  manifest <- .read_ld_manifest(manifest_path, require_r_reference = FALSE)
  labels <- .canonical_chr_order[.canonical_chr_order %in% unique(vapply(manifest$shards, function(e) e$chr, character(1)))]
  if (!length(labels)) {
    stop("prepare_ld_npz_for_r: manifest has no supported LD shards", call. = FALSE)
  }

  reference_cache <- if (is.null(manifest$r_reference_cache)) {
    "reference_cache.rds"
  } else {
    .validate_plain_relative_filename(manifest$r_reference_cache, sprintf("%s: r_reference_cache", manifest_path))
  }
  reference_panel <- .load_reference_from_ld_bims(npz_root, manifest$shards, labels)
  .validate_ld_manifest_reference_panel(reference_panel, manifest$shards, npz_root)
  reference_cache_path <- file.path(npz_root, reference_cache)
  .save_reference_cache_atomic(reference_panel, reference_cache_path)

  manifest$r_reference_cache <- reference_cache
  manifest$r_reference_cache_md5 <- .md5_file(reference_cache_path)
  if (isTRUE(extract_npz)) {
    for (entry in manifest$shards) {
      .extract_npz_ld_payload_cache(file.path(npz_root, entry$file), entry$file_md5)
    }
  }
  .write_ld_manifest(manifest, manifest_path)
  invisible(manifest)
}

num_snp.LDShard <- function(x, ...) x$num_snp
num_snp.LDPanel <- function(x, ...) x$num_snp
shards.LDPanel <- function(x, ...) .ld_panel_selected_shards(x, NULL)
shard_offsets.LDPanel <- function(x, ...) x$shard_offsets
reference.LDPanel <- function(x, ...) x$reference
default_chrX_sex.LDPanel <- function(x, ...) x$default_chrX_sex

a1freq.LDShard <- function(x, chrX_sex = NULL, ...) x$a1freq
a1freq.LDPanel <- function(x, chrX_sex = NULL, ...) {
  parts <- lapply(.ld_panel_selected_shards(x, chrX_sex), function(s) s$a1freq)
  unlist(parts, use.names = FALSE)
}

select_shards.LDPanel <- function(x, shards, ...) {
  available <- vapply(x$shard_groups, function(g) g[[1]]$label, character(1))
  selected <- .validate_requested_shards(shards, available, "LDPanel.select_shards")
  by_label <- stats::setNames(x$shard_groups, available)
  .new_ld_panel(
    unname(by_label[selected]),
    default_chrX_sex = x$default_chrX_sex,
    reference = if (is.null(x$reference)) NULL else select_shards(x$reference, selected)
  )
}

multiply_r2.LDPanel <- function(ld_panel, M, chrX_sex = NULL, ...) {
  if (inherits(M, "sparseMatrix")) {
    return(.multiply_r2_sparse(ld_panel, M, chrX_sex))
  }
  arr <- as.matrix(M)
  vector_input <- is.null(dim(M))
  if (vector_input) {
    if (length(M) != ld_panel$num_snp) {
      stop(sprintf("multiply_r2: vector length mismatch: expected %d, got %d", ld_panel$num_snp, length(M)), call. = FALSE)
    }
    arr <- matrix(as.numeric(M), ncol = 1L)
  } else if (length(dim(M)) == 2L) {
    if (nrow(arr) != ld_panel$num_snp) {
      stop(sprintf("multiply_r2: matrix row count mismatch: expected %d, got %d", ld_panel$num_snp, nrow(arr)), call. = FALSE)
    }
    storage.mode(arr) <- "numeric"
  } else {
    stop("multiply_r2: M must be a vector or 2D matrix", call. = FALSE)
  }

  out <- matrix(0, nrow = nrow(arr), ncol = ncol(arr))
  for (entry in .ld_panel_selected_shards_with_offsets(ld_panel, chrX_sex)) {
    rows <- seq.int(entry$start0 + 1L, entry$stop0)
    out[rows, ] <- as.matrix(entry$shard$ld_r2 %*% arr[rows, , drop = FALSE])
  }
  if (vector_input) {
    return(as.numeric(out[, 1L]))
  }
  out
}

fast_prune <- function(logpvec, ld_panel, r2_threshold = 0.2, chrX_sex = NULL) {
  if (!inherits(ld_panel, "LDPanel")) {
    stop("fast_prune: ld_panel must be an LDPanel", call. = FALSE)
  }
  threshold <- as.numeric(r2_threshold)
  if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) || threshold < 0) {
    stop("fast_prune: r2_threshold must be a finite non-negative number", call. = FALSE)
  }
  dims <- dim(logpvec)
  if (!is.null(dims) && !(length(dims) == 2L && (dims[[1L]] == 1L || dims[[2L]] == 1L))) {
    stop("fast_prune: logpvec must be a vector", call. = FALSE)
  }
  vec <- as.numeric(logpvec)
  if (length(vec) != ld_panel$num_snp) {
    stop(sprintf("fast_prune: vector length mismatch: expected %d, got %d", ld_panel$num_snp, length(vec)), call. = FALSE)
  }
  out <- vec
  for (entry in .ld_panel_selected_shards_with_offsets(ld_panel, chrX_sex)) {
    rows <- seq.int(entry$start0 + 1L, entry$stop0)
    out[rows] <- .fast_prune_shard(out[rows], entry$shard$ld_r2, threshold)
  }
  if (!is.null(dims)) {
    return(matrix(out, nrow = dims[[1L]], ncol = dims[[2L]]))
  }
  out
}

print.LDPanel <- function(x, ...) {
  labels <- vapply(x$shard_groups, function(g) g[[1]]$label, character(1))
  cat(sprintf("<LDPanel: %d SNPs across %d shard group(s): %s>\n", x$num_snp, length(labels), paste(labels, collapse = ", ")))
  invisible(x)
}

.new_ld_shard <- function(chr, sex, num_snp, ld_r, a1freq, reference_checksum,
                          num_monomorphic_snps = 0L, retain_ld_r = TRUE) {
  if (!inherits(ld_r, "dgCMatrix")) {
    ld_r <- methods::as(ld_r, "dgCMatrix")
  }
  ld_r2 <- .build_ld_r2(ld_r)
  structure(
    list(
      label = as.character(chr),
      chr = as.character(chr),
      sex = sex,
      num_snp = as.integer(num_snp),
      ld_r = if (isTRUE(retain_ld_r)) ld_r else NULL,
      ld_r2 = ld_r2,
      a1freq = as.numeric(a1freq),
      reference_checksum = as.character(reference_checksum),
      num_monomorphic_snps = as.integer(num_monomorphic_snps)
    ),
    class = "LDShard"
  )
}

.new_ld_panel <- function(shard_groups, default_chrX_sex = "female", reference = NULL) {
  if (!length(shard_groups)) {
    stop("LDPanel requires at least one shard group", call. = FALSE)
  }
  default_chrX_sex <- .validate_chrx_sex(default_chrX_sex, "default_chrX_sex")
  pos <- 0L
  offsets <- data.frame(
    shard_label = character(length(shard_groups)),
    start0 = integer(length(shard_groups)),
    stop0 = integer(length(shard_groups)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(shard_groups)) {
    group <- shard_groups[[i]]
    if (!length(group)) {
      stop("LDPanel shard groups must be non-empty", call. = FALSE)
    }
    first <- group[[1L]]
    for (shard in group) {
      if (!identical(shard$chr, first$chr)) {
        stop("all LD shards in a group must share chr", call. = FALSE)
      }
      if (!identical(as.integer(shard$num_snp), as.integer(first$num_snp))) {
        stop("all LD shards in a group must share num_snp", call. = FALSE)
      }
    }
    if (identical(first$chr, "X")) {
      present <- vapply(group, function(s) s$sex, character(1))
      if (!(default_chrX_sex %in% present)) {
        stop(sprintf(
          "default_chrX_sex must name a loaded chrX LD shard; got %s, present %s",
          sQuote(default_chrX_sex), paste(sQuote(sort(present)), collapse = ", ")
        ), call. = FALSE)
      }
    }
    offsets$shard_label[[i]] <- first$label
    offsets$start0[[i]] <- pos
    pos <- pos + first$num_snp
    offsets$stop0[[i]] <- pos
  }
  if (!is.null(reference)) {
    ref_shards <- shards(reference)
    if (length(ref_shards) != length(shard_groups)) {
      stop("LDPanel reference shard count must match LD shard groups", call. = FALSE)
    }
    for (i in seq_along(ref_shards)) {
      if (!identical(ref_shards[[i]]$label, shard_groups[[i]][[1]]$label) ||
          !identical(as.integer(ref_shards[[i]]$num_snp), as.integer(shard_groups[[i]][[1]]$num_snp))) {
        stop("LDPanel reference shards must match LD shard groups", call. = FALSE)
      }
    }
  }
  structure(
    list(
      shard_groups = shard_groups,
      num_snp = pos,
      shard_offsets = offsets,
      default_chrX_sex = default_chrX_sex,
      reference = reference
    ),
    class = "LDPanel"
  )
}

.build_ld_r2 <- function(ld_r) {
  .new_dgCMatrix_from_slots(
    i = ld_r@i,
    p = ld_r@p,
    x = ld_r@x * ld_r@x,
    dim = ld_r@Dim,
    dimnames = ld_r@Dimnames
  )
}

.multiply_r2_sparse <- function(ld_panel, M, chrX_sex) {
  if (length(dim(M)) != 2L) {
    stop("multiply_r2: sparse M must be a 2D matrix", call. = FALSE)
  }
  if (nrow(M) != ld_panel$num_snp) {
    stop(sprintf("multiply_r2: matrix row count mismatch: expected %d, got %d", ld_panel$num_snp, nrow(M)), call. = FALSE)
  }
  parts <- vector("list", length(ld_panel$shard_groups))
  k <- 0L
  for (entry in .ld_panel_selected_shards_with_offsets(ld_panel, chrX_sex)) {
    k <- k + 1L
    rows <- seq.int(entry$start0 + 1L, entry$stop0)
    parts[[k]] <- entry$shard$ld_r2 %*% M[rows, , drop = FALSE]
  }
  do.call(rbind, parts)
}

.fast_prune_shard <- function(values, ld_r2, threshold) {
  out <- as.numeric(values)
  finite_idx <- which(is.finite(out))
  if (!length(finite_idx)) {
    return(out)
  }
  order_idx <- finite_idx[order(-abs(out[finite_idx]), finite_idx)]
  pruned <- rep.int(FALSE, length(out))
  retained <- rep.int(FALSE, length(out))
  # PERF: loop retained because greedy pruning is order-dependent; vectorization
  # is not used because each retained SNP changes which later SNPs are eligible.
  for (idx in order_idx) {
    if (pruned[[idx]]) {
      next
    }
    retained[[idx]] <- TRUE
    start <- ld_r2@p[[idx]] + 1L
    stop <- ld_r2@p[[idx + 1L]]
    if (start > stop) {
      next
    }
    row_idx <- ld_r2@i[start:stop] + 1L
    values_r2 <- ld_r2@x[start:stop]
    neighbors <- row_idx[values_r2 >= threshold]
    neighbors <- neighbors[neighbors != idx & !retained[neighbors]]
    if (length(neighbors)) {
      pruned[neighbors] <- TRUE
      out[neighbors] <- NaN
    }
  }
  out
}

.ld_panel_selected_shards <- function(panel, chrX_sex) {
  lapply(.ld_panel_selected_shards_with_offsets(panel, chrX_sex), function(entry) entry$shard)
}

.ld_panel_selected_shards_with_offsets <- function(panel, chrX_sex) {
  out <- vector("list", length(panel$shard_groups))
  for (i in seq_along(panel$shard_groups)) {
    group <- panel$shard_groups[[i]]
    first <- group[[1L]]
    if (!identical(first$chr, "X")) {
      selected <- first
    } else {
      sex <- if (is.null(chrX_sex)) panel$default_chrX_sex else .validate_chrx_sex(chrX_sex, "chrX_sex")
      by_sex <- stats::setNames(group, vapply(group, function(s) s$sex, character(1)))
      selected <- by_sex[[sex]]
      if (is.null(selected)) {
        present <- names(by_sex)
        stop(sprintf(
          "chrX_sex must name a loaded chrX LD shard; got %s, present %s",
          sQuote(sex), paste(sQuote(sort(present)), collapse = ", ")
        ), call. = FALSE)
      }
    }
    out[[i]] <- list(
      shard = selected,
      start0 = panel$shard_offsets$start0[[i]],
      stop0 = panel$shard_offsets$stop0[[i]]
    )
  }
  out
}

.load_ld_shard_groups <- function(root, manifest, ref_shards, retain_ld_r) {
  entries <- manifest$shards
  groups <- vector("list", length(ref_shards))
  for (i in seq_along(ref_shards)) {
    ref_shard <- ref_shards[[i]]
    if (identical(ref_shard$label, "X")) {
      selected <- Filter(function(e) identical(e$chr, "X"), entries)
      if (!length(selected)) {
        stop("LD manifest has no chrX shard for reference shard X", call. = FALSE)
      }
    } else {
      selected <- Filter(function(e) identical(e$chr, ref_shard$label) && is.null(e$sex), entries)
      if (length(selected) != 1L) {
        stop(sprintf("LD manifest must contain exactly one autosomal shard for chr %s", ref_shard$label), call. = FALSE)
      }
    }
    group <- vector("list", length(selected))
    seen_sex <- character()
    for (j in seq_along(selected)) {
      entry <- selected[[j]]
      shard_path <- file.path(root, entry$file)
      suffix <- if (is.null(entry$sex)) "" else sprintf(" (%s)", entry$sex)
      .ld_info(sprintf("statgen.load_ld: loading shard %s%s from %s", entry$chr, suffix, shard_path))
      loaded <- .read_ld_npz_shard(
        shard_path,
        check_payload_structure = FALSE,
        retain_ld_r = retain_ld_r,
        expected_file_md5 = entry$file_md5
      )
      .validate_manifest_entry_agreement(entry, loaded$metadata, shard_path)
      .validate_ld_reference_compatibility(loaded$shard, ref_shard, shard_path)
      key <- if (is.null(loaded$shard$sex)) "<null>" else loaded$shard$sex
      if (key %in% seen_sex) {
        stop(sprintf("LD manifest has duplicate shard for chr %s sex %s", loaded$shard$chr, sQuote(key)), call. = FALSE)
      }
      seen_sex <- c(seen_sex, key)
      group[[j]] <- loaded$shard
    }
    groups[[i]] <- group
  }
  groups
}

.ld_info <- function(message_text) {
  if (identical(get_verbosity(), "info")) {
    message(message_text)
  }
  invisible(NULL)
}
