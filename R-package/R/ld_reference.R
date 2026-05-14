.load_ld_manifest_reference <- function(root, manifest, shards) {
  reference_cache_path <- .ld_r_reference_cache_path(root, manifest)
  .require_ld_file(reference_cache_path)
  load_reference_cache(reference_cache_path, shards = shards)
}

.ld_r_reference_cache_path <- function(root, manifest) {
  file.path(root, .validate_plain_relative_filename(manifest$r_reference_cache, "r_reference_cache"))
}

.save_reference_cache_atomic <- function(panel, path) {
  tmp <- tempfile(".reference_cache_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  save_reference_cache(panel, tmp)
  if (!file.rename(tmp, path)) {
    stop(sprintf("Failed to replace reference cache: %s", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_ld_reference_compatibility <- function(shard, ref_shard, path) {
  if (!identical(shard$chr, ref_shard$label)) {
    stop(sprintf("%s: LD chr %s does not match reference shard %s", path, sQuote(shard$chr), sQuote(ref_shard$label)), call. = FALSE)
  }
  if (!identical(as.integer(shard$num_snp), as.integer(ref_shard$num_snp))) {
    stop(sprintf("%s: LD num_snp does not match reference shard", path), call. = FALSE)
  }
  if (!identical(shard$reference_checksum, ref_shard$checksum)) {
    stop(sprintf("%s: LD reference_checksum does not match reference shard", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_ld_reference_entry <- function(ref_shard, entry, path) {
  if (!identical(as.integer(ref_shard$num_snp), as.integer(entry$num_snp))) {
    stop(sprintf("%s: reference cache num_snp does not match manifest", path), call. = FALSE)
  }
  if (!identical(ref_shard$checksum, entry$reference_checksum)) {
    stop(sprintf("%s: reference cache reference_checksum does not match manifest", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_bundled_reference_bim <- function(path, entry) {
  .require_ld_file(path)
  bim <- .parse_bim(path)
  keep <- bim$chr %in% .canonical_chr_order
  bim <- bim[keep, , drop = FALSE]
  if (!nrow(bim)) {
    stop(sprintf("%s: bundled reference BIM has no supported rows", path), call. = FALSE)
  }
  if (any(bim$chr != entry$chr)) {
    stop(sprintf("%s: bundled reference BIM chr does not match manifest", path), call. = FALSE)
  }
  checksum <- .checksum_from_arrays(bim$chr, bim$bp, bim$a1, bim$a2)
  if (!identical(nrow(bim), as.integer(entry$num_snp))) {
    stop(sprintf("%s: bundled reference BIM num_snp does not match manifest", path), call. = FALSE)
  }
  if (!identical(checksum, entry$reference_checksum)) {
    stop(sprintf("%s: bundled reference BIM checksum does not match manifest", path), call. = FALSE)
  }
  invisible(NULL)
}

.validate_ld_manifest_reference_panel <- function(reference_panel, entries, root) {
  by_label <- stats::setNames(shards(reference_panel), vapply(shards(reference_panel), function(s) s$label, character(1)))
  seen <- character()
  for (entry in entries) {
    ref_shard <- by_label[[entry$chr]]
    if (is.null(ref_shard)) {
      stop(sprintf("%s: reference cache missing shard %s", root, entry$chr), call. = FALSE)
    }
    .validate_ld_reference_entry(ref_shard, entry, root)
    if (!(entry$chr %in% seen)) {
      .require_ld_file(file.path(root, entry$reference_bim))
      seen <- c(seen, entry$chr)
    }
  }
  invisible(NULL)
}

.load_reference_from_ld_bims <- function(root, entries, selected_labels) {
  by_chr <- list()
  for (entry in entries) {
    if (is.null(by_chr[[entry$chr]])) {
      by_chr[[entry$chr]] <- entry$reference_bim
    }
  }
  out <- vector("list", length(selected_labels))
  for (i in seq_along(selected_labels)) {
    label <- selected_labels[[i]]
    bim_path <- file.path(root, by_chr[[label]])
    .require_ld_file(bim_path)
    bim <- .parse_bim(bim_path)
    keep <- bim$chr == label
    out[[i]] <- .new_reference_shard(
      label, bim$chr[keep], bim$snp[keep], bim$bp[keep], bim$a1[keep], bim$a2[keep],
      bim$a1_hash64[keep], bim$a2_hash64[keep]
    )
  }
  .new_reference_panel(out)
}
