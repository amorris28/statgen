.read_npz_ld_payload <- function(path, scratch_parent = NULL, expected_file_md5 = NULL) {
  extracted_dir <- .ld_npz_extracted_dir(path)
  if (.is_valid_extracted_npz_cache(path, extracted_dir, expected_file_md5)) {
    return(.read_extracted_npz_ld_payload(path, extracted_dir))
  }

  listing <- utils::unzip(path, list = TRUE)
  names_in_zip <- listing$Name
  required_files <- paste0(.ld_required_arrays, ".npy")
  missing <- setdiff(required_files, names_in_zip)
  if (length(missing)) {
    stop(sprintf("%s: missing required arrays: %s", path, paste(sub("\\.npy$", "", missing), collapse = ", ")), call. = FALSE)
  }
  if (is.null(scratch_parent)) {
    scratch_parent <- .npz_scratch_parent(path)
  }
  dir.create(scratch_parent, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("statgen_npz_", tmpdir = scratch_parent)
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(path, files = required_files, exdir = tmp)
  arrays <- lapply(.ld_required_arrays, function(name) .read_npy_array(file.path(tmp, paste0(name, ".npy")), path, name))
  .npz_arrays_to_payload(arrays, path)
}

.ld_npz_extracted_dir <- function(path) {
  paste0(path, ".d")
}

.ld_npz_extracted_marker <- function(cache_dir) {
  file.path(cache_dir, ".statgen-extracted-npz.json")
}

.is_valid_extracted_npz_cache <- function(path, cache_dir, expected_file_md5 = NULL) {
  if (!dir.exists(cache_dir)) {
    return(FALSE)
  }
  required_files <- file.path(cache_dir, paste0(.ld_required_arrays, ".npy"))
  if (any(!file.exists(required_files) | dir.exists(required_files))) {
    return(FALSE)
  }
  marker_path <- .ld_npz_extracted_marker(cache_dir)
  if (!file.exists(marker_path) || dir.exists(marker_path)) {
    return(FALSE)
  }
  marker <- tryCatch(
    jsonlite::fromJSON(marker_path, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (!is.list(marker) ||
      !identical(marker$format, "statgen_extracted_npz") ||
      !identical(marker$source_file, basename(path))) {
    return(FALSE)
  }
  source_md5 <- if (is.null(expected_file_md5)) .md5_file(path) else expected_file_md5
  identical(marker$source_file_md5, source_md5)
}

.read_extracted_npz_ld_payload <- function(path, cache_dir) {
  arrays <- lapply(.ld_required_arrays, function(name) .read_npy_array(file.path(cache_dir, paste0(name, ".npy")), path, name))
  .npz_arrays_to_payload(arrays, path)
}

.npz_arrays_to_payload <- function(arrays, path) {
  names(arrays) <- .ld_required_arrays
  metadata_raw <- as.raw(arrays$metadata$data)
  metadata <- tryCatch(
    jsonlite::fromJSON(rawToChar(metadata_raw), simplifyVector = FALSE),
    error = function(e) stop(sprintf("%s: metadata must be UTF-8 JSON bytes", path), call. = FALSE)
  )
  list(
    data = arrays$data$data,
    indices = arrays$indices$data,
    indptr = arrays$indptr$data,
    shape = arrays$shape$data,
    a1freq = arrays$a1freq$data,
    metadata = metadata
  )
}

.extract_npz_ld_payload_cache <- function(path, expected_file_md5) {
  .require_ld_file(path)
  source_md5 <- .md5_file(path)
  if (!identical(source_md5, expected_file_md5)) {
    stop(sprintf("%s: file_md5 does not match manifest", path), call. = FALSE)
  }
  listing <- utils::unzip(path, list = TRUE)
  required_files <- paste0(.ld_required_arrays, ".npy")
  missing <- setdiff(required_files, listing$Name)
  if (length(missing)) {
    stop(sprintf("%s: missing required arrays: %s", path, paste(sub("\\.npy$", "", missing), collapse = ", ")), call. = FALSE)
  }

  cache_dir <- .ld_npz_extracted_dir(path)
  tmp <- tempfile(".statgen_npz_d_", tmpdir = dirname(path))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(path, files = required_files, exdir = tmp)
  extracted_files <- file.path(tmp, required_files)
  if (any(!file.exists(extracted_files) | dir.exists(extracted_files))) {
    stop(sprintf("%s: failed to extract required LD arrays", path), call. = FALSE)
  }
  marker <- list(
    format = "statgen_extracted_npz",
    schema_version = "1.0",
    source_file = basename(path),
    source_file_md5 = source_md5
  )
  marker_text <- jsonlite::toJSON(marker, auto_unbox = TRUE, pretty = TRUE)
  writeLines(marker_text, .ld_npz_extracted_marker(tmp), useBytes = TRUE)

  if (file.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
  }
  if (!file.rename(tmp, cache_dir)) {
    stop(sprintf("Failed to replace extracted LD cache: %s", cache_dir), call. = FALSE)
  }
  invisible(cache_dir)
}

.npz_scratch_parent <- function(path) {
  env_root <- Sys.getenv("STATGEN_SCRATCH", unset = "")
  if (nzchar(env_root)) {
    dir.create(env_root, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(env_root)) {
      return(env_root)
    }
    stop(sprintf("Cannot create STATGEN_SCRATCH directory: %s", env_root), call. = FALSE)
  }
  root <- dirname(path)
  if (dir.exists(root) && file.access(root, 2L) == 0L) {
    return(root)
  }
  warning("STATGEN_SCRATCH is unset and LD .npz directory is not writable; using system temporary directory for extraction", call. = FALSE)
  tempdir()
}

.read_npy_array <- function(file, archive_path, array_name) {
  con <- file(file, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 6L)
  if (!identical(magic, as.raw(c(0x93, charToRaw("NUMPY"))))) {
    stop(sprintf("%s: %s is not an NPY array", archive_path, array_name), call. = FALSE)
  }
  version <- readBin(con, "integer", n = 2L, size = 1L, signed = FALSE, endian = "little")
  if (length(version) != 2L || !(version[[1L]] %in% c(1L, 2L, 3L))) {
    stop(sprintf("%s: unsupported NPY version for %s", archive_path, array_name), call. = FALSE)
  }
  header_len <- if (version[[1L]] == 1L) {
    readBin(con, "integer", n = 1L, size = 2L, signed = FALSE, endian = "little")
  } else {
    readBin(con, "integer", n = 1L, size = 4L, signed = FALSE, endian = "little")
  }
  header_raw <- readBin(con, "raw", n = header_len)
  header <- rawToChar(header_raw)
  descr <- .parse_npy_header_value(header, "descr")
  fortran <- .parse_npy_header_value(header, "fortran_order")
  shape <- .parse_npy_shape(header)
  if (!identical(fortran, "False")) {
    stop(sprintf("%s: %s must not be Fortran-ordered", archive_path, array_name), call. = FALSE)
  }
  expected_descr <- switch(
    array_name,
    data = "<f4",
    indices = "<i4",
    indptr = "<i4",
    shape = "<i8",
    a1freq = "<f4",
    metadata = "|u1",
    stop(sprintf("%s: unsupported array %s", archive_path, array_name), call. = FALSE)
  )
  if (!identical(descr, expected_descr)) {
    stop(sprintf("%s: %s must have dtype %s, got %s", archive_path, array_name, expected_descr, descr), call. = FALSE)
  }
  count <- prod(shape)
  if (!is.finite(count) || count < 0 || count > .Machine$integer.max) {
    stop(sprintf("%s: %s shape is unsupported", archive_path, array_name), call. = FALSE)
  }
  data <- switch(
    descr,
    "<f4" = readBin(con, "numeric", n = as.integer(count), size = 4L, endian = "little"),
    "<i4" = readBin(con, "integer", n = as.integer(count), size = 4L, signed = TRUE, endian = "little"),
    "|u1" = readBin(con, "integer", n = as.integer(count), size = 1L, signed = FALSE, endian = "little"),
    "<i8" = .read_npy_int64_as_integer(con, as.integer(count), archive_path, array_name)
  )
  if (length(data) != count) {
    stop(sprintf("%s: short read for %s", archive_path, array_name), call. = FALSE)
  }
  list(data = data, shape = shape, descr = descr)
}

.parse_npy_header_value <- function(header, key) {
  pattern <- sprintf("'%s': ([^,}]+)", key)
  hit <- regmatches(header, regexpr(pattern, header, perl = TRUE))
  if (!length(hit) || hit == "") {
    stop(sprintf("Malformed NPY header: missing %s", key), call. = FALSE)
  }
  value <- sub(pattern, "\\1", hit, perl = TRUE)
  value <- trimws(value)
  sub("^'([^']*)'$", "\\1", value)
}

.parse_npy_shape <- function(header) {
  hit <- regmatches(header, regexpr("'shape': \\(([^)]*)\\)", header, perl = TRUE))
  if (!length(hit) || hit == "") {
    stop("Malformed NPY header: missing shape", call. = FALSE)
  }
  inside <- sub(".*'shape': \\(([^)]*)\\).*", "\\1", hit, perl = TRUE)
  parts <- trimws(strsplit(inside, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    return(integer())
  }
  out <- suppressWarnings(as.integer(parts))
  if (anyNA(out) || any(out < 0L)) {
    stop("Malformed NPY header: invalid shape", call. = FALSE)
  }
  out
}

.read_npy_int64_as_integer <- function(con, n, archive_path, array_name) {
  if (!n) {
    return(integer())
  }
  raw <- readBin(con, "raw", n = n * 8L)
  if (length(raw) != n * 8L) {
    stop(sprintf("%s: short read for %s", archive_path, array_name), call. = FALSE)
  }
  bytes <- matrix(as.integer(raw), nrow = 8L)
  negative <- bytes[8L, ] >= 128L
  high <- bytes[5L, ] != 0L | bytes[6L, ] != 0L | bytes[7L, ] != 0L | bytes[8L, ] != 0L
  exceeds_r_integer <- bytes[4L, ] >= 128L
  if (any(negative | high | exceeds_r_integer)) {
    stop(sprintf("%s: %s int64 values exceed the supported R integer range", archive_path, array_name), call. = FALSE)
  }
  values <- bytes[1L, ] + bytes[2L, ] * 256L + bytes[3L, ] * 65536L + bytes[4L, ] * 16777216L
  as.integer(values)
}
