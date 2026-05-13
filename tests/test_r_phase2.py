"""R phase 2 acceptance tests: sumstats and annotations."""

import json
import math

import numpy as np
import pytest

from statgen.annotations import load_annotations
from statgen.reference import load_reference
from statgen.sumstats import load_sumstats
from tests.conftest import FIXTURES_DIR, R_PACKAGE_DIR, run_rscript, skipif_no_rscript

SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
SUMSTATS = FIXTURES_DIR / "sumstats/traits.tsv.gz"
ANNOTATIONS = [FIXTURES_DIR / "annotations/anno1.bed", FIXTURES_DIR / "annotations/anno2.bed"]


def _source_phase2_script(expr: str) -> str:
    files = [
        "utils.R",
        "verbosity.R",
        "hash.R",
        "variant_match.R",
        "bfile_utils.R",
        "reference.R",
        "sumstats.R",
        "annotations.R",
    ]
    sources = " ".join(
        f"source({json.dumps(str(R_PACKAGE_DIR / 'R' / filename))});"
        for filename in files
    )
    return f"{sources} {expr}"


def _parse_float_csv(text: str) -> list[float]:
    out = []
    for token in text.split(","):
        if token == "NaN":
            out.append(math.nan)
        elif token == "Inf":
            out.append(math.inf)
        else:
            out.append(float(token))
    return out


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_matches_python_fixture_alignment():
    reference = load_reference(SHARDED_REF)
    with pytest.warns(RuntimeWarning):
        py = load_sumstats(SUMSTATS, reference)

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(SUMSTATS))}, ref)); "
            "cat(num_snp(s), '\\n'); "
            "cat(paste(vapply(shards(s), function(x) x$label, character(1)), collapse=','), '\\n'); "
            "cat(paste(format(logpvec(s), digits=17), collapse=','), '\\n'); "
            "cat(paste(format(zvec(s), digits=17), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_present(s)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "8"
    assert lines[1] == "1,X"
    np.testing.assert_allclose(_parse_float_csv(lines[2]), py.logpvec, equal_nan=True)
    np.testing.assert_allclose(_parse_float_csv(lines[3]), py.zvec, equal_nan=True)
    assert lines[4] == ",".join(str(int(x)) for x in py.is_present)


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_cache_subset_and_create(tmp_path):
    cache = tmp_path / "sumstats.rds"
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(SUMSTATS))}, ref)); "
            f"save_sumstats_cache(s, {json.dumps(str(cache))}); "
            f"loaded <- suppressWarnings(load_sumstats_cache({json.dumps(str(cache))}, shards = 'X')); "
            "created <- create_sumstats(ref, p = rep(0.5, num_snp(ref)), z = rep(1, num_snp(ref)), n = rep(100, num_snp(ref)), beta = rep(0.1, num_snp(ref)), se = rep(0.2, num_snp(ref)), eaf = rep(0.3, num_snp(ref)), info = rep(0.9, num_snp(ref))); "
            "cat(num_snp(loaded), '\\n'); "
            "cat(paste(format(logpvec(loaded), digits=17), collapse=','), '\\n'); "
            "cat(paste(as.integer(is_present(created)), collapse=','), '\\n'); "
            "cat(paste(format(nvec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(beta_vec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(se_vec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(eaf_vec(created)[1:3], digits=17), collapse=','), '\\n'); "
            "cat(paste(format(info_vec(created)[1:3], digits=17), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "3"
    np.testing.assert_allclose(
        _parse_float_csv(lines[1]),
        np.array([-math.log10(0.003), -math.log10(0.6), math.nan]),
        equal_nan=True,
    )
    assert lines[2] == "1,1,1,1,1,1,1,1"
    np.testing.assert_allclose(_parse_float_csv(lines[3]), np.array([100, 100, 100]))
    np.testing.assert_allclose(_parse_float_csv(lines[4]), np.array([0.1, 0.1, 0.1]))
    np.testing.assert_allclose(_parse_float_csv(lines[5]), np.array([0.2, 0.2, 0.2]))
    np.testing.assert_allclose(_parse_float_csv(lines[6]), np.array([0.3, 0.3, 0.3]))
    np.testing.assert_allclose(_parse_float_csv(lines[7]), np.array([0.9, 0.9, 0.9]))


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_constructor_rejects_inconsistent_optional_fields():
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            "s <- create_sumstats(ref, p = rep(0.5, num_snp(ref)), z = rep(1, num_snp(ref))); "
            "bad <- s; bad$shards[[2]]$zvec <- NULL; "
            "ok <- FALSE; tryCatch(.new_sumstats_panel(bad$shards), error = function(e) ok <<- grepl('consistent zvec presence', e$message)); "
            "ok2 <- FALSE; tryCatch(zvec(bad), error = function(e) ok2 <<- grepl('inconsistent zvec presence', e$message)); "
            "stopifnot(ok, ok2)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_annotations_match_python_fixture_matrix_and_cache(tmp_path):
    reference = load_reference(SHARDED_REF)
    py = load_annotations(ANNOTATIONS, reference)
    expected_col_major = py.annomat.toarray().astype(int).ravel(order="F")
    cache = tmp_path / "annotations.rds"

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"a <- load_annotations(c({json.dumps(str(ANNOTATIONS[0]))}, {json.dumps(str(ANNOTATIONS[1]))}), ref); "
            f"save_annotations_cache(a, {json.dumps(str(cache))}); "
            f"loaded <- load_annotations_cache({json.dumps(str(cache))}); "
            "selected <- select_annotations(a, c('anno2', 'anno1')); "
            "cat(num_snp(a), '\\n'); "
            "cat(num_annot(a), '\\n'); "
            "cat(paste(annonames(a), collapse=','), '\\n'); "
            "cat(class(annomat(a))[[1]], '\\n'); "
            "cat(paste(as.integer(as.vector(as.matrix(annomat(loaded)))), collapse=','), '\\n'); "
            "cat(paste(annonames(selected), collapse=','), '\\n'); "
            "cat(paste(colnames(annomat(selected)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[:4] == ["8", "2", "anno1,anno2", "lgCMatrix"]
    # R stores matrices in column-major order; this flattening matches the
    # Python expected_col_major vector built from the dense annotation matrix.
    assert lines[4] == ",".join(str(x) for x in expected_col_major.tolist())
    assert lines[5] == "anno2,anno1"
    assert lines[6] == "anno2,anno1"


@pytest.mark.r
@skipif_no_rscript
def test_r_annotation_union_uses_reference_compatible_panels():
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            "a <- create_annotation(ref, rep(c(1, 0), length.out = num_snp(ref)), 'a'); "
            "b <- create_annotation(ref, rep(c(0, 1), length.out = num_snp(ref)), 'b'); "
            "u <- union_annotations(a, b); "
            "bad <- b; bad$shards[[1]]$reference_checksum <- 'bad'; "
            "warned <- FALSE; ok <- FALSE; withCallingHandlers(tryCatch(union_annotations(a, bad), error = function(e) ok <<- grepl('compatible reference alignment', e$message)), warning = function(w) warned <<- TRUE); "
            "msg <- NULL; tryCatch(union_annotations(a, bad), error = function(e) msg <<- e$message); "
            "cat(paste(annonames(u), collapse=','), '\\n'); "
            "cat(paste(colnames(annomat(u)), collapse=','), '\\n'); "
            "cat(as.character(ok), '\\n'); "
            "cat(as.character(warned), '\\n'); "
            "cat(as.character(grepl('reference_checksum mismatch', msg)), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["a,b", "a,b", "TRUE", "FALSE", "TRUE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_annotations_allow_leading_comments_but_reject_late_comments(tmp_path):
    leading = tmp_path / "leading.bed"
    leading.write_text("# header\n\n1\t99\t200\n")
    late = tmp_path / "late.bed"
    late.write_text("1\t99\t200\n# late\n")

    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}, shards = '1'); "
            f"a <- load_annotations({json.dumps(str(leading))}, ref); "
            f"ok <- FALSE; tryCatch(load_annotations({json.dumps(str(late))}, ref), error = function(e) ok <<- grepl('blank or comment line after BED data row', e$message)); "
            "cat(num_snp(a), num_annot(a), '\\n'); "
            "cat(paste(as.integer(as.vector(as.matrix(annomat(a)))), collapse=','), '\\n'); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines == ["5 1", "1,1,0,0,0", "TRUE"]


@pytest.mark.r
@skipif_no_rscript
def test_r_phase2_rejects_malformed_inputs(tmp_path):
    bad_sumstats = tmp_path / "bad.tsv"
    bad_sumstats.write_text("chr\tbp\ta1\ta2\tp\n1\t100\tA\tG\t1.5\n")
    bad_nan_sumstats = tmp_path / "bad_nan.tsv"
    bad_nan_sumstats.write_text("chr\tbp\ta1\ta2\tp\n1\t100\tA\tG\tNaN\n")
    bad_bed = tmp_path / "bad.bed"
    bad_bed.write_text("1\t100\n")
    bad_cache = tmp_path / "bad_sumstats_cache.rds"
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}); "
            f"saveRDS(list(metadata = list(schema = 'sumstats_cache/0.1', n_shards = 0L, shard_labels = character(), shard_checksums = character(), shard_start0 = numeric(), shard_stop0 = numeric()), logpvec = numeric()), {json.dumps(str(bad_cache))}); "
            f"ok1 <- FALSE; tryCatch(load_sumstats({json.dumps(str(bad_sumstats))}, ref), error = function(e) ok1 <<- grepl('p must be finite numeric', e$message)); "
            f"ok1_nan <- FALSE; tryCatch(load_sumstats({json.dumps(str(bad_nan_sumstats))}, ref), error = function(e) ok1_nan <<- grepl('p must be finite numeric', e$message)); "
            f"ok2 <- FALSE; tryCatch(load_annotations({json.dumps(str(bad_bed))}, ref), error = function(e) ok2 <<- grepl('at least 3', e$message)); "
            "ok3 <- FALSE; tryCatch(load_annotations('', ref), error = function(e) ok3 <<- grepl('bed_paths must be a non-empty character vector', e$message)); "
            "ok4 <- FALSE; tryCatch(load_annotations(NA_character_, ref), error = function(e) ok4 <<- grepl('bed_paths must be a non-empty character vector', e$message)); "
            f"ok5 <- FALSE; tryCatch(load_sumstats_cache({json.dumps(str(bad_cache))}), error = function(e) ok5 <<- grepl('n_shards must be at least 1', e$message)); "
            "stopifnot(ok1, ok1_nan, ok2, ok3, ok4, ok5)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_sumstats_swapped_allele_warning_and_duplicate_key_error(tmp_path):
    swapped = tmp_path / "swapped.tsv"
    swapped.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\t0.01\n"
        "1\t200\tT\tC\t0.02\n"
    )
    duplicate = tmp_path / "duplicate.tsv"
    duplicate.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tA\tG\t0.01\n"
        "1\t100\tA\tG\t0.02\n"
    )
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(SHARDED_REF))}, shards = '1'); "
            "msg <- NULL; "
            f"s <- withCallingHandlers(load_sumstats({json.dumps(str(swapped))}, ref), warning = function(w) {{ msg <<- paste(msg, conditionMessage(w), sep='|'); invokeRestart('muffleWarning') }}); "
            "stopifnot(grepl('would match the reference if a1/a2 were swapped', msg)); "
            "stopifnot(identical(as.integer(is_present(s)), c(1L, 0L, 0L, 0L, 0L))); "
            f"ok <- FALSE; tryCatch(load_sumstats({json.dumps(str(duplicate))}, ref), error = function(e) ok <<- grepl('duplicate sumstats matching key', e$message)); "
            "stopifnot(ok)"
        )
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_variant_matching_preserves_integer64_hashes_above_double_precision(tmp_path):
    ref_path = tmp_path / "high_hash.bim"
    ref_path.write_text(
        "1\trsHigh1\t0\t100\tACGT\tCCCCCCCC\n"
        "1\trsHigh2\t0\t200\tTGCA\tGGGGGGGG\n"
    )
    sumstats_path = tmp_path / "high_hash.tsv"
    sumstats_path.write_text(
        "chr\tbp\ta1\ta2\tp\n"
        "1\t100\tACGT\tCCCCCCCC\t0.01\n"
        "1\t200\tTGCA\tGGGGGGGG\t0.2\n"
    )
    result = run_rscript(
        _source_phase2_script(
            f"ref <- load_reference({json.dumps(str(ref_path))}); "
            "limit <- bit64::as.integer64('9007199254740992'); "
            "stopifnot(all(a1_hash64(ref) > limit), all(a2_hash64(ref) > limit)); "
            f"s <- suppressWarnings(load_sumstats({json.dumps(str(sumstats_path))}, ref)); "
            "h <- .allele_hash64(c('ACGT', 'TGCA', 'CCCCCCCC', 'GGGGGGGG')); "
            "direct <- .match_variant_keys(c(100L, 200L), h[1:2], h[3:4], c(200L, 100L), h[c(2, 1)], h[c(4, 3)], '1', 'sumstats'); "
            "dup_ok <- FALSE; tryCatch(.count_variant_key_intersections(c(100L), h[1], h[3], c(100L, 100L), h[c(1, 1)], h[c(3, 3)]), error = function(e) dup_ok <<- grepl('duplicate variant matching key', e$message)); "
            "cat(paste(as.integer(is_present(s)), collapse=','), '\\n'); "
            "cat(paste(format(logpvec(s), digits=17), collapse=','), '\\n'); "
            "cat(paste(direct, collapse=','), '\\n'); "
            "cat(as.character(dup_ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "1,1"
    np.testing.assert_allclose(_parse_float_csv(lines[1]), np.array([2.0, -math.log10(0.2)]))
    assert lines[2] == "2,1"
    assert lines[3] == "TRUE"
