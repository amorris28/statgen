"""R phase 4/5 acceptance tests: direct LD NPZ loading and R sidecar preparation."""

import json
import shutil

import numpy as np
import pytest

from tests.conftest import FIXTURES_DIR, R_PACKAGE_DIR, run_rscript, skipif_no_rscript
from tests.test_end_to_end_alignment import (
    _assert_full_outputs,
    _assert_preflight_outputs,
    _build_full_artifacts,
    _build_preflight_artifacts,
    _full_sample_ids,
    _full_x_subset,
    _python_preflight_outputs,
    _read_labeled_numeric_output,
)


LD_PY = FIXTURES_DIR / "ld/python"


def _valid_npz_payload(metadata=None):
    if metadata is None:
        metadata = {
            "object_type": "ld_shard",
            "schema_version": "1.0",
            "format": "statgen_ld_npz_csc32",
            "chr": "1",
            "sex": None,
            "num_snp": 1,
            "nnz": 1,
            "matrix": "symmetric",
            "diagonal": "explicit_unit",
            "value": "r",
            "reference_checksum": "checksum",
            "reference_bim": "reference_chr1.bim",
            "num_monomorphic_snps": 0,
            "sparse_layout": "csc",
            "index_base": 0,
        }
    return {
        "data": np.array([1.0], dtype=np.float32),
        "indices": np.array([0], dtype=np.int32),
        "indptr": np.array([0, 1], dtype=np.int32),
        "shape": np.array([1, 1], dtype=np.int64),
        "a1freq": np.array([0.2], dtype=np.float32),
        "metadata": np.frombuffer(json.dumps(metadata).encode("utf-8"), dtype=np.uint8),
    }


def _source_phase5_script(expr: str) -> str:
    files = [
        "utils.R",
        "verbosity.R",
        "hash.R",
        "variant_match.R",
        "bfile_utils.R",
        "reference.R",
        "ld_schema.R",
        "ld_reference.R",
        "ld_csc.R",
        "ld_npz.R",
        "ld.R",
    ]
    sources = " ".join(
        f"source({json.dumps(str(R_PACKAGE_DIR / 'R' / filename))});"
        for filename in files
    )
    return f"{sources} {expr}"


def _source_phase5_genotype_script(expr: str) -> str:
    files = [
        "utils.R",
        "verbosity.R",
        "hash.R",
        "variant_match.R",
        "bfile_utils.R",
        "reference.R",
        "genotype.R",
        "ld_schema.R",
        "ld_reference.R",
        "ld_csc.R",
        "ld_npz.R",
        "ld.R",
    ]
    sources = " ".join(
        f"source({json.dumps(str(R_PACKAGE_DIR / 'R' / filename))});"
        for filename in files
    )
    return f"{sources} {expr}"


def _r_quote(path) -> str:
    return json.dumps(str(path))


def _copy_ld_py(tmp_path):
    root = tmp_path / "ld_npz"
    shutil.copytree(LD_PY, root)
    return root


def _rewrite_npz_metadata(path, mutate):
    with np.load(path) as loaded:
        payload = {name: loaded[name] for name in loaded.files}
    metadata = json.loads(payload["metadata"].tobytes().decode("utf-8"))
    mutate(metadata)
    payload["metadata"] = np.frombuffer(json.dumps(metadata).encode("utf-8"), dtype=np.uint8)
    np.savez(path, **payload)


@pytest.mark.r
@skipif_no_rscript
def test_r_ld_npz_preparation_loads_and_operates(tmp_path):
    root = _copy_ld_py(tmp_path)
    result = run_rscript(
        _source_phase5_script(
            "set_verbosity('quiet'); "
            f"manifest <- prepare_ld_npz_for_r({json.dumps(str(root))}); "
            f"report <- validate_ld_distribution({json.dumps(str(root))}, TRUE); "
            f"report_default <- validate_ld_distribution({json.dumps(str(root))}); "
            f"ld <- load_ld({json.dumps(str(root))}); "
            "off <- shard_offsets(ld); "
            "cat(as.character(report$ok), as.character(report_default$ok), default_chrX_sex(ld), num_snp(ld), '\\n'); "
            "cat(paste(vapply(shards(ld), function(s) s$label, character(1)), collapse=','), '\\n'); "
            "cat(typeof(off$start0), typeof(off$stop0), '\\n'); "
            "cat(paste(round(a1freq(ld), 6), collapse=','), '\\n'); "
            "cat(paste(round(a1freq(ld, chrX_sex='combined'), 6), collapse=','), '\\n'); "
            "cat(paste(round(multiply_r2(ld, seq_len(num_snp(ld))), 6), collapse=','), '\\n'); "
            "cat(paste(round(multiply_r2(ld, seq_len(num_snp(ld)), chrX_sex='male'), 6), collapse=','), '\\n'); "
            "cat(paste(round(multiply_r2(ld, seq_len(num_snp(ld)), chrX_sex='combined'), 6), collapse=','), '\\n'); "
            "mat_out <- multiply_r2(ld, cbind(seq_len(num_snp(ld)), 2 * seq_len(num_snp(ld)))); "
            "cat(paste(dim(mat_out), collapse='x'), paste(round(mat_out[, 1], 6), collapse=','), paste(round(mat_out[, 2], 6), collapse=','), '\\n'); "
            "cat(paste(round(fast_prune(c(9,9,7,1,2,6,5,4), ld, 0.35), 6), collapse=','), '\\n'); "
            "cat(paste(round(fast_prune(c(9,9,7,1,2,6,5,4), ld, 0.35, chrX_sex='male'), 6), collapse=','), '\\n'); "
            "cat(paste(round(fast_prune(c(9,9,7,1,2,6,5,4), ld, 0.35, chrX_sex='combined'), 6), collapse=','), '\\n'); "
            "ld_x <- select_shards(ld, 'X'); "
            "cat(num_snp(ld_x), paste(vapply(shards(ld_x), function(s) s$label, character(1)), collapse=','), num_snp(reference(ld_x)), paste(round(a1freq(ld_x, chrX_sex='male'), 6), collapse=','), '\\n'); "
            "sparse_out <- multiply_r2(ld, Matrix::sparseMatrix(i=1:num_snp(ld), j=rep(1L, num_snp(ld)), x=seq_len(num_snp(ld)), dims=c(num_snp(ld), 1L))); "
            "dense_out <- multiply_r2(ld, Matrix::Matrix(matrix(seq_len(num_snp(ld)), ncol=1L), sparse=FALSE)); "
            "cat(as.character(inherits(sparse_out, 'sparseMatrix')), as.character(is.matrix(dense_out)), '\\n'); "
            f"ld2 <- load_ld({json.dumps(str(root))}, retain_ld_r = FALSE); "
            "cat(as.character(is.null(ld2$shard_groups[[1]][[1]]$ld_r)), round(multiply_r2(ld2, seq_len(num_snp(ld2)))[1], 6), "
            "paste(round(fast_prune(c(9,9,7,1,2,6,5,4), ld2, 0.35, chrX_sex='combined'), 6), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "TRUE TRUE female 8"
    assert lines[1] == "1,X"
    assert lines[2] == "integer integer"
    np.testing.assert_allclose(
        [float(x) for x in lines[3].split(",")],
        [0.30, 0.40, 0.20, 0.35, 0.15, 0.28, 0.32, 0.22],
        atol=1e-6,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[4].split(",")],
        [0.30, 0.40, 0.20, 0.35, 0.15, 0.25, 0.30, 0.20],
        atol=1e-6,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[5].split(",")],
        [2.89, 3.81, 3.09, 6.95, 6.96, 8.9575, 10.515, 8.8575],
        atol=1e-6,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[6].split(",")],
        [2.89, 3.81, 3.09, 6.95, 6.96, 8.1175, 10.435, 9.4175],
        atol=1e-6,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[7].split(",")],
        [2.89, 3.81, 3.09, 6.95, 6.96, 8.52, 10.44, 9.12],
        atol=1e-6,
    )
    mat_parts = lines[8].split()
    assert mat_parts[0] == "8x2"
    np.testing.assert_allclose(
        [float(x) for x in mat_parts[1].split(",")],
        [2.89, 3.81, 3.09, 6.95, 6.96, 8.9575, 10.515, 8.8575],
        atol=1e-6,
    )
    np.testing.assert_allclose(
        [float(x) for x in mat_parts[2].split(",")],
        [5.78, 7.62, 6.18, 13.9, 13.92, 17.915, 21.03, 17.715],
        atol=1e-6,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[9].split(",")],
        [9, np.nan, 7, np.nan, 2, 6, np.nan, 4],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[10].split(",")],
        [9, np.nan, 7, np.nan, 2, 6, 5, 4],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        [float(x) for x in lines[11].split(",")],
        [9, np.nan, 7, np.nan, 2, 6, np.nan, 4],
        equal_nan=True,
    )
    assert lines[12] == "3 X 3 0.22,0.28,0.18"
    assert lines[13] == "TRUE TRUE"
    retain_parts = lines[14].split(maxsplit=2)
    assert retain_parts[:2] == ["TRUE", "2.89"]
    np.testing.assert_allclose(
        [float(x) for x in retain_parts[2].split(",")],
        [9, np.nan, 7, np.nan, 2, 6, np.nan, 4],
        equal_nan=True,
    )


@pytest.mark.r
@skipif_no_rscript
def test_r_preflight_chr1_a1freq_and_ld_r_match_fetched_genotypes(tmp_path):
    prefix, ld_root = _build_preflight_artifacts(tmp_path)
    expected = _python_preflight_outputs(prefix, ld_root)
    out_path = tmp_path / "r_preflight.tsv"

    result = run_rscript(
        _source_phase5_genotype_script(
            "set_verbosity('quiet'); "
            f"manifest <- prepare_ld_npz_for_r({_r_quote(ld_root)}); "
            f"report <- validate_ld_distribution({_r_quote(ld_root)}, TRUE); "
            f"ld <- load_ld({_r_quote(ld_root)}); "
            f"g <- load_genotype({_r_quote(prefix)}, reference(ld)); "
            "G <- fetch_genotypes(g, seq_len(num_snp(g))); "
            "obs <- is.finite(G); G0 <- G; G0[!obs] <- 0; "
            "freq <- colSums(G0) / (2 * colSums(obs)); "
            "ld_r <- ld$shard_groups[[1]][[1]]$ld_r; "
            "vec <- seq(-1.5, 2.5, length.out = num_snp(g)); "
            "mat <- cbind(vec, rev(vec)); "
            "mvec <- multiply_r2(ld, vec); "
            "mmat <- multiply_r2(ld, mat); "
            "explicit_vec <- as.numeric((ld_r ^ 2) %*% vec); "
            "explicit_mat <- as.numeric((ld_r ^ 2) %*% mat); "
            "r01_keep <- is.finite(G[, 1]) & is.finite(G[, 2]); "
            "r02_keep <- is.finite(G[, 1]) & is.finite(G[, 3]); "
            "direct <- c(cor(G[r01_keep, 1], G[r01_keep, 2]), cor(G[r02_keep, 1], G[r02_keep, 3])); "
            f"out_file <- {_r_quote(out_path)}; "
            "write_line <- function(name, values) { "
            "cat(name, paste(format(as.numeric(values), digits = 12, scientific = FALSE), collapse = ' '), '\\n', file = out_file, append = TRUE) }; "
            "write_line('manifest_count', length(manifest$shards)); "
            "write_line('validate_ok', as.integer(report$ok)); "
            "write_line('is_present', as.integer(is_present(g))); "
            "write_line('ld_a1freq', a1freq(ld)); "
            "write_line('fetch_a1freq', freq); "
            "write_line('ld_r_selected', c(as.numeric(ld_r[1, 2]), as.numeric(ld_r[1, 3]))); "
            "write_line('direct_r_selected', direct); "
            "write_line('multiply_vec', mvec); "
            "write_line('explicit_vec', explicit_vec); "
            "write_line('multiply_mat', as.numeric(mmat)); "
            "write_line('explicit_mat', explicit_mat); "
            "invisible(NULL)"
        ),
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    actual = _read_labeled_numeric_output(out_path)
    assert actual["manifest_count"].tolist() == [1.0]
    assert actual["validate_ok"].tolist() == [1.0]
    for key in [
        "is_present",
        "ld_a1freq",
        "fetch_a1freq",
        "ld_r_selected",
        "direct_r_selected",
        "multiply_vec",
        "explicit_vec",
        "multiply_mat",
        "explicit_mat",
    ]:
        np.testing.assert_allclose(actual[key], expected[key], atol=1e-5)
    _assert_preflight_outputs(actual)


@pytest.mark.r
@skipif_no_rscript
def test_r_full_e2e_partial_overlap_a1freq_ld_r_and_chrx_subject_mapping(tmp_path):
    artifacts = _build_full_artifacts(tmp_path)

    def run_case(ld_root, genotype_prefix, out_path):
        script_path = out_path.with_suffix(".R")
        script_code = _source_phase5_genotype_script(
                "set_verbosity('quiet'); "
                f"manifest <- prepare_ld_npz_for_r({_r_quote(ld_root)}); "
                f"report <- validate_ld_distribution({_r_quote(ld_root)}, TRUE); "
                f"ld <- load_ld({_r_quote(ld_root)}); "
                f"g <- load_genotype({_r_quote(genotype_prefix)}, reference(ld)); "
                "present_idx <- which(is_present(g)); "
                "G <- fetch_genotypes(g, present_idx); "
                "freq_f <- rep(NaN, num_snp(g)); freq_m <- rep(NaN, num_snp(g)); "
                "goff <- shard_offsets(g); x_row <- which(goff$shard_label == 'X'); x_start <- goff$start0[[x_row]] + 1L; "
                "mask_x_f <- is_subject_present(g, 'X') & is_female(g); "
                "mask_x_m <- is_subject_present(g, 'X') & is_male(g); "
                "for (col in seq_along(present_idx)) { "
                "gi <- present_idx[[col]]; "
                "if (gi >= x_start) { "
                "vf <- G[mask_x_f, col]; of <- is.finite(vf); if (any(of)) { freq_f[[gi]] <- sum(vf[of]) / (2 * sum(of)); }; "
                "vm <- G[mask_x_m, col]; om <- is.finite(vm); if (any(om)) { freq_m[[gi]] <- sum(pmin(vm[om], 1)) / sum(om); }; "
                "} else { "
                "v <- G[, col]; o <- is.finite(v); if (any(o)) { f <- sum(v[o]) / (2 * sum(o)); freq_f[[gi]] <- f; freq_m[[gi]] <- f; }; "
                "} "
                "}; "
                "find_idx <- function(marker) which(snp(reference(ld)) == marker)[[1]]; "
                "pair_for <- function(chrom) { ids <- grep(paste0('^e2e_full_chr', chrom, '_.*_shared_'), snp(reference(ld)), value = TRUE); c(find_idx(ids[[1]]), find_idx(ids[[2]])) }; "
                "pairs <- list(pair_for('1'), pair_for('2'), pair_for('X')); "
                "bounds <- shard_offsets(reference(ld)); "
                "shard_row <- function(idx) which(idx > bounds$start0 & idx <= bounds$stop0)[[1]]; "
                "ld_value <- function(pair, sex_label) { row <- shard_row(pair[[1]]); local <- pair - bounds$start0[[row]]; group <- ld$shard_groups[[row]]; "
                "shard <- if (identical(bounds$shard_label[[row]], 'X')) group[[match(sex_label, vapply(group, function(s) s$sex, character(1)))]] else group[[1]]; "
                "as.numeric(shard$ld_r[local[[1]], local[[2]]]) }; "
                "direct_pair <- function(pair, sex_label) { P <- fetch_genotypes(g, pair); row <- shard_row(pair[[1]]); label <- bounds$shard_label[[row]]; "
                "if (identical(label, 'X')) { mask <- is_subject_present(g, 'X'); if (identical(sex_label, 'female')) { mask <- mask & is_female(g) } else { mask <- mask & is_male(g) } } else { mask <- rep(TRUE, num_sample(g)) }; "
                "x <- P[mask, 1]; y <- P[mask, 2]; keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]; "
                "if (identical(label, 'X') && identical(sex_label, 'male')) { x <- pmin(x, 1); y <- pmin(y, 1) }; cor(x, y) }; "
                "ld_r_selected <- c(ld_value(pairs[[1]], 'female'), ld_value(pairs[[2]], 'female'), ld_value(pairs[[3]], 'female'), ld_value(pairs[[3]], 'male')); "
                "direct_r_selected <- c(direct_pair(pairs[[1]], 'female'), direct_pair(pairs[[2]], 'female'), direct_pair(pairs[[3]], 'female'), direct_pair(pairs[[3]], 'male')); "
                "vec <- seq(-2.0, 3.0, length.out = num_snp(ld)); "
                "mf <- multiply_r2(ld, vec, chrX_sex = 'female'); mm <- multiply_r2(ld, vec, chrX_sex = 'male'); "
                "af <- a1freq(ld, chrX_sex = 'female'); am <- a1freq(ld, chrX_sex = 'male'); "
                f"out_file <- {_r_quote(out_path)}; "
                "write_line <- function(name, values) { "
                "cat(name, paste(format(as.numeric(values), digits = 12, scientific = FALSE), collapse = ' '), '\\n', file = out_file, append = TRUE) }; "
                "write_line('manifest_count', length(manifest$shards)); "
                "write_line('validate_ok', as.integer(report$ok)); "
                "write_line('present', as.integer(is_present(g))); "
                "write_line('chrx_subject_present', as.integer(is_subject_present(g, 'X'))); "
                "write_line('ld_a1freq_female_present', af[present_idx]); "
                "write_line('calc_a1freq_female_present', freq_f[present_idx]); "
                "write_line('ld_a1freq_male_present', am[present_idx]); "
                "write_line('calc_a1freq_male_present', freq_m[present_idx]); "
                "write_line('ld_r_selected', ld_r_selected); "
                "write_line('direct_r_selected', direct_r_selected); "
                "write_line('multiply_female_head', mf[1:10]); "
                "write_line('multiply_male_head', mm[1:10]); "
                "invisible(NULL)"
        )
        script_path.write_text(script_code.replace("; ", ";\n"), encoding="utf-8")
        result = run_rscript(
            f"source({_r_quote(script_path)})",
            timeout=180,
        )
        assert result.returncode == 0, result.stderr
        actual = _read_labeled_numeric_output(out_path)
        assert actual["manifest_count"].tolist() == [4.0]
        assert actual["validate_ok"].tolist() == [1.0]
        _assert_full_outputs(actual)
        return actual

    nonsharded = run_case(
        artifacts["ld_full_root"],
        artifacts["genotype_prefix"],
        tmp_path / "r_full_nonsharded.tsv",
    )
    sharded = run_case(
        artifacts["ld_sharded_root"],
        artifacts["genotype_sharded"],
        tmp_path / "r_full_sharded.tsv",
    )

    assert nonsharded["chrx_subject_present"].sum() == 40
    expected_mask = np.zeros_like(sharded["chrx_subject_present"])
    x_samples = set(_full_x_subset(_full_sample_ids()))
    genotype_iids = [f"S{i + 1:02d}" for i in range(40)]
    for i, iid in enumerate(genotype_iids):
        expected_mask[i] = 1.0 if iid in x_samples else 0.0
    np.testing.assert_array_equal(sharded["chrx_subject_present"], expected_mask)


@pytest.mark.r
@skipif_no_rscript
def test_r_prepare_ld_npz_for_r_writes_reference_cache_and_md5(tmp_path):
    root = _copy_ld_py(tmp_path)
    result = run_rscript(
        _source_phase5_script(
            f"manifest <- prepare_ld_npz_for_r({json.dumps(str(root))}, extract_npz = TRUE); "
            f"manifest2 <- prepare_ld_npz_for_r({json.dumps(str(root))}); "
            "cat(manifest$runtime_format, manifest$r_reference_cache, grepl('^[0-9a-f]{32}$', manifest$r_reference_cache_md5), '\\n'); "
            "cat(file.exists(file.path("
            f"{json.dumps(str(root))}, manifest$r_reference_cache)), "
            "file.exists(file.path("
            f"{json.dumps(str(root))}, manifest$reference_cache)), "
            "identical(manifest$r_reference_cache, manifest2$r_reference_cache), "
            "grepl('^[0-9a-f]{32}$', manifest2$r_reference_cache_md5), '\\n'); "
            f"scratch_blocker <- file.path({json.dumps(str(root))}, 'scratch-blocker'); "
            "invisible(file.create(scratch_blocker)); "
            "old_scratch <- Sys.getenv('STATGEN_SCRATCH', unset = NA_character_); "
            "Sys.setenv(STATGEN_SCRATCH = scratch_blocker); "
            "on.exit(if (is.na(old_scratch)) Sys.unsetenv('STATGEN_SCRATCH') else Sys.setenv(STATGEN_SCRATCH = old_scratch), add = TRUE); "
            "cache_dir <- file.path("
            f"{json.dumps(str(root))}, paste0(manifest$shards[[1]]$file, '.d')); "
            "marker <- jsonlite::fromJSON(file.path(cache_dir, '.statgen-extracted-npz.json')); "
            "cat(dir.exists(cache_dir), file.exists(file.path(cache_dir, 'data.npy')), marker$source_file, marker$source_file_md5 == manifest$shards[[1]]$file_md5, '\\n'); "
            f"ref <- load_ld_reference({json.dumps(str(root))}); "
            f"ld <- load_ld({json.dumps(str(root))}); "
            "cat(num_snp(ref), paste(vapply(shards(ref), function(s) s$label, character(1)), collapse=','), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    lines = [line.strip() for line in result.stdout.strip().splitlines()]
    assert lines[0] == "python_npz_csc32 reference_cache.rds TRUE"
    assert lines[1] == "TRUE TRUE TRUE TRUE"
    assert lines[2] == "TRUE TRUE ld_chr1.npz TRUE"
    assert lines[3] == "8 1,X"


@pytest.mark.r
@skipif_no_rscript
def test_r_prepare_ld_npz_for_r_rejects_invalid_existing_sidecar_filename(tmp_path):
    root = _copy_ld_py(tmp_path)
    manifest_path = root / "ld_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["r_reference_cache"] = "subdir/reference_cache.rds"
    manifest["r_reference_cache_md5"] = "0" * 32
    manifest_path.write_text(json.dumps(manifest) + "\n")

    result = run_rscript(
        _source_phase5_script(
            "ok <- FALSE; "
            f"tryCatch(prepare_ld_npz_for_r({json.dumps(str(root))}), "
            "error = function(e) ok <<- grepl('plain relative filename', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_load_ld_emits_monomorphic_snp_warning(tmp_path):
    root = _copy_ld_py(tmp_path)
    _rewrite_npz_metadata(root / "ld_chr1.npz", lambda metadata: metadata.update(num_monomorphic_snps=1))
    result = run_rscript(
        _source_phase5_script(
            f"prepare_ld_npz_for_r({json.dumps(str(root))}); "
            "warned <- FALSE; "
            "ld <- withCallingHandlers("
            f"load_ld({json.dumps(str(root))}), "
            "warning = function(w) { warned <<- grepl('monomorphic SNPs', w$message); invokeRestart('muffleWarning') }); "
            "cat(as.character(warned), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_npz_reader_rejects_malformed_archives(tmp_path):
    missing_arrays = tmp_path / "missing_arrays.npz"
    wrong_dtype = tmp_path / "wrong_dtype.npz"
    fortran_order = tmp_path / "fortran_order.npz"
    malformed = tmp_path / "malformed.npz"

    np.savez(missing_arrays, data=np.array([1.0], dtype=np.float32))

    payload = _valid_npz_payload()
    payload["data"] = np.array([1.0], dtype=np.float64)
    np.savez(wrong_dtype, **payload)

    payload = _valid_npz_payload()
    payload["data"] = np.asfortranarray(np.eye(2, dtype=np.float32))
    np.savez(fortran_order, **payload)

    malformed.write_bytes(b"not a zip archive")

    result = run_rscript(
        _source_phase5_script(
            "ok_missing <- FALSE; "
            "ok_dtype <- FALSE; "
            "ok_fortran <- FALSE; "
            "ok_malformed <- FALSE; "
            f"tryCatch(.read_npz_ld_payload({json.dumps(str(missing_arrays))}, scratch_parent = {json.dumps(str(tmp_path))}), "
            "error = function(e) ok_missing <<- grepl('missing required arrays', e$message)); "
            f"tryCatch(.read_npz_ld_payload({json.dumps(str(wrong_dtype))}, scratch_parent = {json.dumps(str(tmp_path))}), "
            "error = function(e) ok_dtype <<- grepl('must have dtype', e$message)); "
            f"tryCatch(.read_npz_ld_payload({json.dumps(str(fortran_order))}, scratch_parent = {json.dumps(str(tmp_path))}), "
            "error = function(e) ok_fortran <<- grepl('Fortran-ordered', e$message)); "
            f"tryCatch(.read_npz_ld_payload({json.dumps(str(malformed))}, scratch_parent = {json.dumps(str(tmp_path))}), "
            "error = function(e) ok_malformed <<- TRUE); "
            "cat(as.character(ok_missing), as.character(ok_dtype), as.character(ok_fortran), as.character(ok_malformed), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE TRUE TRUE TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_load_ld_rejects_missing_r_reference_sidecar():
    result = run_rscript(
        _source_phase5_script(
            "ok <- FALSE; "
            f"tryCatch(load_ld({json.dumps(str(LD_PY))}), error = function(e) ok <<- grepl('missing r_reference_cache', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_ld_metadata_rejects_csc32_nnz_upper_bound():
    result = run_rscript(
        _source_phase5_script(
            "meta <- list("
            "object_type = 'ld_shard', schema_version = .ld_manifest_schema, format = .ld_npz_format, "
            "chr = '1', sex = NULL, num_snp = 1, nnz = 2^31, "
            "matrix = 'symmetric', diagonal = 'explicit_unit', value = 'r', "
            "reference_checksum = 'checksum', reference_bim = 'reference_chr1.bim', "
            "num_monomorphic_snps = 0, sparse_layout = 'csc', index_base = 0); "
            "ok <- FALSE; "
            "tryCatch(.validate_ld_shard_metadata(meta, 'shard.npz', .ld_npz_format), "
            "error = function(e) ok <<- grepl('CSC32 metadata nnz must be < 2\\\\^31', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_ld_manifest_rejects_chrx_reference_identity_mismatch(tmp_path):
    manifest = {
        "object_type": "ld_panel_manifest",
        "schema_version": "1.0",
        "runtime_format": "python_npz_csc32",
        "reference_cache": "reference_cache.npz",
        "reference_cache_md5": "0" * 32,
        "r_reference_cache": "reference_cache.rds",
        "r_reference_cache_md5": "0" * 32,
        "shards": [
            {
                "chr": "X",
                "sex": "female",
                "file": "ld_chrX_female.npz",
                "file_md5": "1" * 32,
                "num_snp": 3,
                "nnz": 7,
                "reference_checksum": "a" * 32,
                "reference_bim": "reference_chrX.bim",
            },
            {
                "chr": "X",
                "sex": "male",
                "file": "ld_chrX_male.npz",
                "file_md5": "2" * 32,
                "num_snp": 4,
                "nnz": 7,
                "reference_checksum": "b" * 32,
                "reference_bim": "reference_chrX.bim",
            },
        ],
    }
    (tmp_path / "ld_manifest.json").write_text(json.dumps(manifest) + "\n")

    result = run_rscript(
        _source_phase5_script(
            "ok <- FALSE; "
            f"tryCatch(load_ld({json.dumps(str(tmp_path))}), "
            "error = function(e) ok <<- grepl('reference_bim, num_snp, and reference_checksum', e$message)); "
            "cat(as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "TRUE"


@pytest.mark.r
@skipif_no_rscript
def test_r_npy_int64_reader_rejects_values_above_r_integer_range():
    result = run_rscript(
        _source_phase5_script(
            "tmp <- tempfile(); "
            "con <- file(tmp, 'wb'); "
            "writeBin(as.raw(c(255,255,255,127,0,0,0,0)), con); "
            "close(con); "
            "con <- file(tmp, 'rb'); "
            "value <- .read_npy_int64_as_integer(con, 1, 'archive.npz', 'shape'); "
            "close(con); "
            "tmp2 <- tempfile(); "
            "con <- file(tmp2, 'wb'); "
            "writeBin(as.raw(c(0,0,0,128,0,0,0,0)), con); "
            "close(con); "
            "ok <- FALSE; "
            "con <- file(tmp2, 'rb'); "
            "tryCatch(.read_npy_int64_as_integer(con, 1, 'archive.npz', 'shape'), "
            "error = function(e) ok <<- grepl('exceed', e$message)); "
            "close(con); "
            "cat(value, as.character(ok), '\\n')"
        )
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "2147483647 TRUE"
