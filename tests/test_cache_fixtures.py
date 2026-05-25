import json
from pathlib import Path

import numpy as np
import pytest

from statgen.annotations import load_annotation, load_annotations, load_annotations_cache
from statgen.genotype import load_genotype, load_genotype_cache
from statgen.reference import load_reference, load_reference_cache
from statgen.sumstats import load_sumstats, load_sumstats_cache
from tests.conftest import (
    FIXTURES_DIR,
    matlab_data_lines,
    run_octave,
    run_rscript,
    skipif_no_octave,
    skipif_no_rscript,
)

CACHE_DIR = FIXTURES_DIR / "cache"
SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
SUMSTATS = FIXTURES_DIR / "sumstats/traits_complete.tsv.gz"
ANNOTATIONS = [FIXTURES_DIR / "annotations/anno1.bed", FIXTURES_DIR / "annotations/anno2.bed"]
CONTINUOUS_ANNOT = FIXTURES_DIR / "annotations/continuous.annot"
CONTINUOUS_META = FIXTURES_DIR / "annotations/continuous.meta"
ANNOTATIONS_REL = [Path("tests/fixtures/annotations/anno1.bed"), Path("tests/fixtures/annotations/anno2.bed")]
CONTINUOUS_ANNOT_REL = Path("tests/fixtures/annotations/continuous.annot")
CONTINUOUS_META_REL = Path("tests/fixtures/annotations/continuous.meta")
GENOTYPE = "tests/fixtures/genotype/sharded/@"

EXPECTED_CACHE_FIXTURES = [
    ("reference", "python", "reference_cache_0_1", "npz", "py"),
    ("reference", "matlab", "reference_cache_0_1", "mat", "m"),
    ("reference", "r", "reference_cache_0_1", "rds", "R"),
    ("sumstats", "python", "sumstats_cache_0_1", "npz", "py"),
    ("sumstats", "matlab", "sumstats_cache_0_1", "mat", "m"),
    ("sumstats", "r", "sumstats_cache_0_1", "rds", "R"),
    ("annotations", "python", "annotations_cache_0_2", "npz", "py"),
    ("annotations", "matlab", "annotations_cache_0_2", "mat", "m"),
    ("annotations", "r", "annotations_cache_0_2", "rds", "R"),
    ("annotations", "python", "annotations_cache_0_1", "npz", "py"),
    ("annotations", "matlab", "annotations_cache_0_1", "mat", "m"),
    ("annotations", "r", "annotations_cache_0_1", "rds", "R"),
    ("genotype", "python", "genotype_cache_0_1", "npz", "py"),
    ("genotype", "matlab", "genotype_cache_0_1", "mat", "m"),
    ("genotype", "r", "genotype_cache_0_1", "rds", "R"),
]


def test_cache_fixture_artifacts_have_provenance_generators():
    for object_name, runtime, schema, cache_ext, generator_ext in EXPECTED_CACHE_FIXTURES:
        stem = f"{object_name}_{runtime}_{schema}"
        assert (CACHE_DIR / f"{stem}.{cache_ext}").is_file()
        assert (CACHE_DIR / f"generate_{stem}.{generator_ext}").is_file()


def _source_reference():
    return load_reference(SHARDED_REF)


def _source_annotations(reference):
    binary = load_annotations(ANNOTATIONS_REL, reference)
    continuous = load_annotation(
        CONTINUOUS_ANNOT_REL,
        reference,
        header=True,
        value_columns=["score", "weight"],
        annotation_metadata_path=CONTINUOUS_META_REL,
    )
    return binary.union_annotations(continuous)


def test_python_cache_fixtures_load_and_match_sources():
    ref = _source_reference()
    ref_cached = load_reference_cache(CACHE_DIR / "reference_python_reference_cache_0_1.npz")
    assert [s.label for s in ref_cached.shards] == [s.label for s in ref.shards]
    np.testing.assert_array_equal(ref_cached.bp, ref.bp)
    np.testing.assert_array_equal(ref_cached.snp, ref.snp)
    np.testing.assert_array_equal(ref_cached.a1, ref.a1)
    np.testing.assert_array_equal(ref_cached.a2, ref.a2)
    np.testing.assert_array_equal(ref_cached.a1_hash64, ref.a1_hash64)
    np.testing.assert_array_equal(ref_cached.a2_hash64, ref.a2_hash64)

    sumstats = load_sumstats(SUMSTATS, ref)
    sumstats_cached = load_sumstats_cache(CACHE_DIR / "sumstats_python_sumstats_cache_0_1.npz")
    np.testing.assert_allclose(sumstats_cached.logpvec, sumstats.logpvec, equal_nan=True)
    np.testing.assert_allclose(sumstats_cached.zvec, sumstats.zvec, equal_nan=True)
    np.testing.assert_allclose(sumstats_cached.nvec, sumstats.nvec, equal_nan=True)
    assert sumstats_cached.beta_vec is None
    assert sumstats_cached.se_vec is None
    assert sumstats_cached.eaf_vec is None
    assert sumstats_cached.info_vec is None

    annotations = _source_annotations(ref)
    annotations_cached = load_annotations_cache(CACHE_DIR / "annotations_python_annotations_cache_0_2.npz")
    assert list(annotations_cached.annonames) == list(annotations.annonames)
    np.testing.assert_array_equal(annotations_cached.is_binary, annotations.is_binary)
    assert list(annotations_cached.annotation_metadata) == list(annotations.annotation_metadata)
    np.testing.assert_allclose(annotations_cached.annomat.toarray(), annotations.annomat.toarray())

    binary_annotations = load_annotations(ANNOTATIONS_REL, ref)
    annotations_old = load_annotations_cache(CACHE_DIR / "annotations_python_annotations_cache_0_1.npz")
    assert list(annotations_old.annonames) == list(binary_annotations.annonames)
    np.testing.assert_array_equal(annotations_old.is_binary, np.array([True, True]))
    assert list(annotations_old.annotation_metadata) == ["", ""]
    np.testing.assert_allclose(annotations_old.annomat.toarray(), binary_annotations.annomat.toarray())

    genotype = load_genotype(GENOTYPE, ref)
    genotype_cached = load_genotype_cache(CACHE_DIR / "genotype_python_genotype_cache_0_1.npz")
    np.testing.assert_array_equal(genotype_cached.is_present, genotype.is_present)
    np.testing.assert_allclose(genotype_cached.ploidy_male, genotype.ploidy_male, equal_nan=True)
    np.testing.assert_allclose(genotype_cached.ploidy_female, genotype.ploidy_female, equal_nan=True)
    np.testing.assert_array_equal(genotype_cached.source_row0, genotype.source_row0)
    np.testing.assert_array_equal(genotype_cached.fid, genotype.fid)
    np.testing.assert_array_equal(genotype_cached.iid, genotype.iid)
    np.testing.assert_array_equal(
        genotype_cached.fetch_genotypes_int8([0, 5, 7]),
        genotype.fetch_genotypes_int8([0, 5, 7]),
    )


def _r_script(expr: str) -> str:
    cache_helper = CACHE_DIR / "cache_fixture_helpers.R"
    return (
        f"statgen_cache_dir <- {json.dumps(str(CACHE_DIR))}; "
        f"source({json.dumps(str(cache_helper))}); "
        "repo_root <- statgen_repo_root(); "
        "statgen_source_r_package(repo_root); "
        + expr
    )


@pytest.mark.r
@skipif_no_rscript
def test_r_cache_fixtures_load_and_match_sources():
    result = run_rscript(
        _r_script(
            "ref <- statgen_cache_fixture_reference(repo_root); "
            f"ref_cached <- load_reference_cache({json.dumps(str(CACHE_DIR / 'reference_r_reference_cache_0_1.rds'))}); "
            "stopifnot(identical(bp(ref_cached), bp(ref))); "
            "stopifnot(identical(snp(ref_cached), snp(ref))); "
            "stopifnot(identical(a1(ref_cached), a1(ref))); "
            "stopifnot(identical(a2(ref_cached), a2(ref))); "
            f"s <- load_sumstats({json.dumps(str(SUMSTATS))}, ref); "
            f"s_cached <- load_sumstats_cache({json.dumps(str(CACHE_DIR / 'sumstats_r_sumstats_cache_0_1.rds'))}); "
            "stopifnot(isTRUE(all.equal(logpvec(s_cached), logpvec(s), check.attributes = FALSE))); "
            "stopifnot(isTRUE(all.equal(zvec(s_cached), zvec(s), check.attributes = FALSE))); "
            "stopifnot(isTRUE(all.equal(nvec(s_cached), nvec(s), check.attributes = FALSE))); "
            "a <- statgen_cache_fixture_annotations(repo_root, ref); "
            f"a_cached <- load_annotations_cache({json.dumps(str(CACHE_DIR / 'annotations_r_annotations_cache_0_2.rds'))}); "
            "stopifnot(identical(annonames(a_cached), annonames(a))); "
            "stopifnot(identical(is_binary(a_cached), is_binary(a))); "
            "stopifnot(identical(annotation_metadata(a_cached), annotation_metadata(a))); "
            "stopifnot(isTRUE(all.equal(as.matrix(annomat(a_cached)), as.matrix(annomat(a)), check.attributes = FALSE))); "
            f"a_old <- load_annotations_cache({json.dumps(str(CACHE_DIR / 'annotations_r_annotations_cache_0_1.rds'))}); "
            f"a_binary <- load_annotations(c({json.dumps(str(ANNOTATIONS_REL[0]))}, {json.dumps(str(ANNOTATIONS_REL[1]))}), ref); "
            "stopifnot(identical(annonames(a_old), annonames(a_binary))); "
            "stopifnot(identical(is_binary(a_old), c(TRUE, TRUE))); "
            "stopifnot(identical(annotation_metadata(a_old), c('', ''))); "
            "stopifnot(isTRUE(all.equal(as.matrix(annomat(a_old)), as.matrix(annomat(a_binary)), check.attributes = FALSE))); "
            f"g <- load_genotype({json.dumps(GENOTYPE)}, ref); "
            f"g_cached <- load_genotype_cache({json.dumps(str(CACHE_DIR / 'genotype_r_genotype_cache_0_1.rds'))}); "
            "stopifnot(identical(is_present(g_cached), is_present(g))); "
            "stopifnot(isTRUE(all.equal(ploidy_male(g_cached), ploidy_male(g), check.attributes = FALSE))); "
            "stopifnot(isTRUE(all.equal(ploidy_female(g_cached), ploidy_female(g), check.attributes = FALSE))); "
            "stopifnot(identical(source_row0(g_cached), source_row0(g))); "
            "stopifnot(identical(fid(g_cached), fid(g))); "
            "stopifnot(identical(iid(g_cached), iid(g))); "
            "stopifnot(identical(fetch_genotypes_int8(g_cached, c(1L, 6L, 8L)), fetch_genotypes_int8(g, c(1L, 6L, 8L)))); "
            "cat('OK\\n')"
        ),
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "OK"


@pytest.mark.octave
@skipif_no_octave
def test_matlab_cache_fixtures_load_and_match_sources():
    script = (
        "fixture_dir = 'tests/fixtures'; "
        "cache_dir = 'tests/fixtures/cache'; "
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "ref_cached = statgen.load_reference_cache([cache_dir '/reference_matlab_reference_cache_0_1.mat']); "
        "ok = isequal(ref_cached.bp, ref.bp) && isequal(ref_cached.snp, ref.snp) && isequal(ref_cached.a1, ref.a1) && isequal(ref_cached.a2, ref.a2); "
        "s = statgen.load_sumstats([fixture_dir '/sumstats/traits_complete.tsv.gz'], ref); "
        "s_cached = statgen.load_sumstats_cache([cache_dir '/sumstats_matlab_sumstats_cache_0_1.mat']); "
        "ok = ok && isequaln(s_cached.logpvec, s.logpvec) && isequaln(s_cached.zvec, s.zvec) && isequaln(s_cached.nvec, s.nvec); "
        "binary = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        "continuous = statgen.load_annotation([fixture_dir '/annotations/continuous.annot'], ref, 'header', true, 'value_columns', {'score', 'weight'}, 'annotation_metadata_path', [fixture_dir '/annotations/continuous.meta']); "
        "a = binary.union_annotations(continuous); "
        "a_cached = statgen.load_annotations_cache([cache_dir '/annotations_matlab_annotations_cache_0_2.mat']); "
        "ok = ok && isequal(a_cached.annonames, a.annonames) && isequal(a_cached.is_binary, a.is_binary) && isequal(a_cached.annotation_metadata, a.annotation_metadata); "
        "ok = ok && isequal(full(a_cached.annomat), full(a.annomat)); "
        "a_old = statgen.load_annotations_cache([cache_dir '/annotations_matlab_annotations_cache_0_1.mat']); "
        "ok = ok && isequal(a_old.annonames, binary.annonames) && isequal(a_old.is_binary, [true; true]) && isequal(a_old.annotation_metadata, {'', ''}'); "
        "ok = ok && isequal(full(a_old.annomat), full(binary.annomat)); "
        "g = statgen.load_genotype('tests/fixtures/genotype/sharded/@', ref); "
        "g_cached = statgen.load_genotype_cache([cache_dir '/genotype_matlab_genotype_cache_0_1.mat']); "
        "ok = ok && isequal(g_cached.is_present, g.is_present) && isequaln(g_cached.ploidy_male, g.ploidy_male) && isequaln(g_cached.ploidy_female, g.ploidy_female); "
        "ok = ok && isequal(g_cached.source_row0, g.source_row0) && isequal(g_cached.fid, g.fid) && isequal(g_cached.iid, g.iid); "
        "ok = ok && isequal(g_cached.fetch_genotypes_int8([1 6 8]), g.fetch_genotypes_int8([1 6 8])); "
        "fprintf('%d\\n', ok);"
    )
    result = run_octave(script, timeout=60)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1"]
