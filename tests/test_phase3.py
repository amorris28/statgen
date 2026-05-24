import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from statgen.annotations import (
    AnnotationPanel,
    AnnotationShard,
    create_annotation,
    create_annotations,
    load_annotation,
    load_annotations,
    load_annotations_cache,
    save_annotations_cache,
)
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, run_octave, skipif_no_octave

SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
ANNO1 = FIXTURES_DIR / "annotations/anno1.bed"
ANNO2 = FIXTURES_DIR / "annotations/anno2.bed"


EXPECTED_MASK = np.array(
    [
        [1, 0],
        [1, 0],
        [1, 0],
        [1, 0],
        [1, 0],
        [0, 1],
        [0, 1],
        [0, 1],
    ],
    dtype=np.uint8,
)


def _write_bed(path: Path, rows: list[tuple[str, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for chr_label, start, end in rows:
            f.write(f"{chr_label}\t{start}\t{end}\n")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _load_cache_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def test_load_annotations_fixture_masks_and_names():
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)

    assert [s.label for s in a.shards] == ["1", "X"]
    assert a.num_snp == 8
    assert a.num_annot == 2
    assert list(a.annonames) == ["anno1", "anno2"]
    np.testing.assert_array_equal(a.is_binary, np.array([True, True]))
    assert all(json.loads(x)["source_file"] for x in a.annotation_metadata)
    assert sparse.isspmatrix_csr(a.annomat)
    np.testing.assert_array_equal(a.annomat.toarray(), EXPECTED_MASK)


def test_load_annotations_accepts_single_path_scalar():
    reference = load_reference(SHARDED_REF)
    a = load_annotations(ANNO1, reference)
    assert list(a.annonames) == ["anno1"]
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), EXPECTED_MASK[:, 0])


def test_empty_bed_file_fails(tmp_path):
    reference = load_reference(SHARDED_REF)
    empty = tmp_path / "empty.bed"
    empty.write_text("")
    with pytest.raises(ValueError, match="BED file is empty"):
        load_annotations(empty, reference)


def test_duplicate_bed_basenames_fail(tmp_path):
    reference = load_reference(SHARDED_REF)
    p1 = tmp_path / "set1" / "dup.bed"
    p2 = tmp_path / "set2" / "dup.bed"
    _write_bed(p1, [("1", 99, 100)])
    _write_bed(p2, [("X", 99, 100)])
    with pytest.raises(ValueError, match="duplicate annotation names"):
        load_annotations([p1, p2], reference)


def test_bed_metadata_and_comment_lines_are_ignored(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "with_headers.bed"
    _write_text(
        bed,
        "# track name=CodingExon description=test\n"
        "# browser position chr1:1-1000\n"
        "# comment\n"
        "\n"
        "1\t99\t200\n"
        "1\t299\t400\n",
    )
    a = load_annotations(bed, reference)
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=np.uint8))


def test_bed_comment_or_blank_lines_are_skipped_anywhere(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "comments_anywhere.bed"
    _write_text(
        bed,
        "1\t99\t100\n"
        "# late comment skipped by pandas comment handling\n"
        "\n"
        "1\t399\t400\n",
    )
    a = load_annotations(bed, reference)
    np.testing.assert_array_equal(
        a.annomat.toarray().reshape(-1),
        np.array([1, 0, 0, 1, 0, 0, 0, 0], dtype=np.uint8),
    )


def test_metadata_only_bed_file_fails(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "metadata_only.bed"
    _write_text(
        bed,
        "# track name=foo\n"
        "# browser position chr1:1-100\n"
        "# note\n"
        "\n",
    )
    with pytest.raises(ValueError, match="BED file is empty"):
        load_annotations(bed, reference)


def test_non_tab_delimited_bed_rows_fail(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "space_delimited.bed"
    _write_text(bed, "1  99   200\n1  299  400\n")
    with pytest.raises(ValueError, match="BED must have at least 3 tab-separated columns"):
        load_annotations(bed, reference)


def test_boundary_membership_is_start_inclusive_end_exclusive(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "boundary.bed"
    _write_bed(
        bed,
        [
            ("1", 99, 100),  # includes bp=100 (pos0=99)
            ("1", 199, 199),  # empty interval, includes nothing
            ("1", 399, 400),  # includes bp=400 (pos0=399)
            ("1", 500, 501),  # excludes bp=500 (pos0=499)
        ],
    )

    a = load_annotations([bed], reference)
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), np.array([1, 0, 0, 1, 0, 0, 0, 0], dtype=np.uint8))


def test_adjacent_interval_merge_matches_premerged_intervals(tmp_path):
    reference = load_reference(SHARDED_REF)
    merged = tmp_path / "merged.bed"
    split = tmp_path / "split.bed"
    _write_bed(merged, [("1", 99, 200), ("1", 299, 400)])
    _write_bed(split, [("1", 99, 150), ("1", 150, 200), ("1", 299, 350), ("1", 350, 400)])

    a_merged = load_annotations([merged], reference)
    a_split = load_annotations([split], reference)
    np.testing.assert_array_equal(a_merged.annomat.toarray(), a_split.annomat.toarray())


def test_overlapping_binary_bed_intervals_paint_union(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "overlap.bed"
    _write_bed(
        bed,
        [
            ("1", 99, 150),
            ("1", 120, 200),
            ("1", 299, 350),
        ],
    )

    a = load_annotations(bed, reference)

    np.testing.assert_array_equal(
        a.annomat.toarray().reshape(-1),
        np.array([1, 1, 1, 0, 0, 0, 0, 0], dtype=np.uint8),
    )


def test_chromosome_absent_in_bed_yields_zero_mask(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "chr2_only.bed"
    _write_bed(bed, [("2", 0, 1000)])

    a = load_annotations([bed], reference)
    np.testing.assert_array_equal(a.annomat.toarray(), np.zeros((reference.num_snp, 1), dtype=np.uint8))


def test_annotation_rejects_ambiguous_chr_labels(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "chr_prefixed.bed"
    _write_text(bed, "chr1\t99\t100\n")

    with pytest.raises(ValueError, match="chr-style labels"):
        load_annotations(bed, reference)


def test_load_annotations_metadata_sidecars_and_direct_metadata(tmp_path):
    reference = load_reference(SHARDED_REF)
    sidecar = tmp_path / "anno1.meta"
    sidecar.write_text('{"name":"anno1"}\nsecond line', encoding="utf-8")

    a = load_annotations(
        [ANNO1, ANNO2],
        reference,
        annotation_metadata_paths=[sidecar, None],
    )
    assert a.annotation_metadata[0] == '{"name":"anno1"}\nsecond line'
    assert json.loads(a.annotation_metadata[1])["source_file"] == str(ANNO2)

    b = load_annotations([ANNO1, ANNO2], reference, annotation_metadata=["m1", "m2"])
    assert b.annotation_metadata.tolist() == ["m1", "m2"]

    with pytest.raises(ValueError, match="at most one"):
        load_annotations(
            ANNO1,
            reference,
            annotation_metadata=["m"],
            annotation_metadata_paths=[sidecar],
        )


def test_load_annotation_numeric_headerless_single_column(tmp_path):
    reference = load_reference(SHARDED_REF)
    path = tmp_path / "score.annot"
    _write_text(
        path,
        "1\t99\t200\t0.5\n"
        "1\t299\t400\t-1.25\n"
        "X\t99\t101\t2.0\n",
    )

    a = load_annotation(path, reference)

    assert list(a.annonames) == ["score"]
    np.testing.assert_array_equal(a.is_binary, np.array([False]))
    dense = a.annomat.toarray().reshape(-1)
    np.testing.assert_allclose(dense, np.array([0.5, 0.5, -1.25, -1.25, 0.0, 2.0, 0.0, 0.0]))
    meta = json.loads(a.annotation_metadata[0])
    assert meta["source_column0"] == 3
    assert meta["source_column_name"] is None


def test_load_annotation_headered_multi_column_selection_and_sidecar(tmp_path):
    reference = load_reference(SHARDED_REF)
    path = tmp_path / "wide.annot"
    _write_text(
        path,
        "chrom\tstart0\tend0\tscore\tweight\n"
        "1\t99\t200\t0.5\t10\n"
        "1\t299\t400\t1.5\t20\n"
        "X\t99\t301\t2.5\t30\n",
    )
    sidecar = tmp_path / "wide.meta"
    sidecar.write_text("chrom meta\nstart meta\nend meta\nscore meta\nweight meta", encoding="utf-8")

    a = load_annotation(
        path,
        reference,
        header=True,
        value_columns=["weight", "score"],
        annotation_metadata_path=sidecar,
    )

    assert list(a.annonames) == ["weight", "score"]
    assert a.annotation_metadata.tolist() == ["weight meta", "score meta"]
    np.testing.assert_array_equal(a.is_binary, np.array([False, False]))
    np.testing.assert_allclose(
        a.annomat.toarray(),
        np.array(
            [
                [10.0, 0.5],
                [10.0, 0.5],
                [20.0, 1.5],
                [20.0, 1.5],
                [0.0, 0.0],
                [30.0, 2.5],
                [30.0, 2.5],
                [30.0, 2.5],
            ]
        ),
    )


def test_load_annotation_integer_value_columns_python_zero_based(tmp_path):
    reference = load_reference(SHARDED_REF)
    path = tmp_path / "integer_columns.annot"
    _write_text(
        path,
        "1\t99\t200\t0.5\t10\n"
        "1\t299\t400\t1.5\t20\n"
        "X\t99\t301\t2.5\t30\n",
    )

    a = load_annotation(
        path,
        reference,
        value_columns=[4, 3],
        annotation_names=["weight", "score"],
    )

    assert list(a.annonames) == ["weight", "score"]
    np.testing.assert_array_equal(a.is_binary, np.array([False, False]))
    np.testing.assert_allclose(
        a.annomat.toarray(),
        np.array(
            [
                [10.0, 0.5],
                [10.0, 0.5],
                [20.0, 1.5],
                [20.0, 1.5],
                [0.0, 0.0],
                [30.0, 2.5],
                [30.0, 2.5],
                [30.0, 2.5],
            ]
        ),
    )
    assert json.loads(a.annotation_metadata[0])["source_column0"] == 4
    assert json.loads(a.annotation_metadata[1])["source_column0"] == 3


def test_load_annotation_binary_sidecar_and_name_override(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "raw.bed"
    meta = tmp_path / "raw.meta"
    _write_text(bed, "1\t99\t200\n")
    meta.write_text("stored\nmetadata\n", encoding="utf-8")

    a = load_annotation(
        bed,
        reference,
        annotation_names=["renamed"],
        annotation_metadata_path=meta,
    )

    assert list(a.annonames) == ["renamed"]
    np.testing.assert_array_equal(a.is_binary, np.array([True]))
    assert a.annotation_metadata.tolist() == ["stored\nmetadata\n"]
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), np.array([1, 1, 0, 0, 0, 0, 0, 0]))


def test_load_annotation_validation_errors(tmp_path):
    reference = load_reference(SHARDED_REF)

    overlap = tmp_path / "overlap.annot"
    _write_text(overlap, "1\t99\t200\t1\n1\t150\t250\t2\n")
    with pytest.raises(ValueError, match="numeric annotation intervals overlap"):
        load_annotation(overlap, reference)

    nonfinite = tmp_path / "nonfinite.annot"
    _write_text(nonfinite, "1\t99\t200\tNaN\n")
    with pytest.raises(ValueError, match="finite numeric"):
        load_annotation(nonfinite, reference)

    wide = tmp_path / "wide.annot"
    _write_text(wide, "1\t99\t200\t1\t2\n")
    with pytest.raises(ValueError, match="requires explicit value_columns"):
        load_annotation(wide, reference)
    with pytest.raises(ValueError, match="requires annotation_names"):
        load_annotation(wide, reference, value_columns=[3])

    with pytest.raises(ValueError, match="named value_columns are invalid"):
        load_annotation(wide, reference, value_columns=["score"], annotation_names=["score"])

    one = tmp_path / "one.annot"
    _write_text(one, "1\t99\t200\t1\n")
    with pytest.raises(ValueError, match="string vector"):
        load_annotation(one, reference, annotation_metadata="scalar")

    bad_sidecar = tmp_path / "bad.meta"
    bad_sidecar.write_text("c1\nc2\nc3", encoding="utf-8")
    with pytest.raises(ValueError, match="line count mismatch"):
        load_annotation(one, reference, annotation_metadata_path=bad_sidecar)


def test_select_annotations_preserves_order_and_rejects_unknown():
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)

    sel = a.select_annotations(["anno2", "anno1"])
    assert list(sel.annonames) == ["anno2", "anno1"]
    np.testing.assert_array_equal(sel.is_binary, np.array([True, True]))
    np.testing.assert_array_equal(sel.annomat.toarray(), EXPECTED_MASK[:, [1, 0]])

    with pytest.raises(ValueError, match="unknown annotation"):
        a.select_annotations(["anno3"])
    with pytest.raises(ValueError, match="names must be unique"):
        a.select_annotations(["anno1", "anno1"])


def test_union_annotations_enforces_compatibility_and_name_collisions():
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1], reference)
    b = load_annotations([ANNO2], reference)

    u = a.union_annotations(b)
    assert list(u.annonames) == ["anno1", "anno2"]
    np.testing.assert_array_equal(u.is_binary, np.array([True, True]))
    np.testing.assert_array_equal(u.annomat.toarray(), EXPECTED_MASK)

    with pytest.raises(ValueError, match="name collision"):
        a.union_annotations(a)

    x_ref = reference.select_shards(["X"])
    x_only = create_annotation(x_ref, np.ones(x_ref.num_snp), "xonly")
    with pytest.raises(ValueError, match="shard count mismatch"):
        a.union_annotations(x_only)

    with pytest.raises(ValueError, match="requires another AnnotationPanel"):
        a.union_annotations({"annonames": ["other"], "shards": a.shards})

    other_bad_checksum = AnnotationPanel(
        [
            AnnotationShard._from_arrays(s.label, "deadbeef" * 4, s.annomat[:, :1])
            for s in a.shards
        ],
        ["other"],
        is_binary=np.array([True]),
        annotation_metadata=np.array([""], dtype=object),
    )
    with pytest.raises(ValueError, match="reference_checksum mismatch"):
        a.union_annotations(other_bad_checksum)


def test_create_annotations_and_create_annotation_validation():
    reference = load_reference(SHARDED_REF)

    panel = create_annotations(
        reference,
        EXPECTED_MASK,
        ["anno1", "anno2"],
        annotation_metadata=["m1", "m2"],
    )
    np.testing.assert_array_equal(panel.annomat.toarray(), EXPECTED_MASK)
    np.testing.assert_array_equal(panel.is_binary, np.array([True, True]))
    assert panel.annotation_metadata.tolist() == ["m1", "m2"]

    continuous = EXPECTED_MASK.astype(float)
    continuous[0, 0] = 2.5
    cont_panel = create_annotations(
        reference,
        continuous,
        ["cont1", "cont2"],
        is_binary=[False, True],
    )
    np.testing.assert_array_equal(cont_panel.is_binary, np.array([False, True]))
    np.testing.assert_allclose(cont_panel.annomat.toarray(), continuous)

    with pytest.raises(ValueError, match="shape mismatch"):
        create_annotations(reference, EXPECTED_MASK[:, :1], ["anno1", "anno2"])
    with pytest.raises(ValueError, match="non-binary"):
        create_annotations(
            reference,
            np.where(EXPECTED_MASK == 1, 2, 0),
            ["anno1", "anno2"],
            is_binary=[True, True],
        )
    with pytest.raises(ValueError, match="must be unique"):
        create_annotations(reference, EXPECTED_MASK, ["dup", "dup"])
    with pytest.raises(ValueError, match="annotation_metadata length mismatch"):
        create_annotations(reference, EXPECTED_MASK, ["anno1", "anno2"], annotation_metadata=["one"])

    col = EXPECTED_MASK[:, 0]
    single = create_annotation(reference, col, "anno1", annotation_metadata="single meta")
    np.testing.assert_array_equal(single.annomat.toarray(), col.reshape(-1, 1))
    assert single.annotation_metadata.tolist() == ["single meta"]

    with pytest.raises(ValueError, match="annovec length mismatch"):
        create_annotation(reference, np.array([1, 0], dtype=np.uint8), "bad")
    with pytest.raises(ValueError, match="must be binary"):
        create_annotation(reference, np.r_[2, np.zeros(reference.num_snp - 1)], "bad", is_binary=True)


def test_create_annotation_can_represent_all_snps():
    reference = load_reference(SHARDED_REF)
    panel = create_annotation(reference, np.ones(reference.num_snp, dtype=np.uint8), "all_snps")
    assert list(panel.annonames) == ["all_snps"]
    np.testing.assert_array_equal(panel.annomat.toarray(), np.ones((reference.num_snp, 1), dtype=np.uint8))

    panel_named = create_annotation(reference, np.ones(reference.num_snp, dtype=np.uint8), "custom_all")
    assert list(panel_named.annonames) == ["custom_all"]
    np.testing.assert_array_equal(panel_named.annomat.toarray(), np.ones((reference.num_snp, 1), dtype=np.uint8))


def test_create_annotations_sparse_explicit_zeros_do_not_flip_to_one():
    reference = load_reference(SHARDED_REF)
    # Handcrafted CSR with an explicit zero stored at (0, 0) and a true one at (1, 0).
    # Stored zeros must remain logical zeros after normalization.
    annomat = sparse.csr_matrix(
        (
            np.array([0.0, 1.0], dtype=float),
            np.array([0, 0], dtype=np.int32),
            np.array([0, 1, 2, 2, 2, 2, 2, 2, 2], dtype=np.int32),
        ),
        shape=(reference.num_snp, 1),
    )

    panel = create_annotations(reference, annotation_matrix=annomat, annotation_names=["anno"])
    dense = panel.annomat.toarray().reshape(-1)
    assert dense[0] == 0
    assert dense[1] == 1
    assert np.count_nonzero(dense) == 1


def test_cache_roundtrip_subset_and_compatibility(tmp_path):
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)
    cache = tmp_path / "annotations_cache.npz"

    a.save_cache(cache)
    loaded = load_annotations_cache(cache)

    assert sparse.isspmatrix_csr(loaded.annomat)
    np.testing.assert_array_equal(loaded.annomat.toarray(), EXPECTED_MASK)
    np.testing.assert_array_equal(loaded.is_binary, np.array([True, True]))
    assert all(json.loads(x)["source_file"] for x in loaded.annotation_metadata)
    assert reference.is_object_compatible(loaded) is True

    x_only = load_annotations_cache(cache, shards=["X"])
    assert [sh.label for sh in x_only.shards] == ["X"]
    np.testing.assert_array_equal(x_only.annomat.toarray(), EXPECTED_MASK[5:, :])


def test_cache_validation_and_post_load_compatibility_check(tmp_path):
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)
    cache = tmp_path / "annotations_cache.npz"
    save_annotations_cache(a, cache)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    assert meta["schema"] == "annotations_cache/0.2"
    assert meta["is_binary"] == [True, True]
    assert len(meta["annotation_metadata"]) == 2

    meta["schema"] = "annotations_cache/bad"
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_schema = tmp_path / "bad_schema.npz"
    np.savez_compressed(bad_schema, **arrays)
    with pytest.raises(ValueError, match="Unsupported annotations cache schema"):
        load_annotations_cache(bad_schema)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["shard_checksums"] = meta["shard_checksums"][:-1]
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_len = tmp_path / "bad_len.npz"
    np.savez_compressed(bad_len, **arrays)
    with pytest.raises(ValueError, match="shard_labels and shard_checksums length mismatch"):
        load_annotations_cache(bad_len)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["shard_checksums"][0] = "deadbeef" * 4
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_chk = tmp_path / "bad_chk.npz"
    np.savez_compressed(bad_chk, **arrays)
    loaded_bad = load_annotations_cache(bad_chk)
    assert reference.is_object_compatible(loaded_bad) is False

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["schema"] = "annotations_cache/0.1"
    meta.pop("is_binary")
    meta.pop("annotation_metadata")
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    old_cache = tmp_path / "old_cache.npz"
    np.savez_compressed(old_cache, **arrays)
    loaded_old = load_annotations_cache(old_cache)
    np.testing.assert_array_equal(loaded_old.is_binary, np.array([True, True]))
    assert loaded_old.annotation_metadata.tolist() == ["", ""]
    assert loaded_old.annomat.dtype == np.float64


def test_cross_language_cache_not_supported(tmp_path):
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)

    py_cache = tmp_path / "annotations_cache_py.npz"
    save_annotations_cache(a, py_cache)
    octave_result = run_octave(
        _octave_script(
            f"try; statgen.load_annotations_cache('{py_cache}'); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
        )
    )
    assert octave_result.returncode == 0, octave_result.stderr
    assert octave_result.stdout.strip() == "FAIL"

    mat_cache = tmp_path / "annotations_cache_mat.mat"
    octave_result = run_octave(
        _octave_script(
            f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
            f"a = statgen.load_annotations({{[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}}, ref); "
            f"statgen.save_annotations_cache(a, '{mat_cache}'); "
            "fprintf('OK\\n');"
        )
    )
    assert octave_result.returncode == 0, octave_result.stderr
    assert octave_result.stdout.strip() == "OK"
    with pytest.raises(Exception):
        load_annotations_cache(mat_cache)


def _octave_script(expr: str) -> str:
    matlab_dir = str(MATLAB_DIR)
    fixture_dir = str(FIXTURES_DIR)
    return f"addpath('{matlab_dir}'); fixture_dir = '{fixture_dir}'; " + expr


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotations_fixture_masks_and_names():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        "fprintf('%d\\n', a.num_snp); "
        "fprintf('%d\\n', a.num_annot); "
        "fprintf('%s,%s\\n', a.annonames{1}, a.annonames{2}); "
        "fprintf('%d,%d\\n', nnz(a.annomat(:,1)), nnz(a.annomat(:,2))); "
        "M = full(a.annomat); "
        "fprintf('%d', M(1,1)); fprintf('%d', M(5,1)); fprintf('%d', M(6,1)); fprintf('%d', M(6,2)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "8"
    assert lines[1] == "2"
    assert lines[2] == "anno1,anno2"
    assert lines[3] == "5,3"
    assert lines[4] == "1101"


@pytest.mark.octave
@skipif_no_octave
def test_octave_annotations_cache_roundtrip_and_subset(tmp_path):
    cache = tmp_path / "annotations_cache.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        f"statgen.save_annotations_cache(a, '{cache}'); "
        f"b = statgen.load_annotations_cache('{cache}'); "
        f"x = statgen.load_annotations_cache('{cache}', {{'X'}}); "
        "fprintf('%d\\n', ref.is_object_compatible(b)); "
        "fprintf('%d\\n', b.num_snp); "
        "fprintf('%d\\n', x.num_snp); "
        "fprintf('%s\\n', x.shards{1}.label); "
        f"L = load('{cache}'); "
        "fprintf('%d %d %d %d %d\\n', isfield(L, 'metadata'), isfield(L, 'annomat'), isfield(L, 'is_binary'), isfield(L, 'annotation_metadata'), isfield(L, 'cache_shards'));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
    assert lines[1] == "8"
    assert lines[2] == "3"
    assert lines[4] == "1 1 1 1 0"
    assert lines[3] == "X"


@pytest.mark.octave
@skipif_no_octave
def test_octave_annotations_cache_roundtrip_preserves_continuous_metadata(tmp_path):
    annot = tmp_path / "score_cache.annot"
    cache = tmp_path / "mixed_annotations_cache.mat"
    _write_text(
        annot,
        "1\t99\t200\t0.5\n"
        "1\t299\t400\t-1.25\n"
        "X\t99\t301\t2.0\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotation('{annot}', ref); "
        f"statgen.save_annotations_cache(a, '{cache}'); "
        f"b = statgen.load_annotations_cache('{cache}'); "
        "M = full(b.annomat); "
        "fprintf('%d\\n', b.is_binary(1)); "
        "fprintf('%d\\n', ~isempty(strfind(b.annotation_metadata{1}, 'source_column0'))); "
        "fprintf('%d\\n', ~isempty(strfind(b.annotation_metadata{1}, 'score_cache.annot'))); "
        "fprintf('%.2f %.2f %.2f\\n', M(1,1), M(3,1), M(8,1));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "0"
    assert lines[1] == "1"
    assert lines[2] == "1"
    assert lines[3] == "0.50 -1.25 2.00"


@pytest.mark.octave
@skipif_no_octave
def test_octave_annotations_cache_rejects_zero_shards(tmp_path):
    cache = tmp_path / "annotations_cache.mat"
    bad_cache = tmp_path / "annotations_zero_shards.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        f"statgen.save_annotations_cache(a, '{cache}'); "
        f"L = load('{cache}'); "
        "L.metadata.n_shards = 0; "
        "L.metadata.shard_labels = {}; "
        "L.metadata.shard_checksums = {}; "
        "L.metadata.shard_start0 = []; "
        "L.metadata.shard_stop0 = []; "
        f"save('{bad_cache}', '-struct', 'L'); "
        f"ok = 0; try; statgen.load_annotations_cache('{bad_cache}'); catch ME; ok = ~isempty(strfind(ME.message, 'n_shards must be positive')); end; "
        "fprintf('%d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_empty_bed_fails(tmp_path):
    empty = tmp_path / "empty.bed"
    empty.write_text("")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations('{empty}', ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_create_annotation_can_represent_all_snps():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.create_annotation(ref, ones(ref.num_snp,1), 'all_snps'); "
        "b = statgen.create_annotation(ref, ones(ref.num_snp,1), 'custom_all'); "
        "fprintf('%s\\n', a.annonames{1}); "
        "fprintf('%d\\n', nnz(a.annomat)); "
        "fprintf('%s\\n', b.annonames{1}); "
        "fprintf('%d\\n', nnz(b.annomat));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "all_snps"
    assert lines[1] == "8"
    assert lines[2] == "custom_all"
    assert lines[3] == "8"


@pytest.mark.octave
@skipif_no_octave
def test_octave_duplicate_bed_basenames_fail(tmp_path):
    p1 = tmp_path / "set1" / "dup.bed"
    p2 = tmp_path / "set2" / "dup.bed"
    _write_bed(p1, [("1", 99, 100)])
    _write_bed(p2, [("X", 99, 100)])

    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations({{'{p1}','{p2}'}}, ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotations_metadata_sidecars_and_generated_fallback(tmp_path):
    sidecar = tmp_path / "anno1.meta"
    sidecar.write_text('{"name":"anno1"}\nsecond line', encoding="utf-8")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations({{[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}}, ref, 'annotation_metadata_paths', {{'{sidecar}', ''}}); "
        "fprintf('%d\\n', strcmp(a.annotation_metadata{1}, sprintf('{\"name\":\"anno1\"}\\nsecond line'))); "
        "fprintf('%d\\n', ~isempty(strfind(a.annotation_metadata{2}, 'source_file'))); "
        "fprintf('%d\\n', ~isempty(strfind(a.annotation_metadata{2}, 'anno2.bed'))); "
        "b = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref, 'annotation_metadata', {'m1','m2'}); "
        "fprintf('%s,%s\\n', b.annotation_metadata{1}, b.annotation_metadata{2});"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
    assert lines[1] == "1"
    assert lines[2] == "1"
    assert lines[3] == "m1,m2"


@pytest.mark.octave
@skipif_no_octave
def test_octave_bed_metadata_and_comment_lines_are_ignored(tmp_path):
    bed = tmp_path / "with_headers_octave.bed"
    _write_text(
        bed,
        "# track name=CodingExon description=test\n"
        "# browser position chr1:1-1000\n"
        "# comment\n"
        "\n"
        "1\t99\t200\n"
        "1\t299\t400\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{bed}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%d', v(1)); fprintf('%d', v(2)); fprintf('%d', v(3)); fprintf('%d', v(4)); "
        "fprintf('%d', v(5)); fprintf('%d', v(6)); fprintf('%d', v(7)); fprintf('%d', v(8)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "11110000"


@pytest.mark.octave
@skipif_no_octave
def test_octave_metadata_only_bed_file_fails(tmp_path):
    bed = tmp_path / "metadata_only_octave.bed"
    _write_text(
        bed,
        "# track name=foo\n"
        "# browser position chr1:1-100\n"
        "# note\n"
        "\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations('{bed}', ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_non_tab_delimited_bed_rows_fail(tmp_path):
    bed = tmp_path / "space_delimited_octave.bed"
    _write_text(bed, "1  99   200\n1  299  400\n")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations('{bed}', ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_boundary_membership_start_inclusive_end_exclusive(tmp_path):
    bed = tmp_path / "boundary_octave.bed"
    _write_bed(
        bed,
        [
            ("1", 99, 100),
            ("1", 199, 199),
            ("1", 399, 400),
            ("1", 500, 501),
        ],
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{bed}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%d', v(1)); fprintf('%d', v(2)); fprintf('%d', v(3)); fprintf('%d', v(4)); "
        "fprintf('%d', v(5)); fprintf('%d', v(6)); fprintf('%d', v(7)); fprintf('%d', v(8)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "10010000"


@pytest.mark.octave
@skipif_no_octave
def test_octave_overlapping_binary_bed_intervals_paint_union(tmp_path):
    bed = tmp_path / "overlap_octave.bed"
    _write_bed(bed, [("1", 99, 150), ("1", 120, 200), ("1", 299, 350)])
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{bed}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%d', v(1)); fprintf('%d', v(2)); fprintf('%d', v(3)); fprintf('%d', v(4)); "
        "fprintf('%d', v(5)); fprintf('%d', v(6)); fprintf('%d', v(7)); fprintf('%d', v(8)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "11100000"


@pytest.mark.octave
@skipif_no_octave
def test_octave_adjacent_interval_merge_matches_premerged(tmp_path):
    merged = tmp_path / "merged_octave.bed"
    split = tmp_path / "split_octave.bed"
    _write_bed(merged, [("1", 99, 200), ("1", 299, 400)])
    _write_bed(split, [("1", 99, 150), ("1", 150, 200), ("1", 299, 350), ("1", 350, 400)])
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{merged}', ref); "
        f"b = statgen.load_annotations('{split}', ref); "
        "eq = isequal(full(a.annomat), full(b.annomat)); "
        "fprintf('%d\\n', eq);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotation_numeric_headerless_single_column(tmp_path):
    path = tmp_path / "score.annot"
    _write_text(
        path,
        "1\t99\t200\t0.5\n"
        "1\t299\t400\t-1.25\n"
        "X\t99\t301\t2.0\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotation('{path}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%s\\n', a.annonames{1}); "
        "fprintf('%d\\n', a.is_binary(1)); "
        "fprintf('%.2f %.2f %.2f %.2f\\n', v(1), v(3), v(6), v(8));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "score"
    assert lines[1] == "0"
    assert lines[2] == "0.50 -1.25 2.00 2.00"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotation_headered_multi_column_selection_and_sidecar(tmp_path):
    path = tmp_path / "wide.annot"
    _write_text(
        path,
        "chrom\tstart0\tend0\tscore\tweight\n"
        "1\t99\t200\t0.5\t10\n"
        "1\t299\t400\t1.5\t20\n"
        "X\t99\t301\t2.5\t30\n",
    )
    sidecar = tmp_path / "wide.meta"
    sidecar.write_text("chrom meta\nstart meta\nend meta\nscore meta\nweight meta", encoding="utf-8")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotation('{path}', ref, 'header', true, 'value_columns', {{'weight','score'}}, 'annotation_metadata_path', '{sidecar}'); "
        "M = full(a.annomat); "
        "fprintf('%s,%s\\n', a.annonames{1}, a.annonames{2}); "
        "fprintf('%s,%s\\n', a.annotation_metadata{1}, a.annotation_metadata{2}); "
        "fprintf('%d,%d\\n', a.is_binary(1), a.is_binary(2)); "
        "fprintf('%.1f %.1f %.1f %.1f\\n', M(1,1), M(1,2), M(8,1), M(8,2));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "weight,score"
    assert lines[1] == "weight meta,score meta"
    assert lines[2] == "0,0"
    assert lines[3] == "10.0 0.5 30.0 2.5"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotation_matlab_one_based_value_columns(tmp_path):
    path = tmp_path / "integer_columns.annot"
    _write_text(
        path,
        "1\t99\t200\t0.5\t10\n"
        "1\t299\t400\t1.5\t20\n"
        "X\t99\t301\t2.5\t30\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotation('{path}', ref, 'value_columns', [5 4], 'annotation_names', {{'weight','score'}}); "
        "M = full(a.annomat); "
        "fprintf('%s,%s\\n', a.annonames{1}, a.annonames{2}); "
        "fprintf('%.1f %.1f %.1f %.1f\\n', M(1,1), M(1,2), M(8,1), M(8,2));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "weight,score"
    assert lines[1] == "10.0 0.5 30.0 2.5"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotation_validation_errors(tmp_path):
    overlap = tmp_path / "overlap.annot"
    _write_text(overlap, "1\t99\t200\t1\n1\t150\t250\t2\n")
    nonfinite = tmp_path / "nonfinite.annot"
    _write_text(nonfinite, "1\t99\t200\tNaN\n")
    wide = tmp_path / "wide.annot"
    _write_text(wide, "1\t99\t200\t1\t2\n")
    sidecar = tmp_path / "bad.meta"
    sidecar.write_text("c1\nc2\nc3", encoding="utf-8")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotation('{overlap}', ref); fprintf('NOFAIL1\\n'); catch; fprintf('FAIL1\\n'); end; "
        f"try; statgen.load_annotation('{nonfinite}', ref); fprintf('NOFAIL2\\n'); catch; fprintf('FAIL2\\n'); end; "
        f"try; statgen.load_annotation('{wide}', ref); fprintf('NOFAIL3\\n'); catch; fprintf('FAIL3\\n'); end; "
        f"try; statgen.load_annotation('{wide}', ref, 'value_columns', 4); fprintf('NOFAIL4\\n'); catch; fprintf('FAIL4\\n'); end; "
        f"try; statgen.load_annotation('{wide}', ref, 'annotation_metadata_path', '{sidecar}'); fprintf('NOFAIL5\\n'); catch; fprintf('FAIL5\\n'); end;"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines() == ["FAIL1", "FAIL2", "FAIL3", "FAIL4", "FAIL5"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_union_annotations_happy_path_and_name_collision():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations([fixture_dir '/annotations/anno1.bed'], ref); "
        "b = statgen.load_annotations([fixture_dir '/annotations/anno2.bed'], ref); "
        "u = a.union_annotations(b); "
        "fprintf('%s,%s\\n', u.annonames{1}, u.annonames{2}); "
        "fprintf('%d,%d\\n', nnz(u.annomat(:,1)), nnz(u.annomat(:,2))); "
        "try; a.union_annotations(a); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "anno1,anno2"
    assert lines[1] == "5,3"
    assert lines[2] == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_union_annotations_preserves_continuous_columns_and_rejects_non_panels():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations([fixture_dir '/annotations/anno1.bed'], ref); "
        "cont = statgen.create_annotation(ref, 2 * ones(ref.num_snp, 1), 'continuous', false, 'cont meta'); "
        "u = a.union_annotations(cont); "
        "fprintf('%s,%s\\n', u.annonames{1}, u.annonames{2}); "
        "fprintf('%d,%d\\n', u.is_binary(1), u.is_binary(2)); "
        "fprintf('%s\\n', u.annotation_metadata{2}); "
        "M = full(u.annomat); fprintf('%.1f\\n', M(1,2)); "
        "fallback = struct(); "
        "fallback.annonames = {'fallback'}; "
        "fallback.is_binary = false; "
        "fallback.annotation_metadata = {'fallback meta'}; "
        "fallback.shards = a.shards; "
        "try; a.union_annotations(fallback); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end;"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "anno1,continuous"
    assert lines[1] == "1,0"
    assert lines[2] == "cont meta"
    assert lines[3] == "2.0"
    assert lines[4] == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_create_annotations_and_create_annotation_validation():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "A = [ones(ref.num_snp,1), zeros(ref.num_snp,1)]; "
        "ok = statgen.create_annotations(ref, A, {'a','b'}); "
        "fprintf('%d\\n', ok.num_annot); "
        "fprintf('%d,%d\\n', ok.is_binary(1), ok.is_binary(2)); "
        "try; statgen.create_annotations(ref, ones(ref.num_snp,1), {'a','b'}); fprintf('NOFAIL1\\n'); catch; fprintf('FAIL1\\n'); end; "
        "cont = statgen.create_annotations(ref, 2*ones(ref.num_snp,2), {'c','d'}, [false; false]); fprintf('%d,%d\\n', cont.is_binary(1), cont.is_binary(2)); "
        "try; statgen.create_annotations(ref, 2*ones(ref.num_snp,2), {'a','b'}, [true; true]); fprintf('NOFAIL2\\n'); catch; fprintf('FAIL2\\n'); end; "
        "try; statgen.create_annotations(ref, A, {'dup','dup'}); fprintf('NOFAIL3\\n'); catch; fprintf('FAIL3\\n'); end; "
        "single = statgen.create_annotation(ref, ones(ref.num_snp,1), 'single'); fprintf('%d\\n', single.num_annot); "
        "try; statgen.create_annotation(ref, ones(2,1), 'bad'); fprintf('NOFAIL4\\n'); catch; fprintf('FAIL4\\n'); end; "
        "v = ones(ref.num_snp,1); v(1)=2; try; statgen.create_annotation(ref, v, 'bad', true); fprintf('NOFAIL5\\n'); catch; fprintf('FAIL5\\n'); end;"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "2"
    assert lines[1] == "1,1"
    assert lines[2] == "FAIL1"
    assert lines[3] == "0,0"
    assert lines[4] == "FAIL2"
    assert lines[5] == "FAIL3"
    assert lines[6] == "1"
    assert lines[7] == "FAIL4"
    assert lines[8] == "FAIL5"


@pytest.mark.octave
@skipif_no_octave
def test_octave_select_annotations_and_unknown_name_error():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        "b = a.select_annotations({'anno2','anno1'}); "
        "fprintf('%s,%s\\n', b.annonames{1}, b.annonames{2}); "
        "try; a.select_annotations({'missing'}); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "anno2,anno1"
    assert lines[1] == "FAIL"
