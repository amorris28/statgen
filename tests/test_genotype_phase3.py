import json

import numpy as np
import pytest

from statgen.genotype import load_genotype, load_genotype_cache, save_genotype_cache
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"


def _load_cache_arrays(path):
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def _write_cache(path, arrays, meta):
    arrays = dict(arrays)
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    np.savez_compressed(path, **arrays)


def _copy_sharded_genotype(dst):
    for label in ("1", "X"):
        for suffix in (".bim", ".fam", ".bed", ".ploidy"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            if src.exists():
                (dst / f"{label}{suffix}").write_bytes(src.read_bytes())


def _write_zero_bed(path, *, num_snp, num_sample):
    bytes_per_snp = (num_sample + 3) // 4
    path.write_bytes(b"\x6c\x1b\x01" + b"\x00" * (num_snp * bytes_per_snp))


def test_genotype_cache_roundtrip_preserves_accessors_and_fetch(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)
    cache = tmp_path / "genotype_cache.npz"
    save_genotype_cache(panel, cache)

    loaded = load_genotype_cache(cache)

    assert loaded.source_layout == panel.source_layout
    assert [s.label for s in loaded.shards] == [s.label for s in panel.shards]
    assert loaded.num_snp == panel.num_snp
    assert loaded.num_sample == panel.num_sample
    assert np.array_equal(loaded.is_present, panel.is_present)
    np.testing.assert_allclose(loaded.ploidy_male, panel.ploidy_male, equal_nan=True)
    np.testing.assert_allclose(loaded.ploidy_female, panel.ploidy_female, equal_nan=True)
    assert np.array_equal(loaded.source_row0, panel.source_row0)
    assert np.array_equal(loaded.fid, panel.fid)
    assert np.array_equal(loaded.iid, panel.iid)
    assert np.array_equal(loaded.father_id, panel.father_id)
    assert np.array_equal(loaded.mother_id, panel.mother_id)
    assert np.array_equal(loaded.sex, panel.sex)
    assert np.array_equal(loaded.is_male, panel.is_male)
    assert np.array_equal(loaded.is_female, panel.is_female)
    assert np.array_equal(loaded.is_subject_present("X"), panel.is_subject_present("X"))

    requested = [5, 0, 4, 5, 2]
    assert np.array_equal(
        loaded.fetch_genotypes_int8(requested),
        panel.fetch_genotypes_int8(requested),
    )
    np.testing.assert_allclose(
        loaded.fetch_genotypes(requested),
        panel.fetch_genotypes(requested),
        equal_nan=True,
    )


def test_genotype_cache_load_is_lazy_for_missing_bed_and_override_fetches(tmp_path):
    _copy_sharded_genotype(tmp_path)
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(str(tmp_path / "@"), ref)
    cache = tmp_path / "genotype_cache.npz"
    save_genotype_cache(panel, cache)

    (tmp_path / "1.bed").unlink()
    (tmp_path / "X.bed").unlink()

    loaded = load_genotype_cache(cache)
    assert loaded.num_snp == panel.num_snp

    with pytest.raises(FileNotFoundError):
        loaded.fetch_genotypes_int8([0])

    geno = loaded.fetch_genotypes_int8(
        [0, 5],
        bed_path=str(FIXTURES_DIR / "genotype/sharded/@.bed"),
    )
    assert np.array_equal(geno, np.full((4, 2), 2, dtype=np.int8))


def test_genotype_cache_subset_preserves_panel_axis_and_chrx_mask(tmp_path):
    _copy_sharded_genotype(tmp_path)
    (tmp_path / "X.fam").write_text(
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t1\t-9\n"
    )
    _write_zero_bed(tmp_path / "X.bed", num_snp=3, num_sample=2)

    ref = load_reference(REF_SHARDED)
    panel = load_genotype(str(tmp_path / "@"), ref)
    cache = tmp_path / "genotype_cache.npz"
    save_genotype_cache(panel, cache)

    loaded_x = load_genotype_cache(cache, shards=["X"])

    assert [s.label for s in loaded_x.shards] == ["X"]
    assert loaded_x.source_layout == "sharded"
    assert loaded_x.num_snp == 3
    assert loaded_x.num_sample == 4
    assert np.array_equal(loaded_x.fid, panel.fid)
    assert loaded_x.is_subject_present("X").tolist() == [True, False, True, False]
    assert loaded_x.shards[0].source_subject_row0.tolist() == [1, -1, 0, -1]


def test_genotype_cache_validation_errors(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)
    cache = tmp_path / "genotype_cache.npz"
    save_genotype_cache(panel, cache)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["schema"] = "genotype_cache/bad"
    bad_schema = tmp_path / "bad_schema.npz"
    _write_cache(bad_schema, arrays, meta)
    with pytest.raises(ValueError, match="Unsupported genotype cache schema"):
        load_genotype_cache(bad_schema)

    arrays = _load_cache_arrays(cache)
    arrays["ploidy_male"] = arrays["ploidy_male"][:-1]
    meta = json.loads(bytes(arrays["_meta"]).decode())
    bad_lengths = tmp_path / "bad_lengths.npz"
    _write_cache(bad_lengths, arrays, meta)
    with pytest.raises(ValueError, match="ploidy vector lengths mismatch"):
        load_genotype_cache(bad_lengths)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["shard_start0"][1] = meta["shard_start0"][1] + 1
    bad_offsets = tmp_path / "bad_offsets.npz"
    _write_cache(bad_offsets, arrays, meta)
    with pytest.raises(ValueError, match="shard offsets"):
        load_genotype_cache(bad_offsets)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["source_layout"] = "flat"
    bad_layout = tmp_path / "bad_layout.npz"
    _write_cache(bad_layout, arrays, meta)
    with pytest.raises(ValueError, match="source_layout"):
        load_genotype_cache(bad_layout)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["bed_file_sizes"][0] = meta["bed_file_sizes"][0] + 1
    bad_bed_size = tmp_path / "bad_bed_size.npz"
    _write_cache(bad_bed_size, arrays, meta)
    with pytest.raises(ValueError, match="bed_file_size mismatch"):
        load_genotype_cache(bad_bed_size)

    arrays = _load_cache_arrays(cache)
    arrays.pop("source_row0")
    meta = json.loads(bytes(arrays["_meta"]).decode())
    missing_field = tmp_path / "missing_field.npz"
    _write_cache(missing_field, arrays, meta)
    with pytest.raises(ValueError, match="missing array field 'source_row0'"):
        load_genotype_cache(missing_field)

    arrays = _load_cache_arrays(cache)
    arrays["is_male"] = ~arrays["is_male"]
    meta = json.loads(bytes(arrays["_meta"]).decode())
    bad_sex_mask = tmp_path / "bad_sex_mask.npz"
    _write_cache(bad_sex_mask, arrays, meta)
    with pytest.raises(ValueError, match="is_male does not match sex"):
        load_genotype_cache(bad_sex_mask)

    with pytest.raises(ValueError, match="requested shard '2' is not present"):
        load_genotype_cache(cache, shards=["2"])
