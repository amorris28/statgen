import math

import numpy as np
import pytest

from statgen.genotype import load_genotype
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


_CHR1_CALLS = np.array(
    [
        [2, -1, 1, 0],
        [0, 1, 2, -1],
        [1, 1, 0, 2],
        [2, 0, -1, 1],
        [-1, 2, 0, 1],
    ],
    dtype=np.int8,
)


def _write_bed_calls(path, calls):
    calls = np.asarray(calls, dtype=np.int8)
    code = {
        2: 0b00,
        -1: 0b01,
        1: 0b10,
        0: 0b11,
    }
    bytes_per_snp = (calls.shape[1] + 3) // 4
    payload = bytearray()
    for row in calls:
        row_bytes = [0] * bytes_per_snp
        for j, value in enumerate(row.tolist()):
            row_bytes[j // 4] |= code[int(value)] << (2 * (j % 4))
        payload.extend(row_bytes)
    path.write_bytes(b"\x6c\x1b\x01" + bytes(payload))


def _copy_shard_files(dst, label, *, ploidy=True):
    for suffix in (".bim", ".fam", ".bed"):
        (dst / f"{label}{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/{label}{suffix}").read_bytes()
        )
    src_ploidy = FIXTURES_DIR / f"genotype/sharded/{label}.ploidy"
    if ploidy and src_ploidy.exists():
        (dst / f"{label}.ploidy").write_bytes(src_ploidy.read_bytes())


def test_fetch_committed_fixture_one_multi_cross_shard_and_repeated():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)

    one = panel.fetch_genotypes_int8([0])
    assert one.shape == (4, 1)
    assert one.tolist() == [[2], [2], [2], [2]]

    requested = [5, 0, 4, 5, 2]
    geno = panel.fetch_genotypes_int8(requested)
    assert geno.shape == (4, len(requested))
    assert np.array_equal(geno, np.full((4, len(requested)), 2, dtype=np.int8))


def test_fetch_decodes_all_two_bit_states_and_preserves_order(tmp_path):
    for suffix in (".bim", ".fam"):
        (tmp_path / f"chr1{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes()
        )
    _write_bed_calls(tmp_path / "chr1.bed", _CHR1_CALLS)

    ref = load_reference(REF_SHARDED, shards=["1"])
    panel = load_genotype(tmp_path / "chr1", ref)

    requested = [1, 0, 1, 4]
    expected = _CHR1_CALLS[requested].T
    geno_int8 = panel.fetch_genotypes_int8(requested)
    assert np.array_equal(geno_int8, expected)

    geno = panel.fetch_genotypes(requested)
    assert geno.dtype == np.float64
    assert math.isnan(geno[3, 0])
    expected_float = expected.astype(np.float64)
    expected_float[expected == -1] = np.nan
    assert np.array_equal(geno, expected_float, equal_nan=True)


def test_fetch_chrx_subset_expands_to_panel_sample_axis(tmp_path):
    _copy_shard_files(tmp_path, "1")
    _copy_shard_files(tmp_path, "X")
    (tmp_path / "X.fam").write_text(
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t1\t-9\n"
    )
    _write_bed_calls(
        tmp_path / "X.bed",
        np.array(
            [
                [0, 1],
                [2, -1],
                [1, 0],
            ],
            dtype=np.int8,
        ),
    )

    ref = load_reference(REF_SHARDED)
    panel = load_genotype(str(tmp_path / "@"), ref)

    geno = panel.fetch_genotypes_int8([5])
    assert geno[:, 0].tolist() == [1, -1, 0, -1]


def test_absent_snp_fails_before_reading_bed_payload(tmp_path):
    extra_bim = (
        "1\trs1001\t0\t100\tA\tG\n"
        "1\trs1002\t0\t200\tC\tT\n"
        "1\trs1003\t0\t300\tA\tC\n"
        "1\trs1004\t0\t400\tG\tA\n"
        "1\trs1005\t0\t500\tT\tC\n"
        "1\trs1006\t0\t600\tA\tC\n"
    )
    (tmp_path / "extra.bim").write_text(extra_bim)
    ref = load_reference(tmp_path / "extra.bim")
    panel = load_genotype(G_SHARDED, ref)

    with pytest.raises(ValueError, match="not present in the genotype source"):
        panel.fetch_genotypes_int8([5], bed_path=tmp_path / "missing.bed")


def test_flat_bed_override_works_for_nonsharded_panel():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_NONSHARDED, ref)

    geno = panel.fetch_genotypes_int8(
        [0, 5],
        bed_path=FIXTURES_DIR / "genotype/nonsharded/all.bed",
    )
    assert np.array_equal(geno, np.full((4, 2), 2, dtype=np.int8))


def test_at_bed_override_works_for_sharded_panel_and_single_shard_subset():
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)
    pattern = str(FIXTURES_DIR / "genotype/sharded/@.bed")

    geno = panel.fetch_genotypes_int8([0, 5], bed_path=pattern)
    assert np.array_equal(geno, np.full((4, 2), 2, dtype=np.int8))

    chrx = panel.select_shards(["X"])
    x_geno = chrx.fetch_genotypes_int8([0], bed_path=pattern)
    assert np.array_equal(x_geno, np.full((4, 1), 2, dtype=np.int8))


def test_at_override_against_nonsharded_panel_fails_before_reading(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_NONSHARDED, ref)

    with pytest.raises(ValueError, match="@ override incompatible with non-sharded"):
        panel.fetch_genotypes_int8([0], bed_path=str(tmp_path / "@.bed"))


def test_flat_override_against_unequal_sharded_sources_fails_before_reading(tmp_path):
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)

    with pytest.raises(ValueError, match="flat bed_path override requires equal"):
        panel.fetch_genotypes_int8([0, 5], bed_path=tmp_path / "flat.bed")


def test_bed_override_size_and_header_validation(tmp_path):
    ref = load_reference(REF_SHARDED, shards=["1"])
    panel = load_genotype(G_SHARDED, ref)
    expected_size = 3 + 5

    bad_magic = tmp_path / "bad_magic.bed"
    bad_magic.write_bytes(b"\x00\x00\x01" + b"\x00" * (expected_size - 3))
    with pytest.raises(ValueError, match="invalid PLINK BED magic"):
        panel.fetch_genotypes_int8([0], bed_path=bad_magic)

    sample_major = tmp_path / "sample_major.bed"
    sample_major.write_bytes(b"\x6c\x1b\x00" + b"\x00" * (expected_size - 3))
    with pytest.raises(ValueError, match="SNP-major mode"):
        panel.fetch_genotypes_int8([0], bed_path=sample_major)

    trailing = tmp_path / "trailing.bed"
    trailing.write_bytes(b"\x6c\x1b\x01" + b"\x00" * (expected_size - 2))
    with pytest.raises(ValueError, match="BED file size mismatch"):
        panel.fetch_genotypes_int8([0], bed_path=trailing)

