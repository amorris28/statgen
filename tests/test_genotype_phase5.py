import numpy as np
import pytest

from tests.conftest import FIXTURES_DIR, matlab_data_lines, run_octave, skipif_no_octave


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


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_fixture_mixed_shard_repeated_and_float_wrapper():
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        "gi = g.fetch_genotypes_int8([6 1 5 6 3]); "
        "gf = g.fetch_genotypes([6 1 5 6 3]); "
        "fprintf('%d %d\\n', size(gi, 1), size(gi, 2)); "
        "fprintf('%d ', gi); fprintf('\\n'); "
        "fprintf('%s %.0f\\n', class(gf), sum(isnan(gf(:))));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "4 5",
        "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2",
        "double 0",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_decodes_all_two_bit_states_and_preserves_order(tmp_path):
    for suffix in (".bim", ".fam"):
        (tmp_path / f"chr1{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes()
        )
    _write_bed_calls(tmp_path / "chr1.bed", _CHR1_CALLS)

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}', {{'1'}}); "
        f"g = statgen.load_genotype('{tmp_path / 'chr1'}', ref); "
        "gi = g.fetch_genotypes_int8([2 1 2 5]); "
        "gf = g.fetch_genotypes([2 1 2 5]); "
        "gf(isnan(gf)) = -9; "
        "fprintf('%d ', gi); fprintf('\\n'); "
        "fprintf('%.0f ', gf); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "0 1 2 -1 2 -1 1 0 0 1 2 -1 -1 2 0 1",
        "0 1 2 -9 2 -9 1 0 0 1 2 -9 -9 2 0 1",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_chrx_subset_expands_to_panel_sample_axis(tmp_path):
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

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{tmp_path / '@'}', ref); "
        "gi = g.fetch_genotypes_int8(6); "
        "fprintf('%d ', gi); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 -1 0 -1"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_fetch_bed_path_overrides_and_validation_failures(tmp_path):
    bad_magic = tmp_path / "bad_magic.bed"
    bad_magic.write_bytes(b"\x00\x00\x01" + b"\x00" * 5)
    sample_major = tmp_path / "sample_major.bed"
    sample_major.write_bytes(b"\x6c\x1b\x00" + b"\x00" * 5)
    trailing = tmp_path / "trailing.bed"
    trailing.write_bytes(b"\x6c\x1b\x01" + b"\x00" * 6)

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"gs = statgen.load_genotype('{G_SHARDED}', ref); "
        f"gn = statgen.load_genotype('{G_NONSHARDED}', ref); "
        f"ok = isequal(gs.fetch_genotypes_int8([1 6], '{FIXTURES_DIR / 'genotype/sharded/@.bed'}'), int8(2 * ones(4, 2))); "
        f"ok = ok && isequal(gn.fetch_genotypes_int8([1 6], '{FIXTURES_DIR / 'genotype/nonsharded/all.bed'}'), int8(2 * ones(4, 2))); "
        f"e1 = 0; try; gn.fetch_genotypes_int8(1, '{tmp_path / '@.bed'}'); catch; e1 = 1; end; "
        f"e2 = 0; try; gs.fetch_genotypes_int8([1 6], '{tmp_path / 'flat.bed'}'); catch; e2 = 1; end; "
        f"ref1 = statgen.load_reference('{REF_SHARDED}', {{'1'}}); "
        f"g1 = statgen.load_genotype('{G_SHARDED}', ref1); "
        f"e3 = 0; try; g1.fetch_genotypes_int8(1, '{bad_magic}'); catch; e3 = 1; end; "
        f"e4 = 0; try; g1.fetch_genotypes_int8(1, '{sample_major}'); catch; e4 = 1; end; "
        f"e5 = 0; try; g1.fetch_genotypes_int8(1, '{trailing}'); catch; e5 = 1; end; "
        "fprintf('%d %d %d %d %d %d\\n', ok, e1, e2, e3, e4, e5);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1 1 1 1 1"]
