import shutil

import pytest

from tests.conftest import FIXTURES_DIR, matlab_data_lines, run_octave, skipif_no_octave


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
REF_NONSHARDED = FIXTURES_DIR / "reference/nonsharded/all.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


def _copy_sharded_genotype(dst):
    for label in ("1", "X"):
        for suffix in (".bim", ".fam", ".bed", ".ploidy"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            if src.exists():
                (dst / f"{label}{suffix}").write_bytes(src.read_bytes())


def _write_bed(path, num_snp, num_sample):
    bytes_per_snp = (num_sample + 3) // 4
    path.write_bytes(b"\x6c\x1b\x01" + b"\x00" * (num_snp * bytes_per_snp))


def _write_bim(path):
    path.write_text(
        "1\trs1001\t0\t100\tA\tG\n"
        "1\trs1002\t0\t200\tC\tT\n"
        "1\trs1003\t0\t300\tA\tC\n"
        "1\trs1004\t0\t400\tG\tA\n"
        "1\trs1005\t0\t500\tT\tC\n"
    )


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_sharded_metadata_loads():
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        "fprintf('%s\\n', g.source_layout); "
        "fprintf('%d %d %d\\n', numel(g.shards), g.num_snp, g.num_sample); "
        "for i = 1:numel(g.shards); fprintf('%s ', g.shards{i}.label); end; fprintf('\\n'); "
        "fprintf('%.0f ', g.source_row0); fprintf('\\n'); "
        "fprintf('%.0f ', g.ploidy_male); fprintf('\\n'); "
        "fprintf('%.0f ', g.ploidy_female); fprintf('\\n'); "
        "fprintf('%s ', g.fid{:}); fprintf('\\n'); "
        "fprintf('%.0f ', g.sex); fprintf('\\n'); "
        "fprintf('%d ', g.is_male); fprintf('\\n'); "
        "fprintf('%d ', g.is_subject_present('X')); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "sharded",
        "2 8 4",
        "1 X",
        "0 1 2 3 4 0 1 2",
        "2 2 2 2 2 1 1 1",
        "2 2 2 2 2 2 2 2",
        "FAM1 FAM1 FAM2 FAM2",
        "1 2 1 2",
        "1 0 1 0",
        "1 1 1 1",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_nonsharded_metadata_loads():
    script = (
        f"ref = statgen.load_reference('{REF_NONSHARDED}'); "
        f"g = statgen.load_genotype('{G_NONSHARDED}', ref); "
        "fprintf('%s\\n', g.source_layout); "
        "fprintf('%d %d %d\\n', numel(g.shards), g.num_snp, g.num_sample); "
        "fprintf('%s\\n', g.shards{1}.bed_path); "
        "fprintf('%s\\n', g.shards{2}.bed_path); "
        "fprintf('%.0f ', g.source_row0); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "non_sharded"
    assert lines[1] == "2 8 4"
    assert lines[2] == str(FIXTURES_DIR / "genotype/nonsharded/all.bed")
    assert lines[3] == str(FIXTURES_DIR / "genotype/nonsharded/all.bed")
    assert lines[4] == "0 1 2 3 4 5 6 7"


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_chrx_subset_fam_maps_to_panel_axis(tmp_path):
    _copy_sharded_genotype(tmp_path)
    (tmp_path / "X.fam").write_text(
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t1\t-9\n"
    )

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{tmp_path / '@'}', ref); "
        "fprintf('%d ', g.is_subject_present('X')); fprintf('\\n'); "
        "fprintf('%.0f ', g.shards{2}.source_subject_row0); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == ["1 0 1 0", "1 -1 0 -1"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_cache_roundtrip_and_subset(tmp_path):
    cache = tmp_path / "genotype_cache.mat"
    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"g = statgen.load_genotype('{G_SHARDED}', ref); "
        f"statgen.save_genotype_cache(g, '{cache}'); "
        f"h = statgen.load_genotype_cache('{cache}', {{'X'}}); "
        "fprintf('%s %d %d\\n', h.source_layout, h.num_snp, h.num_sample); "
        "fprintf('%s\\n', h.shards{1}.label); "
        "fprintf('%.0f ', h.source_row0); fprintf('\\n'); "
        "fprintf('%.0f ', h.ploidy_male); fprintf('\\n'); "
        "fprintf('%s ', h.fid{:}); fprintf('\\n'); "
        "fprintf('%d ', h.is_subject_present('X')); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines == [
        "sharded 3 4",
        "X",
        "0 1 2",
        "1 1 1",
        "FAM1 FAM1 FAM2 FAM2",
        "1 1 1 1",
    ]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_missing_source_shard_fails(tmp_path):
    for suffix in (".bim", ".bed", ".fam"):
        shutil.copyfile(FIXTURES_DIR / f"genotype/sharded/1{suffix}", tmp_path / f"1{suffix}")

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}'); "
        f"ok = 0; try; statgen.load_genotype('{tmp_path / '@'}', ref); catch; ok = 1; end; "
        "fprintf('%d\\n', ok);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1"]


@pytest.mark.octave
@skipif_no_octave
def test_octave_genotype_fam_and_ploidy_validation_failures(tmp_path):
    _write_bim(tmp_path / "badsex.bim")
    (tmp_path / "badsex.fam").write_text("FAM1\tIND1\t0\t0\t3\t-9\n")
    _write_bed(tmp_path / "badsex.bed", num_snp=5, num_sample=1)

    _write_bim(tmp_path / "badploidy.bim")
    (tmp_path / "badploidy.fam").write_text("FAM1\tIND1\t0\t0\t1\t-9\n")
    _write_bed(tmp_path / "badploidy.bed", num_snp=5, num_sample=1)
    (tmp_path / "badploidy.ploidy").write_text("2\t2\n")

    script = (
        f"ref = statgen.load_reference('{REF_SHARDED}', {{'1'}}); "
        f"ok1 = 0; try; statgen.load_genotype('{tmp_path / 'badsex'}', ref); catch; ok1 = 1; end; "
        f"ok2 = 0; try; statgen.load_genotype('{tmp_path / 'badploidy'}', ref); catch; ok2 = 1; end; "
        "fprintf('%d %d\\n', ok1, ok2);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1"]
