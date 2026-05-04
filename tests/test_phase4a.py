import json
import hashlib
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from statgen.ld import load_ld, validate_ld_distribution
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, matlab_data_lines, run_octave, skipif_no_octave


SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
LD_PY = FIXTURES_DIR / "ld/python"
LD_MAT = FIXTURES_DIR / "ld/matlab"


def _read_manifest(root: Path) -> dict:
    return json.loads((root / "ld_manifest.json").read_text())


def _write_manifest(root: Path, manifest: dict) -> None:
    (root / "ld_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _copy_ld_distribution(tmp_path: Path) -> Path:
    out = tmp_path / "ld"
    shutil.copytree(LD_PY, out)
    return out


def _load_npz_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def _write_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    np.savez_compressed(path, **arrays)


def _replace_npz_metadata(path: Path, updates: dict) -> None:
    arrays = _load_npz_arrays(path)
    meta = json.loads(bytes(arrays["metadata"]).decode())
    meta.update(updates)
    arrays["metadata"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    _write_npz(path, arrays)


def _octave_script(expr: str) -> str:
    return f"addpath('{MATLAB_DIR}'); fixture_dir = '{FIXTURES_DIR}'; " + expr


def test_load_ld_panel_matches_fixture_sparse_payloads():
    reference = load_reference(SHARDED_REF)
    ld = load_ld(LD_PY, reference)

    assert ld.default_chrX_sex == "female"
    assert [s.label for s in ld.shards] == ["1", "X"]
    assert [len(g) for g in ld.shard_groups] == [1, 3]
    assert reference.is_object_compatible(ld) is True

    chr1 = ld.shard_groups[0][0]
    assert chr1.sex is None
    assert sparse.isspmatrix_csc(chr1.ld_r)
    np.testing.assert_allclose(
        chr1.a1freq,
        np.array([0.30, 0.40, 0.20, 0.35, 0.15], dtype=np.float32),
    )
    np.testing.assert_allclose(chr1.ld_r.toarray(), chr1.ld_r.toarray().T)
    assert chr1.ld_r[0, 1] == pytest.approx(0.9)
    assert chr1.ld_r[1, 3] == pytest.approx(-0.5)

    x_by_sex = {s.sex: s for s in ld.shard_groups[1]}
    assert set(x_by_sex) == {"female", "male", "combined"}
    assert x_by_sex["female"].ld_r[0, 1] == pytest.approx(0.65)
    assert x_by_sex["male"].ld_r[0, 1] == pytest.approx(0.55)
    assert x_by_sex["combined"].ld_r[1, 2] == pytest.approx(-0.40)


def test_load_ld_single_shard_with_single_chromosome_reference():
    reference = load_reference(SHARDED_REF, shards=["1"])
    ld = load_ld(LD_PY / "ld_chr1.npz", reference)

    assert [s.label for s in ld.shards] == ["1"]
    assert ld.shards[0].num_snp == 5
    assert reference.is_object_compatible(ld) is True


def test_load_ld_single_shard_reference_shape_errors_are_specific():
    reference = load_reference(SHARDED_REF)
    with pytest.raises(ValueError, match="single-shard reference; reference has 2 shards"):
        load_ld(LD_PY / "ld_chr1.npz", reference)

    x_reference = load_reference(SHARDED_REF, shards=["X"])
    with pytest.raises(ValueError, match="chromosome does not match"):
        load_ld(LD_PY / "ld_chr1.npz", x_reference)


def test_load_ld_default_chrx_sex_validation():
    reference = load_reference(SHARDED_REF)
    ld = load_ld(LD_PY, reference, default_chrX_sex="male")
    assert ld.default_chrX_sex == "male"
    assert ld.shards[1].sex == "male"

    with pytest.raises(ValueError, match="chrX sex"):
        load_ld(LD_PY, reference, default_chrX_sex="unknown")


def test_validate_ld_distribution_checks_manifest_md5_and_payload_structure(tmp_path):
    report = validate_ld_distribution(LD_PY, check_payload_structure=True)
    assert report["ok"] is True

    root = _copy_ld_distribution(tmp_path)
    arrays = _load_npz_arrays(root / "ld_chr1.npz")
    arrays["a1freq"] = arrays["a1freq"][:-1]
    _write_npz(root / "ld_chr1.npz", arrays)

    with pytest.raises(ValueError, match="file_md5"):
        validate_ld_distribution(root)

    manifest = _read_manifest(root)
    manifest["shards"][0]["file_md5"] = _md5_file(root / "ld_chr1.npz")
    _write_manifest(root, manifest)
    with pytest.raises(ValueError, match="a1freq length"):
        validate_ld_distribution(root, check_payload_structure=True)


def test_load_ld_missing_manifest_shard_file_fails(tmp_path):
    reference = load_reference(SHARDED_REF)
    root = _copy_ld_distribution(tmp_path)
    (root / "ld_chr1.npz").unlink()

    with pytest.raises(FileNotFoundError, match="LD file not found"):
        load_ld(root, reference)


def test_invalid_ld_manifest_fails_clearly(tmp_path):
    reference = load_reference(SHARDED_REF)

    bad_object = _copy_ld_distribution(tmp_path / "bad_object")
    manifest = _read_manifest(bad_object)
    manifest["object_type"] = "not_ld_panel_manifest"
    _write_manifest(bad_object, manifest)
    with pytest.raises(ValueError, match="object_type"):
        load_ld(bad_object, reference)

    bad_schema = _copy_ld_distribution(tmp_path / "bad_schema")
    manifest = _read_manifest(bad_schema)
    manifest["schema_version"] = "bad"
    _write_manifest(bad_schema, manifest)
    with pytest.raises(ValueError, match="schema_version"):
        validate_ld_distribution(bad_schema)

    malformed = _copy_ld_distribution(tmp_path / "malformed")
    (malformed / "ld_manifest.json").write_text("{not json\n")
    with pytest.raises(json.JSONDecodeError):
        load_ld(malformed, reference)


def test_load_ld_rejects_manifest_metadata_and_reference_mismatches(tmp_path):
    reference = load_reference(SHARDED_REF)
    root = _copy_ld_distribution(tmp_path)

    manifest = _read_manifest(root)
    manifest["shards"][0]["num_snp"] = 99
    _write_manifest(root, manifest)
    with pytest.raises(ValueError, match="manifest/per-file metadata mismatch"):
        load_ld(root, reference)

    root = _copy_ld_distribution(tmp_path / "checksum")
    _replace_npz_metadata(root / "ld_chr1.npz", {"reference_checksum": "deadbeef" * 4})
    manifest = _read_manifest(root)
    manifest["shards"][0]["reference_checksum"] = "deadbeef" * 4
    _write_manifest(root, manifest)
    with pytest.raises(ValueError, match="reference_checksum"):
        load_ld(root, reference)


def test_load_ld_rejects_unknown_chrx_sex_label(tmp_path):
    reference = load_reference(SHARDED_REF)
    root = _copy_ld_distribution(tmp_path)
    _replace_npz_metadata(root / "ld_chrX_female.npz", {"sex": "unknown"})
    manifest = _read_manifest(root)
    manifest["shards"][1]["sex"] = "unknown"
    _write_manifest(root, manifest)

    with pytest.raises(ValueError, match="chrX sex"):
        load_ld(root, reference)


def test_payload_structure_detects_malformed_sparse_indices(tmp_path):
    root = _copy_ld_distribution(tmp_path)
    arrays = _load_npz_arrays(root / "ld_chr1.npz")
    arrays["indices"] = arrays["indices"].copy()
    arrays["indices"][0] = 999
    _write_npz(root / "ld_chr1.npz", arrays)
    manifest = _read_manifest(root)
    manifest["shards"][0]["file_md5"] = _md5_file(root / "ld_chr1.npz")
    _write_manifest(root, manifest)

    with pytest.raises(ValueError, match="row indices out of bounds"):
        validate_ld_distribution(root, check_payload_structure=True)


def test_python_validator_rejects_matlab_ld_distribution():
    with pytest.raises(ValueError, match="expected runtime_format 'python_npz_csc32'"):
        validate_ld_distribution(LD_MAT)
    with pytest.raises(ValueError, match="use MATLAB/Octave"):
        validate_ld_distribution(LD_MAT / "ld_chr1.mat")


def test_ld_metadata_rejects_bool_integer_fields(tmp_path):
    root = _copy_ld_distribution(tmp_path)
    _replace_npz_metadata(root / "ld_chr1.npz", {"num_snp": True})
    manifest = _read_manifest(root)
    manifest["shards"][0]["file_md5"] = _md5_file(root / "ld_chr1.npz")
    _write_manifest(root, manifest)

    with pytest.raises(ValueError, match="positive integer"):
        validate_ld_distribution(root)


def test_validate_ld_distribution_missing_manifest_shard_file_fails(tmp_path):
    root = _copy_ld_distribution(tmp_path)
    (root / "ld_chr1.npz").unlink()

    with pytest.raises(FileNotFoundError, match="LD file not found"):
        validate_ld_distribution(root)


def test_npz_archives_are_not_pickle_payloads():
    with zipfile.ZipFile(LD_PY / "ld_chr1.npz") as zf:
        assert all(name.endswith(".npy") for name in zf.namelist())
    report = validate_ld_distribution(LD_PY)
    assert report["ok"] is True


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_ld_mat_fixture_sparse_payloads():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "ld = statgen.load_ld([fixture_dir '/ld/matlab'], ref); "
        "fprintf('%s\\n', ld.default_chrX_sex); "
        "fprintf('%d,%d\\n', numel(ld.shard_groups{1}), numel(ld.shard_groups{2})); "
        "fprintf('%d\\n', ref.is_object_compatible(ld)); "
        "chr1 = ld.shard_groups{1}{1}; "
        "fprintf('%.2f,%.2f\\n', full(chr1.ld_r(1,2)), full(chr1.ld_r(2,4))); "
        "fprintf('%.2f\\n', chr1.a1freq(1)); "
        "x = ld.shard_groups{2}{1}; "
        "fprintf('%s,%.2f\\n', x.sex, full(x.ld_r(1,2)));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "female"
    assert lines[1] == "1,3"
    assert lines[2] == "1"
    assert lines[3] == "0.90,-0.50"
    assert lines[4] == "0.30"
    assert lines[5] == "female,0.65"


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_ld_single_mat_shard():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim'], {'1'}); "
        "ld = statgen.load_ld([fixture_dir '/ld/matlab/ld_chr1.mat'], ref); "
        "fprintf('%d\\n', numel(ld.shard_groups)); "
        "fprintf('%s\\n', ld.shards{1}.label); "
        "fprintf('%d\\n', ref.is_object_compatible(ld)); "
        "fprintf('%.2f\\n', full(ld.shards{1}.ld_r(1,2)));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
    assert lines[1] == "1"
    assert lines[2] == "1"
    assert lines[3] == "0.90"


@pytest.mark.octave
@skipif_no_octave
def test_octave_validate_ld_distribution_and_bad_chrx_sex(tmp_path):
    bad_root = tmp_path / "ld_bad"
    shutil.copytree(LD_MAT, bad_root)
    manifest = _read_manifest(bad_root)
    manifest["shards"][1]["sex"] = "unknown"
    _write_manifest(bad_root, manifest)

    script = _octave_script(
        "warning('off', 'statgen:ld:v5mat'); "
        "report = statgen.validate_ld_distribution([fixture_dir '/ld/matlab'], true); "
        "fprintf('%d\\n', report.ok); "
        f"try; ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); statgen.load_ld('{bad_root}', ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = matlab_data_lines(result.stdout)
    assert lines[0] == "1"
    assert lines[-1] == "FAIL"
