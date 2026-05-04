import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from statgen.ld import _write_ld_npz_distribution, load_ld, validate_ld_distribution
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, run_octave, skipif_no_octave


SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"


def _synthetic_ld_specs(reference):
    return [
        {
            "reference_shard": reference.shards[0],
            "a1freq": pd.DataFrame({"a1freq": [0.10, 0.20, 0.30, 0.40, 0.50]}),
            "ld_pairs": pd.DataFrame(
                {
                    "idx1": [0, 2],
                    "idx2": [1, 3],
                    "r": [0.25, -0.50],
                }
            ),
            "build_metadata": {
                "build_tool": "pytest",
                "build_command": "synthetic phase5a",
                "plink_version": None,
                "ld_window_kb": 10000,
                "ld_r2_threshold": 0.05,
                "num_sample": 4,
            },
        },
        {
            "reference_shard": reference.shards[1],
            "sex": "female",
            "a1freq": {"a1freq": [0.22, 0.33, 0.44]},
            "ld_pairs": {"idx1": [0], "idx2": [2], "r": [0.70]},
            "build_metadata": {
                "build_tool": "pytest",
                "build_command": "synthetic phase5a chrX",
                "num_sample": 4,
            },
        },
    ]


def _synthetic_ld_specs_all_chrx(reference):
    specs = _synthetic_ld_specs(reference)
    specs.extend(
        [
            {
                "reference_shard": reference.shards[1],
                "sex": "male",
                "a1freq": [0.11, 0.22, 0.33],
                "ld_pairs": np.array([[0, 1, 0.55], [1, 2, -0.45]], dtype=np.float64),
            },
            {
                "reference_shard": reference.shards[1],
                "sex": "combined",
                "a1freq": [0.15, 0.25, 0.35],
                "ld_pairs": np.array([[0, 2, 0.60]], dtype=np.float64),
                "extra_metadata": {
                    "chrX_combined_rationale": "synthetic phase5a coverage",
                },
            },
        ]
    )
    return specs


def _metadata(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as z:
        return json.loads(bytes(z["metadata"]).decode("utf-8"))


def test_importing_ld_does_not_eagerly_import_pandas():
    code = "import statgen.ld, sys; print('pandas' in sys.modules)"
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"


def test_write_ld_npz_distribution_from_synthetic_tables_validates_and_loads(tmp_path):
    reference = load_reference(SHARDED_REF)
    root = tmp_path / "ld_python"

    manifest = _write_ld_npz_distribution(root, _synthetic_ld_specs(reference))

    assert manifest["runtime_format"] == "python_npz_csc32"
    assert [entry["file"] for entry in manifest["shards"]] == [
        "ld_chr1.npz",
        "ld_chrX_female.npz",
    ]
    report = validate_ld_distribution(root, check_payload_structure=True)
    assert report["ok"] is True

    meta = _metadata(root / "ld_chr1.npz")
    assert meta["format"] == "statgen_ld_npz_csc32"
    assert meta["build_tool"] == "pytest"
    assert meta["num_sample"] == 4

    ld = load_ld(root, reference)
    chr1 = ld.shard_groups[0][0]
    chrx = ld.shard_groups[1][0]
    np.testing.assert_allclose(chr1.a1freq, [0.10, 0.20, 0.30, 0.40, 0.50])
    assert chr1.ld_r[0, 1] == pytest.approx(0.25)
    assert chr1.ld_r[2, 3] == pytest.approx(-0.50)
    assert chr1.ld_r.diagonal().tolist() == [1.0] * 5
    assert chrx.sex == "female"
    assert chrx.ld_r[0, 2] == pytest.approx(0.70)


def test_write_ld_npz_distribution_with_all_chrx_sexes_metadata_and_numpy_pairs(tmp_path):
    reference = load_reference(SHARDED_REF)
    root = tmp_path / "ld_python_all_x"

    manifest = _write_ld_npz_distribution(root, _synthetic_ld_specs_all_chrx(reference))

    assert [entry["file"] for entry in manifest["shards"]] == [
        "ld_chr1.npz",
        "ld_chrX_female.npz",
        "ld_chrX_male.npz",
        "ld_chrX_combined.npz",
    ]

    male_meta = _metadata(root / "ld_chrX_male.npz")
    assert male_meta["plink_version"] is None
    assert male_meta["ld_window_kb"] is None
    assert male_meta["ld_r2_threshold"] is None

    combined_meta = _metadata(root / "ld_chrX_combined.npz")
    assert combined_meta["chrX_combined_rationale"] == "synthetic phase5a coverage"

    ld = load_ld(root, reference)
    by_sex = {shard.sex: shard for shard in ld.shard_groups[1]}
    assert set(by_sex) == {"female", "male", "combined"}
    assert by_sex["male"].ld_r[0, 1] == pytest.approx(0.55)
    assert by_sex["male"].ld_r[1, 2] == pytest.approx(-0.45)
    assert by_sex["combined"].ld_r[0, 2] == pytest.approx(0.60)


def test_ld_npz_writer_rejects_malformed_synthetic_tables(tmp_path):
    reference = load_reference(SHARDED_REF)
    base = _synthetic_ld_specs(reference)[0]

    bad = dict(base)
    bad["ld_pairs"] = pd.DataFrame({"idx1": [0, 1], "idx2": [1, 0], "r": [0.2, 0.2]})
    with pytest.raises(ValueError, match="duplicate unordered pairs"):
        _write_ld_npz_distribution(tmp_path / "dup", [bad])

    bad = dict(base)
    bad["a1freq"] = {"a1freq": [0.1, 0.2]}
    with pytest.raises(ValueError, match="a1freq length mismatch"):
        _write_ld_npz_distribution(tmp_path / "freq", [bad])

    bad = dict(base)
    bad["ld_pairs"] = pd.DataFrame({"idx1": [0], "idx2": [99], "r": [0.2]})
    with pytest.raises(ValueError, match="indices out of bounds"):
        _write_ld_npz_distribution(tmp_path / "bounds", [bad])

    bad = dict(base)
    bad["ld_pairs"] = pd.DataFrame({"i": [0], "j": [1], "r": [0.2]})
    with pytest.raises(ValueError, match="idx1, idx2, r"):
        _write_ld_npz_distribution(tmp_path / "undocumented_alias", [bad])

    bad = dict(base)
    bad.pop("a1freq")
    bad["freq"] = {"a1freq": [0.1, 0.2, 0.3, 0.4, 0.5]}
    with pytest.raises(ValueError, match="unknown fields: freq"):
        _write_ld_npz_distribution(tmp_path / "freq_alias", [bad])

    bad = dict(base)
    bad.pop("ld_pairs")
    bad["ld"] = pd.DataFrame({"idx1": [0], "idx2": [1], "r": [0.2]})
    with pytest.raises(ValueError, match="unknown fields: ld"):
        _write_ld_npz_distribution(tmp_path / "ld_alias", [bad])


def test_ld_npz_writer_enforces_manifest_invariants_when_validation_disabled(tmp_path):
    reference = load_reference(SHARDED_REF)
    base = _synthetic_ld_specs(reference)[0]

    with pytest.raises(ValueError, match="at least one shard"):
        _write_ld_npz_distribution(tmp_path / "empty", [], validate=False)

    duplicate_key = [dict(base), dict(base)]
    with pytest.raises(ValueError, match="duplicate LD shard"):
        _write_ld_npz_distribution(tmp_path / "duplicate_key", duplicate_key, validate=False)

    other = dict(_synthetic_ld_specs(reference)[1])
    duplicate_file = [dict(base, file="same.npz"), dict(other, file="same.npz")]
    with pytest.raises(ValueError, match="duplicate LD shard output file"):
        _write_ld_npz_distribution(tmp_path / "duplicate_file", duplicate_file, validate=False)

    bad = dict(base)
    bad["ld_pairs"] = pd.DataFrame({"idx1": [0], "idx2": [1], "r": [1.5]})
    with pytest.raises(ValueError, match="between -1 and 1"):
        _write_ld_npz_distribution(tmp_path / "bad_r", [bad], validate=False)


@pytest.mark.octave
@skipif_no_octave
def test_octave_npz_to_mat_conversion_validates_and_matches_generated_npz(tmp_path):
    reference = load_reference(SHARDED_REF)
    py_root = tmp_path / "ld_python"
    mat_root = tmp_path / "ld_matlab"
    _write_ld_npz_distribution(py_root, _synthetic_ld_specs(reference))
    py_chr1_meta = _metadata(py_root / "ld_chr1.npz")
    py_x_meta = _metadata(py_root / "ld_chrX_female.npz")

    script = (
        f"ref = statgen.load_reference('{SHARDED_REF}'); "
        f"statgen.convert_ld_npz_to_mat('{py_root}', '{mat_root}'); "
        f"report = statgen.validate_ld_distribution('{mat_root}', true); "
        f"ld = statgen.load_ld('{mat_root}', ref); "
        f"chr1_payload = load('{mat_root / 'ld_chr1.mat'}'); "
        f"x_payload = load('{mat_root / 'ld_chrX_female.mat'}'); "
        "fprintf('%d\\n', report.ok); "
        "fprintf('%d\\n', numel(report.warnings)); "
        "fprintf('%.2f,%.2f\\n', full(ld.shard_groups{1}{1}.ld_r(1,2)), full(ld.shard_groups{1}{1}.ld_r(3,4))); "
        "fprintf('%.2f\\n', ld.shard_groups{2}{1}.a1freq(2)); "
        "fprintf('%s\\n', ld.shard_groups{2}{1}.sex); "
        "fprintf('%d,%d,%s\\n', chr1_payload.metadata.num_snp, chr1_payload.metadata.nnz, chr1_payload.metadata.reference_checksum); "
        "fprintf('%d,%d,%s\\n', x_payload.metadata.num_snp, x_payload.metadata.nnz, x_payload.metadata.reference_checksum);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
    assert lines[1] == "1"
    assert lines[2] == "0.25,-0.50"
    assert lines[3] == "0.33"
    assert lines[4] == "female"
    assert lines[5] == f"{py_chr1_meta['num_snp']},{py_chr1_meta['nnz']},{py_chr1_meta['reference_checksum']}"
    assert lines[6] == f"{py_x_meta['num_snp']},{py_x_meta['nnz']},{py_x_meta['reference_checksum']}"
