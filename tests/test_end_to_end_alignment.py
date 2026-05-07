import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from statgen.genotype import load_genotype
from statgen.ld import load_ld, validate_ld_distribution
from tests.conftest import run_octave, skipif_no_octave


BUILD_SCRIPT = Path(__file__).parents[1] / "script/statgen_build_ld.py"
MANIFEST_SCRIPT = Path(__file__).parents[1] / "script/statgen_create_ld_manifest.py"


def _preflight_variant_rows(num_snp: int = 31):
    allele_pairs = [
        ("G", "A"),
        ("T", "C"),
        ("C", "G"),
        ("A", "T"),
        ("G", "C"),
        ("T", "A"),
    ]
    rows = []
    for i in range(num_snp):
        ref, alt = allele_pairs[i % len(allele_pairs)]
        rows.append(
            {
                "chrom": "1",
                "pos": 1000 + i * 100,
                "id": f"e2e_preflight_chr1_{i + 1:03d}_{alt}_{ref}",
                "ref": ref,
                "alt": alt,
            }
        )
    return rows


def _preflight_calls(num_snp: int = 31, num_sample: int = 29) -> np.ndarray:
    calls = np.fromfunction(
        lambda i, j: (i * 7 + j * 5 + (j % 4)) % 3,
        (num_snp, num_sample),
        dtype=int,
    ).astype(np.int8)

    sentinel = np.array(
        [0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 1, 2, 0, 1, 2,
         0, 1, 2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1],
        dtype=np.int8,
    )
    calls[0] = sentinel
    calls[1] = sentinel.copy()
    calls[1, [5, 17, 24]] = np.array([2, 0, 0], dtype=np.int8)
    calls[2] = 2 - sentinel
    calls[2, [8, 16, 25]] = np.array([0, 2, 1], dtype=np.int8)

    missing_positions = {
        1: [3, 21],
        2: [4, 20],
        5: [0, 11, 22],
        12: [2, 14],
        23: [6, 18, 27],
    }
    for snp_idx, sample_idx in missing_positions.items():
        calls[snp_idx, sample_idx] = -1

    for i, row in enumerate(calls):
        observed = row[row >= 0]
        if observed.size == 0 or np.unique(observed).size == 1:
            raise AssertionError(f"preflight generated monomorphic SNP row {i}")
    return calls


def _write_preflight_vcf(path: Path, variants, calls: np.ndarray) -> list[str]:
    samples = [f"S{i + 1:02d}" for i in range(calls.shape[1])]

    def gt(value: int) -> str:
        if value == -1:
            return "./."
        if value == 0:
            return "0/0"
        if value == 1:
            return "0/1"
        if value == 2:
            return "1/1"
        raise ValueError(value)

    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for row, call_row in zip(variants, calls):
        lines.append(
            "\t".join(
                [
                    row["chrom"],
                    str(row["pos"]),
                    row["id"],
                    row["ref"],
                    row["alt"],
                    ".",
                    "PASS",
                    ".",
                    "GT",
                    *[gt(int(x)) for x in call_row],
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return samples


def _write_preflight_psam(path: Path, samples: list[str]) -> None:
    lines = ["#IID\tSEX"]
    for i, sample in enumerate(samples):
        if i < 14:
            sex = "1"
        elif i < 28:
            sex = "2"
        else:
            sex = "0"
        lines.append(f"{sample}\t{sex}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_checked(cmd: list[str]) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr + result.stdout
    return result


def _require_plink2() -> str:
    plink2 = shutil.which("plink2")
    if plink2 is None:
        pytest.skip("plink2 not installed")
    return plink2


def _create_manifest(out: Path) -> None:
    _run_checked([sys.executable, str(MANIFEST_SCRIPT), "--ld", str(out)])


def _read_bim(path: Path):
    rows = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        chrom, snp, cm, bp, a1, a2 = raw.split()
        rows.append({"chrom": chrom, "id": snp, "bp": int(bp), "a1": a1, "a2": a2})
    return rows


def _assert_bim_preserves_vcf_alleles(bim_path: Path, variants) -> None:
    bim = _read_bim(bim_path)
    assert len(bim) == len(variants)
    for got, expected in zip(bim, variants):
        assert got["chrom"] == expected["chrom"]
        assert got["id"] == expected["id"]
        assert got["bp"] == expected["pos"]
        assert got["a1"] == expected["alt"]
        assert got["a2"] == expected["ref"]


def _allele_frequency_from_fetch(genotypes: np.ndarray) -> np.ndarray:
    observed = np.isfinite(genotypes)
    denom = 2.0 * observed.sum(axis=0)
    return np.nansum(genotypes, axis=0) / denom


def _pearson_pair(genotypes: np.ndarray, i: int, j: int) -> float:
    x = genotypes[:, i]
    y = genotypes[:, j]
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 2:
        raise AssertionError("not enough observed calls for Pearson r")
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def _build_preflight_artifacts(tmp_path: Path) -> tuple[Path, Path]:
    plink2 = _require_plink2()
    variants = _preflight_variant_rows()
    calls = _preflight_calls()
    vcf = tmp_path / "preflight.vcf"
    psam = tmp_path / "preflight.psam"
    prefix = tmp_path / "preflight"
    samples = _write_preflight_vcf(vcf, variants, calls)
    _write_preflight_psam(psam, samples)

    _run_checked(
        [
            plink2,
            "--vcf",
            str(vcf),
            "--psam",
            str(psam),
            "--make-bed",
            "--keep-allele-order",
            "--out",
            str(prefix),
        ]
    )
    _assert_bim_preserves_vcf_alleles(prefix.with_suffix(".bim"), variants)

    ld_root = tmp_path / "ld"
    _run_checked(
        [
            sys.executable,
            str(BUILD_SCRIPT),
            "--plink2",
            plink2,
            "--bfile",
            str(prefix),
            "--out",
            str(ld_root),
            "--shard",
            "1",
        ]
    )
    _create_manifest(ld_root)
    assert validate_ld_distribution(ld_root, check_payload_structure=True)["ok"] is True
    return prefix, ld_root


def _python_preflight_outputs(prefix: Path, ld_root: Path) -> dict[str, np.ndarray]:
    ld = load_ld(ld_root)
    genotype = load_genotype(prefix, ld.reference)
    assert genotype.is_present.tolist() == [True] * ld.reference.num_snp

    fetched = genotype.fetch_genotypes(np.arange(genotype.num_snp))
    assert fetched.shape == (29, 31)
    ld_r = ld.shard_groups[0][0].ld_r
    vec = np.linspace(-1.5, 2.5, genotype.num_snp, dtype=np.float64)
    mat = np.column_stack([vec, vec[::-1]])
    return {
        "is_present": genotype.is_present.astype(np.float64),
        "ld_a1freq": ld.a1freq().astype(np.float64),
        "fetch_a1freq": _allele_frequency_from_fetch(fetched).astype(np.float64),
        "ld_r_selected": np.array([float(ld_r[0, 1]), float(ld_r[0, 2])], dtype=np.float64),
        "direct_r_selected": np.array(
            [_pearson_pair(fetched, 0, 1), _pearson_pair(fetched, 0, 2)],
            dtype=np.float64,
        ),
        "multiply_vec": ld.multiply_r2(vec).astype(np.float64),
        "explicit_vec": ((ld_r.toarray() ** 2) @ vec).astype(np.float64),
        "multiply_mat": ld.multiply_r2(mat).astype(np.float64).reshape(-1, order="F"),
        "explicit_mat": (((ld_r.toarray() ** 2) @ mat).astype(np.float64)).reshape(-1, order="F"),
    }


def _assert_preflight_outputs(outputs: dict[str, np.ndarray]) -> None:
    assert outputs["is_present"].tolist() == [1.0] * 31
    np.testing.assert_allclose(outputs["ld_a1freq"], outputs["fetch_a1freq"], atol=1e-5)
    for expected_r in outputs["direct_r_selected"]:
        assert abs(expected_r) > np.sqrt(0.05)
    np.testing.assert_allclose(outputs["ld_r_selected"], outputs["direct_r_selected"], atol=1e-3)
    np.testing.assert_allclose(
        outputs["multiply_vec"],
        outputs["explicit_vec"],
        atol=1e-5,
    )
    np.testing.assert_allclose(
        outputs["multiply_mat"],
        outputs["explicit_mat"],
        atol=1e-5,
    )


def _matlab_quote(path: Path) -> str:
    return str(path).replace("'", "''")


def _read_labeled_numeric_output(path: Path) -> dict[str, np.ndarray]:
    out = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.strip().split()
        if not fields:
            continue
        out[fields[0]] = np.array([float(x) for x in fields[1:]], dtype=np.float64)
    return out


def test_preflight_chr1_a1freq_and_ld_r_match_fetched_genotypes(tmp_path):
    prefix, ld_root = _build_preflight_artifacts(tmp_path)
    _assert_preflight_outputs(_python_preflight_outputs(prefix, ld_root))


@pytest.mark.octave
@skipif_no_octave
def test_preflight_chr1_octave_matches_python_after_npz_to_mat_conversion(tmp_path):
    prefix, ld_root = _build_preflight_artifacts(tmp_path)
    expected = _python_preflight_outputs(prefix, ld_root)
    mat_root = tmp_path / "ld_matlab"
    out_path = tmp_path / "octave_preflight.tsv"

    script = (
        "warning('off', 'statgen:ld:v5mat'); "
        f"statgen.convert_ld_npz_to_mat('{_matlab_quote(ld_root)}', '{_matlab_quote(mat_root)}', '1', 'format', 'v5'); "
        f"report = statgen.create_ld_mat_manifest('{_matlab_quote(ld_root)}', '{_matlab_quote(mat_root)}', {{'1'}}); "
        f"validate_report = statgen.validate_ld_distribution('{_matlab_quote(mat_root)}', true); "
        f"ld = statgen.load_ld('{_matlab_quote(mat_root)}'); "
        f"g = statgen.load_genotype('{_matlab_quote(prefix)}', ld.reference); "
        "G = g.fetch_genotypes(1:g.num_snp); "
        "obs = isfinite(G); G0 = G; G0(~obs) = 0; "
        "freq = sum(G0, 1) ./ (2 * sum(obs, 1)); "
        "r01_keep = isfinite(G(:,1)) & isfinite(G(:,2)); "
        "r02_keep = isfinite(G(:,1)) & isfinite(G(:,3)); "
        "r01 = corrcoef(G(r01_keep,1), G(r01_keep,2)); "
        "r02 = corrcoef(G(r02_keep,1), G(r02_keep,3)); "
        "ld_r = ld.shard_groups{1}{1}.ld_r; "
        "v = linspace(-1.5, 2.5, g.num_snp)'; "
        "M = [v, flipud(v)]; "
        "mv = ld.multiply_r2(v); "
        "mm = ld.multiply_r2(M); "
        "explicit_v = full((ld_r .^ 2) * v); "
        "explicit_m = full((ld_r .^ 2) * M); "
        f"fid = fopen('{_matlab_quote(out_path)}', 'w'); "
        "fprintf(fid, 'manifest_count %.0f\\n', numel(report.shards)); "
        "fprintf(fid, 'validate_ok %.0f\\n', validate_report.ok); "
        "fprintf(fid, 'is_present '); fprintf(fid, '%.0f ', g.is_present); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_a1freq '); fprintf(fid, '%.10g ', ld.a1freq()); fprintf(fid, '\\n'); "
        "fprintf(fid, 'fetch_a1freq '); fprintf(fid, '%.10g ', freq); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_r_selected %.10g %.10g\\n', full(ld_r(1,2)), full(ld_r(1,3))); "
        "fprintf(fid, 'direct_r_selected %.10g %.10g\\n', r01(1,2), r02(1,2)); "
        "fprintf(fid, 'multiply_vec '); fprintf(fid, '%.10g ', mv); fprintf(fid, '\\n'); "
        "fprintf(fid, 'explicit_vec '); fprintf(fid, '%.10g ', explicit_v); fprintf(fid, '\\n'); "
        "fprintf(fid, 'multiply_mat '); fprintf(fid, '%.10g ', mm(:)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'explicit_mat '); fprintf(fid, '%.10g ', explicit_m(:)); fprintf(fid, '\\n'); "
        "fclose(fid);"
    )
    result = run_octave(script, timeout=60)
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
