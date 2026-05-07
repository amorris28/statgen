import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
REPO_ROOT = Path(__file__).parent.parent
MATLAB_DIR = REPO_ROOT / "matlab"

# Set STATGEN_MATLAB=1 to run octave-marked tests via native MATLAB instead.
_USE_MATLAB = os.environ.get("STATGEN_MATLAB", "0") == "1"
_ENGINE = "matlab" if _USE_MATLAB else "octave"
_ENGINE_AVAILABLE = shutil.which(_ENGINE) is not None
_MATLAB_ENGINE = None
_MATLAB_ENGINE_ERROR = None


@pytest.fixture(scope="session")
def octave_available():
    return _ENGINE_AVAILABLE


skipif_no_octave = pytest.mark.skipif(
    not _ENGINE_AVAILABLE,
    reason=f"{_ENGINE} not installed",
)

skipif_no_matlab = pytest.mark.skipif(
    shutil.which("matlab") is None,
    reason="MATLAB not installed",
)

skipif_matlab_engine = pytest.mark.skipif(
    _USE_MATLAB,
    reason="Octave-specific test",
)


def _get_matlab_engine():
    global _MATLAB_ENGINE, _MATLAB_ENGINE_ERROR
    if _MATLAB_ENGINE is not None or _MATLAB_ENGINE_ERROR is not None:
        return _MATLAB_ENGINE

    try:
        import matlab.engine
    except Exception as exc:  # pragma: no cover - environment dependent
        _MATLAB_ENGINE_ERROR = (
            "MATLAB engine for Python is required when STATGEN_MATLAB=1. "
            f"Import failed: {exc}"
        )
        return None

    try:
        eng = matlab.engine.start_matlab("-nosplash -nodesktop")
        eng.addpath(str(MATLAB_DIR), nargout=0)
        _MATLAB_ENGINE = eng
        return _MATLAB_ENGINE
    except Exception as exc:  # pragma: no cover - environment dependent
        _MATLAB_ENGINE_ERROR = f"Failed to start MATLAB engine: {exc}"
        return None


@pytest.fixture(scope="session", autouse=True)
def _close_matlab_engine_at_end():
    yield
    global _MATLAB_ENGINE
    if _MATLAB_ENGINE is not None:
        try:
            _MATLAB_ENGINE.quit()
        except Exception:
            pass
        _MATLAB_ENGINE = None


def run_octave(expr: str, timeout: int = 30) -> subprocess.CompletedProcess:
    """Run an expression via Octave or MATLAB (controlled by STATGEN_MATLAB=1)."""
    if _USE_MATLAB:
        eng = _get_matlab_engine()
        if eng is None:
            return subprocess.CompletedProcess(
                args=["matlab-engine", expr],
                returncode=1,
                stdout="",
                stderr=_MATLAB_ENGINE_ERROR or "MATLAB engine initialization failed",
            )
        wrapped = (
            "try; "
            + expr
            + "; "
            + "catch ME; fprintf(2, '%s\\n', getReport(ME, 'extended')); rethrow(ME); end"
        )
        try:
            out = eng.evalc(wrapped, nargout=1)
            return subprocess.CompletedProcess(
                args=["matlab-engine", expr],
                returncode=0,
                stdout=out,
                stderr="",
            )
        except Exception as exc:
            return subprocess.CompletedProcess(
                args=["matlab-engine", expr],
                returncode=1,
                stdout="",
                stderr=str(exc),
            )
    else:
        cmd = ["octave", "--no-gui", "--quiet", "--eval",
               f"addpath('{MATLAB_DIR}'); {expr}"]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def matlab_data_lines(stdout: str) -> list[str]:
    """Return non-warning stdout lines from Octave/MATLAB test snippets."""
    out = []
    in_warning_block = False
    for raw in stdout.splitlines():
        line = raw.replace("\b", "").strip()
        if not line:
            continue
        if line.startswith("[Warning:") or line.startswith("Warning:"):
            in_warning_block = line.startswith("[Warning:") and not line.endswith("]")
            continue
        if in_warning_block:
            if line.endswith("]"):
                in_warning_block = False
            continue
        if line.startswith("> In ") or line.startswith("In "):
            continue
        if line == "]":
            continue
        out.append(line)
    return out


GENOTYPE_CHR1_CALLS = np.array(
    [
        [2, -1, 1, 0],
        [0, 1, 2, -1],
        [1, 1, 0, 2],
        [2, 0, -1, 1],
        [-1, 2, 0, 1],
    ],
    dtype=np.int8,
)

GENOTYPE_SAMPLES = [
    ("FAM1", "IND1", 0, 0, 1, -9),
    ("FAM1", "IND2", 0, 0, 2, -9),
    ("FAM2", "IND3", 0, 0, 1, -9),
    ("FAM2", "IND4", 0, 0, 2, -9),
]

GENOTYPE_CHR1_BIM = [
    ("1", "rs1001", 0, 100, "A", "G"),
    ("1", "rs1002", 0, 200, "C", "T"),
    ("1", "rs1003", 0, 300, "A", "C"),
    ("1", "rs1004", 0, 400, "G", "A"),
    ("1", "rs1005", 0, 500, "T", "C"),
]

GENOTYPE_CHR2_BIM = [
    ("2", "rs2001", 0, 100, "A", "G"),
    ("2", "rs2002", 0, 200, "C", "T"),
]


def write_plink_bed(path: Path, *, num_snp: int, num_sample: int) -> None:
    bytes_per_snp = (num_sample + 3) // 4
    path.write_bytes(b"\x6c\x1b\x01" + b"\x00" * (num_snp * bytes_per_snp))


def write_plink_bed_calls(path: Path, calls) -> None:
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


def write_plink_bim(path: Path, rows) -> None:
    path.write_text("".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n" for r in rows))


def write_plink_fam(path: Path, rows) -> None:
    path.write_text("".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n" for r in rows))


def copy_genotype_shard_files(dst: Path, label: str, *, ploidy: bool = True) -> None:
    for suffix in (".bim", ".fam", ".bed"):
        (dst / f"{label}{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/{label}{suffix}").read_bytes()
        )
    src_ploidy = FIXTURES_DIR / f"genotype/sharded/{label}.ploidy"
    if ploidy and src_ploidy.exists():
        (dst / f"{label}.ploidy").write_bytes(src_ploidy.read_bytes())


def copy_sharded_genotype(dst: Path) -> None:
    for label in ("1", "X"):
        copy_genotype_shard_files(dst, label)
