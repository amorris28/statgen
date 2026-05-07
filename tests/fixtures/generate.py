#!/usr/bin/env python3
"""
Generate all committed test fixtures under tests/fixtures/.

Run from anywhere:
    python tests/fixtures/generate.py

Binary fixture arrays use the dtypes documented by their object specs.
Reference checksums are MD5 over "chr:bp:a1:a2\n" lines in shard row order.
"""

import gzip
import hashlib
import io
import json
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import savemat

ROOT = Path(__file__).parent
PYTHON_ROOT = ROOT.parents[1] / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from statgen._ld_schema import md5_file


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bim_checksum(rows: list[tuple]) -> str:
    """MD5 over 'chr:bp:a1:a2\n' lines. rows = [(chr, snp, cm, bp, a1, a2), ...]"""
    text = "".join(f"{r[0]}:{r[3]}:{r[4]}:{r[5]}\n" for r in rows)
    return hashlib.md5(text.encode()).hexdigest()


def write_bim(path: Path, rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def write_bed_annotation(path: Path, intervals: list[tuple]) -> None:
    """Write a BED annotation file (chrom, start, end) — not PLINK .bed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")


def write_plink_bed(path: Path, n_snp: int, n_sample: int) -> None:
    """Write a minimal PLINK .bed file (magic + SNP-major + all homref genotypes)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    bytes_per_snp = (n_sample + 3) // 4
    with open(path, "wb") as f:
        f.write(b"\x6c\x1b\x01")           # magic + SNP-major mode
        f.write(b"\x00" * (n_snp * bytes_per_snp))


def write_ploidy(path: Path, rows: list[tuple]) -> None:
    """Write a .ploidy sidecar: tab-delimited, no header, (male_ploidy, female_ploidy) per row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for male, female in rows:
            f.write(f"{male}\t{female}\n")


def write_fam(path: Path, samples: list[tuple]) -> None:
    """Write a PLINK .fam file. samples = [(fid, iid, pid, mid, sex, pheno), ...]"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for s in samples:
            f.write("\t".join(str(x) for x in s) + "\n")


def ld_sparse_matrix(
    num_snp: int,
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    dtype,
) -> sparse.csc_matrix:
    rows = np.array(idx1 + idx2 + list(range(num_snp)), dtype=np.int32)
    cols = np.array(idx2 + idx1 + list(range(num_snp)), dtype=np.int32)
    data = np.array(r + r + [1.0] * num_snp, dtype=dtype)
    return sparse.coo_matrix((data, (rows, cols)), shape=(num_snp, num_snp)).tocsc()


def ld_metadata(
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    nnz: int,
    runtime_format: str,
    reference_bim: str | None = None,
    extra_meta: dict | None = None,
) -> dict:
    if reference_bim is None:
        reference_bim = f"reference_chr{chr_label}.bim"
    meta = {
        "object_type": "ld_shard",
        "schema_version": "1.0",
        "format": runtime_format,
        "chr": chr_label,
        "sex": sex,
        "num_snp": len(bim_rows),
        "nnz": nnz,
        "matrix": "symmetric",
        "diagonal": "explicit_unit",
        "value": "r",
        "num_monomorphic_snps": 0,
        "reference_checksum": bim_checksum(bim_rows),
        "reference_bim": reference_bim,
        "build_tool": "generate.py",
        "build_command": "synthetic fixture",
        "plink_version": None,
        "ld_window_kb": 10000,
        "ld_r2_threshold": 0.05,
        "num_sample": 100,
    }
    if runtime_format == "statgen_ld_npz_csc32":
        meta.update({"sparse_layout": "csc", "index_base": 0})
    if extra_meta:
        meta.update(extra_meta)
    return meta


def write_deterministic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, array in arrays.items():
            buf = io.BytesIO()
            np.save(buf, array, allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, buf.getvalue())


def write_ld_npz_shard(
    path: Path,
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    a1freq: list[float],
    extra_meta: dict | None = None,
) -> dict:
    mat = ld_sparse_matrix(len(bim_rows), idx1, idx2, r, np.float32)
    meta = ld_metadata(
        chr_label,
        sex,
        bim_rows,
        nnz=int(mat.nnz),
        runtime_format="statgen_ld_npz_csc32",
        extra_meta=extra_meta,
    )
    metadata_bytes = json.dumps(meta, sort_keys=True, separators=(",", ":")).encode()
    write_deterministic_npz(
        path,
        {
            "data": mat.data.astype(np.float32, copy=False),
            "indices": mat.indices.astype(np.int32, copy=False),
            "indptr": mat.indptr.astype(np.int32, copy=False),
            "shape": np.array(mat.shape, dtype=np.int64),
            "a1freq": np.array(a1freq, dtype=np.float32),
            "metadata": np.frombuffer(metadata_bytes, dtype=np.uint8),
        },
    )
    return meta


def matlab_metadata(meta: dict) -> dict:
    return {
        key: (np.array([], dtype=np.float64) if value is None else value)
        for key, value in meta.items()
    }


def write_ld_mat_shard(
    path: Path,
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    a1freq: list[float],
    extra_meta: dict | None = None,
) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    mat = ld_sparse_matrix(len(bim_rows), idx1, idx2, r, np.float64)
    meta = ld_metadata(
        chr_label,
        sex,
        bim_rows,
        nnz=int(mat.nnz),
        runtime_format="statgen_ld_mat_sparse_double",
        extra_meta=extra_meta,
    )
    savemat(
        path,
        {
            "ld_r": mat,
            "a1freq": np.array(a1freq, dtype=np.float64).reshape(-1, 1),
            "metadata": matlab_metadata(meta),
        },
        appendmat=False,
        do_compression=False,
        long_field_names=True,
    )
    deterministic_header = (
        "MATLAB 5.0 MAT-file, Platform: statgen, Created by tests/fixtures/generate.py"
    )
    with open(path, "r+b") as f:
        f.write(deterministic_header.encode("ascii")[:116].ljust(116, b" "))
    return meta


def write_ld_manifest(directory: Path, runtime_format: str, shard_records: list[dict]) -> None:
    manifest = {
        "object_type": "ld_panel_manifest",
        "schema_version": "1.0",
        "runtime_format": runtime_format,
        "shards": shard_records,
    }
    with open(directory / "ld_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")


def write_ld_distribution(directory: Path, runtime_format: str, extension: str, shard_writer) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    reference_bims = {
        "1": "reference_chr1.bim",
        "X": "reference_chrX.bim",
    }
    write_bim(directory / reference_bims["1"], CHR1_BIM)
    write_bim(directory / reference_bims["X"], CHRX_BIM)
    shard_specs = [
        ("1", None, f"ld_chr1.{extension}", CHR1_BIM, CHR1_LD, None),
        ("X", "female", f"ld_chrX_female.{extension}", CHRX_BIM, CHRX_LD["female"], None),
        ("X", "male", f"ld_chrX_male.{extension}", CHRX_BIM, CHRX_LD["male"], None),
        (
            "X",
            "combined",
            f"ld_chrX_combined.{extension}",
            CHRX_BIM,
            CHRX_LD["combined"],
            {"chrX_combined_rationale": "synthetic fixture coverage"},
        ),
    ]
    records = []
    for chr_label, sex, file_name, bim_rows, data, extra_meta in shard_specs:
        path = directory / file_name
        idx1 = data.get("idx1", CHRX_LD_IDX1)
        idx2 = data.get("idx2", CHRX_LD_IDX2)
        meta = shard_writer(
            path,
            chr_label,
            sex,
            bim_rows,
            idx1,
            idx2,
            data["r"],
            data["a1freq"],
            extra_meta,
        )
        records.append({
            "chr": chr_label,
            "sex": sex,
            "file": file_name,
            "file_md5": md5_file(path),
            "num_snp": meta["num_snp"],
            "nnz": meta["nnz"],
            "reference_checksum": meta["reference_checksum"],
            "reference_bim": reference_bims[chr_label],
        })
    write_ld_manifest(directory, runtime_format, records)


# ---------------------------------------------------------------------------
# Variant definitions
# ---------------------------------------------------------------------------

# chr1: 5 SNPs — columns: chr, snp, cm, bp, a1, a2
CHR1_BIM = [
    ("1", "rs1001", 0, 100, "A", "G"),
    ("1", "rs1002", 0, 200, "C", "T"),
    ("1", "rs1003", 0, 300, "A", "C"),
    ("1", "rs1004", 0, 400, "G", "A"),
    ("1", "rs1005", 0, 500, "T", "C"),
]

# chrX: 3 SNPs
CHRX_BIM = [
    ("X", "rsX001", 0, 100, "A", "G"),
    ("X", "rsX002", 0, 200, "C", "T"),
    ("X", "rsX003", 0, 300, "A", "C"),
]

ALL_BIM = CHR1_BIM + CHRX_BIM

SAMPLES = [
    ("FAM1", "IND1", 0, 0, 1, -9),   # male
    ("FAM1", "IND2", 0, 0, 2, -9),   # female
    ("FAM2", "IND3", 0, 0, 1, -9),   # male
    ("FAM2", "IND4", 0, 0, 2, -9),   # female
]

# Ploidy: (male_ploidy, female_ploidy) per BIM row.
# chrX is hemizygous for males (1) and diploid for females (2).
CHRX_PLOIDY = [(1, 2)] * len(CHRX_BIM)
ALL_PLOIDY = [(2, 2)] * len(CHR1_BIM) + CHRX_PLOIDY

# ---------------------------------------------------------------------------
# Fixture data
# ---------------------------------------------------------------------------

# LD for chr1 (upper-triangle pairs, idx1 < idx2)
CHR1_LD = dict(
    idx1=[0, 0, 1, 3],
    idx2=[1, 2, 3, 4],
    r=[0.9, 0.3, -0.5, 0.7],
    a1freq=[0.30, 0.40, 0.20, 0.35, 0.15],
)

# LD for chrX — same sparse structure, sex-specific a1freq and r
CHRX_LD_IDX1 = [0, 1]
CHRX_LD_IDX2 = [1, 2]
CHRX_LD = {
    "female":   dict(r=[0.65, -0.35], a1freq=[0.28, 0.32, 0.22]),
    "male":     dict(r=[0.55, -0.45], a1freq=[0.22, 0.28, 0.18]),
    "combined": dict(r=[0.60, -0.40], a1freq=[0.25, 0.30, 0.20]),
}

# anno1.bed — chr1; overlapping [99,150)+[120,200) and adjacent [299,350)+[350,400)
#   SNP bp→0-based: 100→99, 200→199, 300→299, 400→399, 500→499
#   [99,150) ∪ [120,200) → after merge: [99,200)  covers rs1001(99), rs1002(199)
#   [299,350) + [350,400) → adjacent, merged: [299,400)  covers rs1003(299), rs1004(399)
#   [499,500) → boundary: covers rs1005(499) exactly
ANNO1_BED = [
    ("1", 99,  150),   # overlapping pair ↓
    ("1", 120, 200),
    ("1", 299, 350),   # adjacent pair ↓
    ("1", 350, 400),
    ("1", 499, 500),   # boundary: exactly covers rs1005
]

# anno2.bed — chrX; overlapping [195,205)+[198,302)
#   rsX001(99): covered by [99,100)
#   rsX002(199): covered by merged [195,302)
#   rsX003(299): covered by merged [195,302)
ANNO2_BED = [
    ("X",  99, 100),   # covers rsX001 only
    ("X", 195, 205),   # overlapping pair ↓
    ("X", 198, 302),
]

# Sumstats: chr, bp, a1, a2, z, n, p
# Row notes:
#   1:200  z=NA  (missing z)
#   1:300  p=0   (logp → Inf)
#   1:400  p=-0.1 (invalid p → logp NaN)
#   9:999  absent from reference → NaN row in aligned output
#   X:300 (rsX003) absent from file → NaN row in aligned output
SUMSTATS_ROWS = """\
chr\tbp\ta1\ta2\tz\tn\tp
1\t100\tA\tG\t2.5\t1000\t0.012
1\t200\tC\tT\tNA\t1000\t0.5
1\t300\tA\tC\t1.8\t1000\t0
1\t400\tG\tA\t-1.2\t1000\t-0.1
X\t100\tA\tG\t3.0\t500\t0.003
X\t200\tC\tT\t0.5\t500\t0.6
9\t999\tA\tG\t1.0\t1000\t0.3
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # --- reference ---
    write_bim(ROOT / "reference/sharded/1.bim", CHR1_BIM)
    write_bim(ROOT / "reference/sharded/X.bim", CHRX_BIM)
    write_bim(ROOT / "reference/nonsharded/all.bim", ALL_BIM)

    # --- annotations ---
    write_bed_annotation(ROOT / "annotations/anno1.bed", ANNO1_BED)
    write_bed_annotation(ROOT / "annotations/anno2.bed", ANNO2_BED)

    # --- sumstats ---
    sumstats_path = ROOT / "sumstats/traits.tsv.gz"
    sumstats_path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 makes the gzip header deterministic (no embedded timestamp).
    with open(sumstats_path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            gz.write(SUMSTATS_ROWS.encode())

    # --- LD runtime distributions ---
    shutil.rmtree(ROOT / "ld", ignore_errors=True)
    write_ld_distribution(
        ROOT / "ld/python",
        runtime_format="python_npz_csc32",
        extension="npz",
        shard_writer=write_ld_npz_shard,
    )
    write_ld_distribution(
        ROOT / "ld/matlab",
        runtime_format="matlab_mat_sparse_double",
        extension="mat",
        shard_writer=write_ld_mat_shard,
    )

    # --- genotype (PLINK bfile metadata) ---
    for label, bim_rows in [("1", CHR1_BIM), ("X", CHRX_BIM)]:
        base = ROOT / f"genotype/sharded/{label}"
        write_bim(Path(str(base) + ".bim"), bim_rows)
        write_fam(Path(str(base) + ".fam"), SAMPLES)
        write_plink_bed(
            Path(str(base) + ".bed"),
            n_snp=len(bim_rows),
            n_sample=len(SAMPLES),
        )
    write_ploidy(ROOT / "genotype/sharded/X.ploidy", CHRX_PLOIDY)

    # --- genotype non-sharded bfile ---
    ns_base = ROOT / "genotype/nonsharded/all"
    write_bim(Path(str(ns_base) + ".bim"), ALL_BIM)
    write_fam(Path(str(ns_base) + ".fam"), SAMPLES)
    write_plink_bed(Path(str(ns_base) + ".bed"), n_snp=len(ALL_BIM), n_sample=len(SAMPLES))
    write_ploidy(Path(str(ns_base) + ".ploidy"), ALL_PLOIDY)

    print("Fixtures written to", ROOT)
    print(f"  chr1 reference checksum : {bim_checksum(CHR1_BIM)}")
    print(f"  chrX reference checksum : {bim_checksum(CHRX_BIM)}")
    print(f"  all  reference checksum : {bim_checksum(ALL_BIM)}")


if __name__ == "__main__":
    main()
