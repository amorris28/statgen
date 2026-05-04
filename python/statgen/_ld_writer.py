import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from scipy import sparse

from ._ld_npz import validate_npz_distribution
from ._ld_schema import (
    MANIFEST_SCHEMA,
    NPZ_FORMAT,
    PY_RUNTIME_FORMAT,
    SHARD_SCHEMA,
    md5_file,
    validate_chr_sex,
    validate_positive_int,
    validate_shard_metadata,
)


LD_BUILD_METADATA_DEFAULTS = {
    "build_tool": "statgen",
    "build_command": "synthetic/internal LD artifact writer",
    "plink_version": None,
    "ld_window_kb": None,
    "ld_r2_threshold": None,
    "num_sample": None,
}


def write_ld_npz_distribution(root, shard_specs, *, validate=True) -> dict:
    """Write a Python CSC32 LD distribution from already-aligned shard specs."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    shard_specs = list(shard_specs)
    if not shard_specs:
        raise ValueError("LD distribution must contain at least one shard")

    records = []
    seen_keys = set()
    seen_files = set()
    for spec in shard_specs:
        normalized = _normalize_ld_shard_write_spec(spec)
        key = (normalized["chr_label"], normalized["sex"])
        if key in seen_keys:
            raise ValueError(f"duplicate LD shard for chr {key[0]} sex {key[1]!r}")
        seen_keys.add(key)

        file_name = normalized.get("file") or _default_ld_shard_filename(
            normalized["chr_label"],
            normalized["sex"],
            ".npz",
        )
        file_rel = Path(file_name)
        path = root / file_rel
        if file_rel.is_absolute() or ".." in file_rel.parts:
            raise ValueError("LD shard output file must be a relative path inside the distribution")
        file_key = file_rel.as_posix()
        if file_key in seen_files:
            raise ValueError(f"duplicate LD shard output file: {file_key}")
        seen_files.add(file_key)

        shard_payload = {k: v for k, v in normalized.items() if k != "file"}
        meta = _write_ld_npz_shard(path, **shard_payload)
        records.append(
            {
                "chr": meta["chr"],
                "sex": meta["sex"],
                "file": file_name,
                "file_md5": md5_file(path),
                "num_snp": meta["num_snp"],
                "nnz": meta["nnz"],
                "reference_checksum": meta["reference_checksum"],
            }
        )

    manifest = {
        "object_type": "ld_panel_manifest",
        "schema_version": MANIFEST_SCHEMA,
        "runtime_format": PY_RUNTIME_FORMAT,
        "shards": records,
    }
    _write_json(root / "ld_manifest.json", manifest)
    if validate:
        validate_npz_distribution(root, check_payload_structure=True)
    return manifest


def _normalize_ld_shard_write_spec(spec: dict) -> dict:
    if not isinstance(spec, dict):
        raise ValueError("LD shard write spec must be a dict")
    allowed = {
        "reference_shard",
        "chr",
        "sex",
        "num_snp",
        "reference_checksum",
        "a1freq",
        "ld_pairs",
        "build_metadata",
        "extra_metadata",
        "file",
    }
    unknown = sorted(set(spec).difference(allowed))
    if unknown:
        raise ValueError(f"LD shard write spec contains unknown fields: {', '.join(unknown)}")

    ref = spec.get("reference_shard")
    chr_label = str(spec.get("chr", getattr(ref, "label", "")))
    if not chr_label:
        raise ValueError("LD shard write spec must provide chr or reference_shard")
    sex = spec.get("sex")
    if sex is not None:
        sex = str(sex)
    validate_chr_sex(chr_label, sex, "LD shard write spec")

    num_snp = spec.get("num_snp", getattr(ref, "num_snp", None))
    if num_snp is None:
        raise ValueError("LD shard write spec must provide num_snp or reference_shard")
    validate_positive_int(num_snp, "LD shard write spec num_snp")
    num_snp = int(num_snp)

    reference_checksum = spec.get("reference_checksum", getattr(ref, "checksum", None))
    if not isinstance(reference_checksum, str) or reference_checksum == "":
        raise ValueError("LD shard write spec reference_checksum must be a non-empty string")

    return {
        "chr_label": chr_label,
        "sex": sex,
        "num_snp": num_snp,
        "reference_checksum": reference_checksum,
        "a1freq": _coerce_a1freq(spec.get("a1freq"), num_snp),
        "ld_pairs": spec.get("ld_pairs"),
        "build_metadata": dict(spec.get("build_metadata", {})),
        "extra_metadata": dict(spec.get("extra_metadata", {})),
        "file": spec.get("file"),
    }


def _write_ld_npz_shard(
    path,
    *,
    chr_label,
    sex,
    num_snp,
    reference_checksum,
    a1freq,
    ld_pairs,
    build_metadata=None,
    extra_metadata=None,
) -> dict:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    mat = _ld_pairs_to_csc32(ld_pairs, num_snp)
    nnz = int(mat.nnz)
    if nnz >= 2**31:
        raise ValueError(f"{path}: CSC32 nnz must be < 2^31")

    metadata = _ld_npz_metadata(
        chr_label=chr_label,
        sex=sex,
        num_snp=num_snp,
        nnz=nnz,
        reference_checksum=reference_checksum,
        build_metadata=build_metadata,
        extra_metadata=extra_metadata,
    )
    metadata_bytes = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode("utf-8")
    _write_deterministic_npz(
        path,
        {
            "data": mat.data.astype(np.float32, copy=False),
            "indices": mat.indices.astype(np.int32, copy=False),
            "indptr": mat.indptr.astype(np.int32, copy=False),
            "shape": np.array(mat.shape, dtype=np.int64),
            "a1freq": a1freq.astype(np.float32, copy=False),
            "metadata": np.frombuffer(metadata_bytes, dtype=np.uint8),
        },
    )
    return metadata


def _coerce_a1freq(values, num_snp: int) -> np.ndarray:
    if values is None:
        raise ValueError("LD shard write spec must provide a1freq")
    if isinstance(values, pd.DataFrame):
        if "a1freq" not in values.columns:
            raise ValueError("frequency table must contain a1freq")
        values = values["a1freq"]
    elif isinstance(values, dict):
        if "a1freq" not in values:
            raise ValueError("frequency table must contain a1freq")
        values = values["a1freq"]

    arr = np.asarray(values, dtype=np.float32).reshape(-1)
    if arr.size != num_snp:
        raise ValueError(f"a1freq length mismatch: expected {num_snp}, got {arr.size}")
    if not np.all(np.isfinite(arr)):
        raise ValueError("a1freq must contain finite values")
    if np.any((arr < 0.0) | (arr > 1.0)):
        raise ValueError("a1freq values must be between 0 and 1")
    return arr


def _ld_pairs_to_csc32(ld_pairs, num_snp: int) -> sparse.csc_matrix:
    idx1, idx2, r = _coerce_ld_pair_columns(ld_pairs)
    if idx1.size != idx2.size or idx1.size != r.size:
        raise ValueError("LD pair columns must have matching lengths")
    if idx1.size:
        if idx1.min() < 0 or idx2.min() < 0 or idx1.max() >= num_snp or idx2.max() >= num_snp:
            raise ValueError("LD pair indices out of bounds")
        if np.any(idx1 == idx2):
            raise ValueError("LD pair table must not contain diagonal rows")
        if not np.all(np.isfinite(r)):
            raise ValueError("LD r values must be finite")
        if np.any(np.abs(r) > 1.0):
            raise ValueError("LD r values must be between -1 and 1")

        lo = np.minimum(idx1, idx2)
        hi = np.maximum(idx1, idx2)
        pair_keys = lo.astype(np.int64) * np.int64(num_snp) + hi.astype(np.int64)
        if np.unique(pair_keys).size != pair_keys.size:
            raise ValueError("LD pair table contains duplicate unordered pairs")

    rows = np.concatenate([idx1, idx2, np.arange(num_snp, dtype=np.int32)])
    cols = np.concatenate([idx2, idx1, np.arange(num_snp, dtype=np.int32)])
    data = np.concatenate([r, r, np.ones(num_snp, dtype=np.float32)])
    mat = sparse.coo_matrix((data, (rows, cols)), shape=(num_snp, num_snp), dtype=np.float32)
    mat = mat.tocsc()
    mat.sort_indices()
    if mat.indices.dtype != np.int32 or mat.indptr.dtype != np.int32:
        mat.indices = mat.indices.astype(np.int32, copy=False)
        mat.indptr = mat.indptr.astype(np.int32, copy=False)
    return mat


def _coerce_ld_pair_columns(ld_pairs) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if ld_pairs is None:
        return (
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
            np.array([], dtype=np.float32),
        )
    if isinstance(ld_pairs, pd.DataFrame):
        get = ld_pairs.__getitem__
        columns = set(ld_pairs.columns)
    elif isinstance(ld_pairs, dict):
        get = ld_pairs.__getitem__
        columns = set(ld_pairs)
    else:
        arr = np.asarray(ld_pairs)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError("LD pairs must be a table/dict or an n x 3 array")
        return (
            _coerce_ld_index_array(arr[:, 0], "idx1"),
            _coerce_ld_index_array(arr[:, 1], "idx2"),
            np.asarray(arr[:, 2], dtype=np.float32),
        )

    if {"idx1", "idx2", "r"}.issubset(columns):
        return (
            _coerce_ld_index_array(get("idx1"), "idx1"),
            _coerce_ld_index_array(get("idx2"), "idx2"),
            np.asarray(get("r"), dtype=np.float32),
        )
    raise ValueError("LD pairs must be an n x 3 array or a table/dict with idx1, idx2, r columns")


def _coerce_ld_index_array(values, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(arr)) or not np.all(np.floor(arr) == arr):
        raise ValueError(f"LD pair column {name} must contain integer indices")
    if np.any((arr < np.iinfo(np.int32).min) | (arr > np.iinfo(np.int32).max)):
        raise ValueError(f"LD pair column {name} exceeds int32 range")
    return arr.astype(np.int32)


def _ld_npz_metadata(
    *,
    chr_label,
    sex,
    num_snp,
    nnz,
    reference_checksum,
    build_metadata=None,
    extra_metadata=None,
) -> dict:
    meta = {
        "object_type": "ld_shard",
        "schema_version": SHARD_SCHEMA,
        "format": NPZ_FORMAT,
        "chr": chr_label,
        "sex": sex,
        "num_snp": int(num_snp),
        "nnz": int(nnz),
        "sparse_layout": "csc",
        "index_base": 0,
        "matrix": "symmetric",
        "diagonal": "explicit_unit",
        "value": "r",
        "reference_checksum": reference_checksum,
    }
    build = dict(LD_BUILD_METADATA_DEFAULTS)
    if build_metadata:
        build.update(build_metadata)
    meta.update(build)
    if extra_metadata:
        meta.update(extra_metadata)
    validate_shard_metadata(meta, Path("<generated>"), expected_format=NPZ_FORMAT)
    return meta


def _default_ld_shard_filename(chr_label: str, sex, suffix: str) -> str:
    if chr_label == "X":
        return f"ld_chrX_{sex}{suffix}"
    return f"ld_chr{chr_label}{suffix}"


def _write_deterministic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, array in arrays.items():
            buf = io.BytesIO()
            np.save(buf, array, allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, buf.getvalue())


def _write_json(path: Path, data: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
