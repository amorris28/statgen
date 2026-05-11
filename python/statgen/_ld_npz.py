import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ._ld_schema import (
    NPZ_FORMAT,
    NPZ_REQUIRED,
    PY_RUNTIME_FORMAT,
    md5_file,
    read_manifest,
    require_file,
    runtime_from_suffix,
    validate_manifest_entry_agreement,
    validate_shard_metadata,
)


def validate_npz_distribution(path, check_payload_structure=False) -> dict:
    path = Path(path)
    if not path.is_dir():
        raise ValueError("validate_npz_distribution: path must identify a panel root directory")

    manifest = read_manifest(path / "ld_manifest.json", expected_runtime=PY_RUNTIME_FORMAT)
    reference_cache_path = path / manifest["reference_cache"]
    require_file(reference_cache_path)
    if md5_file(reference_cache_path) != manifest["reference_cache_md5"]:
        raise ValueError(f"{reference_cache_path}: reference_cache_md5 does not match manifest")
    from .reference import load_reference_cache

    reference = load_reference_cache(reference_cache_path)
    ref_by_label = {s.label: s for s in reference.shards}
    seen_reference_bim = set()
    for entry in manifest["shards"]:
        shard_path = path / entry["file"]
        require_file(shard_path)
        file_md5 = md5_file(shard_path)
        if file_md5 != entry.get("file_md5"):
            raise ValueError(f"{shard_path}: file_md5 does not match manifest")
        _ld_r, _a1freq, meta = read_npz_shard(
            shard_path,
            check_payload_structure=check_payload_structure,
        )
        validate_manifest_entry_agreement(entry, meta, shard_path)
        ref_shard = ref_by_label.get(entry["chr"])
        if ref_shard is None:
            raise ValueError(f"{reference_cache_path}: reference cache missing shard {entry['chr']}")
        _validate_reference_cache_entry(ref_shard, entry, reference_cache_path)
        reference_key = (
            entry["reference_bim"],
            entry["chr"],
            int(entry["num_snp"]),
            entry["reference_checksum"],
        )
        if reference_key not in seen_reference_bim:
            _validate_bundled_reference(root=path, entry=entry)
            seen_reference_bim.add(reference_key)

    return {"ok": True}


def validate_npz_shard_file(path, check_payload_structure=False) -> dict:
    path = Path(path)
    runtime = runtime_from_suffix(path)
    if runtime != PY_RUNTIME_FORMAT:
        raise ValueError(f"{path}: unsupported LD runtime {runtime!r}")
    read_npz_shard(path, check_payload_structure=check_payload_structure)
    return {"ok": True}


def _validate_bundled_reference(*, root: Path, entry: dict) -> None:
    from ._ld_reference import validate_bundled_reference_bim

    validate_bundled_reference_bim(root / entry["reference_bim"], entry, target="manifest")


def _validate_reference_cache_entry(ref_shard, entry: dict, path: Path) -> None:
    if ref_shard.num_snp != int(entry["num_snp"]):
        raise ValueError(f"{path}: reference cache num_snp does not match manifest")
    if ref_shard.checksum != entry["reference_checksum"]:
        raise ValueError(f"{path}: reference cache reference_checksum does not match manifest")


def read_npz_shard(path: Path, check_payload_structure: bool):
    require_file(path)
    with np.load(path, allow_pickle=False) as data:
        names = set(data.files)
        missing = sorted(NPZ_REQUIRED.difference(names))
        if missing:
            raise ValueError(f"{path}: missing required arrays: {', '.join(missing)}")

        arr_data = data["data"]
        indices = data["indices"]
        indptr = data["indptr"]
        shape = data["shape"]
        a1freq = data["a1freq"]
        metadata = data["metadata"]

        _validate_npz_array_dtypes(path, arr_data, indices, indptr, shape, a1freq, metadata)
        meta = _decode_npz_metadata(metadata, path)
        validate_shard_metadata(meta, path, expected_format=NPZ_FORMAT)

        num_snp = int(meta["num_snp"])
        nnz = int(meta["nnz"])
        _validate_npz_dimensions(path, arr_data, indices, indptr, shape, a1freq, num_snp, nnz)

        if check_payload_structure:
            _validate_csc_payload(path, arr_data, indices, indptr, num_snp, nnz)

        ld_r = sparse.csc_matrix(
            (arr_data, indices, indptr),
            shape=(num_snp, num_snp),
        )
        if check_payload_structure:
            _validate_sparse_matrix_values(path, ld_r)

        return ld_r, a1freq, meta


def _validate_npz_array_dtypes(path, data, indices, indptr, shape, a1freq, metadata) -> None:
    expected = [
        ("data", data.dtype, np.dtype("float32")),
        ("indices", indices.dtype, np.dtype("int32")),
        ("indptr", indptr.dtype, np.dtype("int32")),
        ("shape", shape.dtype, np.dtype("int64")),
        ("a1freq", a1freq.dtype, np.dtype("float32")),
        ("metadata", metadata.dtype, np.dtype("uint8")),
    ]
    for name, got, want in expected:
        if got != want:
            raise ValueError(f"{path}: {name} must have dtype {want}, got {got}")
    for name, arr in [
        ("data", data),
        ("indices", indices),
        ("indptr", indptr),
        ("shape", shape),
        ("a1freq", a1freq),
        ("metadata", metadata),
    ]:
        if arr.dtype.hasobject:
            raise ValueError(f"{path}: {name} must not be an object array")


def _decode_npz_metadata(metadata: np.ndarray, path: Path) -> dict:
    if metadata.ndim != 1:
        raise ValueError(f"{path}: metadata must be a one-dimensional uint8 array")
    try:
        out = json.loads(bytes(metadata).decode("utf-8"))
    except Exception as exc:
        raise ValueError(f"{path}: metadata must be UTF-8 JSON bytes") from exc
    if not isinstance(out, dict):
        raise ValueError(f"{path}: metadata JSON must be an object")
    return out


def _validate_npz_dimensions(path, data, indices, indptr, shape, a1freq, num_snp, nnz) -> None:
    if shape.ndim != 1 or shape.size != 2:
        raise ValueError(f"{path}: shape must be length 2")
    if tuple(int(x) for x in shape) != (num_snp, num_snp):
        raise ValueError(f"{path}: shape must equal [num_snp, num_snp]")
    if indptr.ndim != 1 or indptr.size != num_snp + 1:
        raise ValueError(f"{path}: indptr length must be num_snp + 1")
    if data.ndim != 1 or indices.ndim != 1:
        raise ValueError(f"{path}: data and indices must be one-dimensional")
    if data.size != indices.size or data.size != nnz:
        raise ValueError(f"{path}: data, indices, and metadata nnz length mismatch")
    if int(indptr[-1]) != nnz:
        raise ValueError(f"{path}: indptr[-1] must equal metadata nnz")
    if a1freq.ndim != 1 or a1freq.size != num_snp:
        raise ValueError(f"{path}: a1freq length must equal num_snp")


def _validate_csc_payload(path, data, indices, indptr, num_snp, nnz) -> None:
    if indptr[0] != 0:
        raise ValueError(f"{path}: indptr[0] must be 0")
    if np.any(indptr[1:] < indptr[:-1]):
        raise ValueError(f"{path}: indptr must be monotonic nondecreasing")
    if indices.size and (indices.min() < 0 or indices.max() >= num_snp):
        raise ValueError(f"{path}: sparse row indices out of bounds")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"{path}: sparse data must be finite")
    if nnz < num_snp:
        raise ValueError(f"{path}: sparse matrix must include explicit unit diagonal")


def _validate_sparse_matrix_values(path, mat) -> None:
    if mat.shape[0] != mat.shape[1]:
        raise ValueError(f"{path}: ld_r must be square")
    diag = mat.diagonal()
    if diag.size != mat.shape[0] or not np.allclose(diag, 1.0, rtol=0, atol=1e-7):
        raise ValueError(f"{path}: ld_r diagonal must be explicit unit")
    mat_csr = mat.tocsr()
    for i in range(mat.shape[0]):
        start, stop = mat_csr.indptr[i], mat_csr.indptr[i + 1]
        if not np.any((mat_csr.indices[start:stop] == i) & np.isclose(mat_csr.data[start:stop], 1.0)):
            raise ValueError(f"{path}: ld_r diagonal must be explicit unit")
    diff = (mat - mat.T).tocoo()
    if diff.nnz and np.max(np.abs(diff.data)) > 1e-6:
        raise ValueError(f"{path}: ld_r must be symmetric")
