import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse


_MANIFEST_SCHEMA = "1.0"
_SHARD_SCHEMA = "1.0"
_PY_RUNTIME_FORMAT = "python_npz_csc32"
_MAT_RUNTIME_FORMAT = "matlab_mat_sparse_double"
_NPZ_FORMAT = "statgen_ld_npz_csc32"
_VALID_CHRX_SEX = {"female", "male", "combined"}
_NPZ_REQUIRED = {"data", "indices", "indptr", "shape", "a1freq", "metadata"}


class LDShard:
    def __init__(self, chr_label, sex, num_snp, ld_r, a1freq, reference_checksum):
        self._chr = str(chr_label)
        self._sex = sex
        self._num_snp = int(num_snp)
        self._ld_r = ld_r.tocsc(copy=False)
        self._a1freq = np.asarray(a1freq, dtype=np.float32).reshape(-1)
        self._reference_checksum = str(reference_checksum)

    @property
    def chr(self) -> str:
        return self._chr

    @property
    def label(self) -> str:
        return self._chr

    @property
    def sex(self):
        return self._sex

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def ld_r(self):
        return self._ld_r

    @property
    def a1freq(self) -> np.ndarray:
        return self._a1freq

    @property
    def reference_checksum(self) -> str:
        return self._reference_checksum

    @property
    def checksum(self) -> str:
        return self._reference_checksum


class LDPanel:
    def __init__(self, shard_groups, default_chrX_sex="female"):
        self._shard_groups = [list(g) for g in shard_groups]
        self._default_chrX_sex = _validate_chrx_sex(default_chrX_sex, "default_chrX_sex")
        self._validate_default_chrX_sex()

    @property
    def shard_groups(self) -> list[list[LDShard]]:
        return [list(g) for g in self._shard_groups]

    @property
    def shards(self) -> list[LDShard]:
        out = []
        for group in self._shard_groups:
            if not group:
                continue
            if group[0].chr == "X":
                by_sex = {s.sex: s for s in group}
                out.append(by_sex[self._default_chrX_sex])
            else:
                out.append(group[0])
        return out

    @property
    def default_chrX_sex(self) -> str:
        return self._default_chrX_sex

    def _validate_default_chrX_sex(self) -> None:
        for group in self._shard_groups:
            if group and group[0].chr == "X":
                present = {s.sex for s in group}
                if self._default_chrX_sex not in present:
                    raise ValueError(
                        "default_chrX_sex must name a loaded chrX LD shard; "
                        f"got {self._default_chrX_sex!r}, present {sorted(present)!r}"
                    )


def load_ld(path, reference, default_chrX_sex=None) -> LDPanel:
    default_chrX_sex = _validate_chrx_sex(
        "female" if default_chrX_sex is None else default_chrX_sex,
        "default_chrX_sex",
    )
    path = Path(path)
    ref_shards = list(reference.shards)
    if not ref_shards:
        raise ValueError("load_ld: reference must contain at least one shard")

    if path.is_dir():
        manifest = _read_manifest(path / "ld_manifest.json", expected_runtime=_PY_RUNTIME_FORMAT)
        groups = _load_panel_root(path, manifest, ref_shards)
    else:
        groups = _load_single_shard(path, ref_shards)

    return LDPanel(groups, default_chrX_sex=default_chrX_sex)


def validate_ld_distribution(path, check_payload_structure=False) -> dict:
    path = Path(path)

    if path.is_dir():
        manifest = _read_manifest(path / "ld_manifest.json", expected_runtime=_PY_RUNTIME_FORMAT)
        runtime = manifest["runtime_format"]
        for entry in manifest["shards"]:
            shard_path = path / entry["file"]
            _require_file(shard_path)
            file_md5 = _md5_file(shard_path)
            if file_md5 != entry.get("file_md5"):
                raise ValueError(f"{shard_path}: file_md5 does not match manifest")
            shard, meta = _read_shard_by_runtime(
                shard_path,
                runtime,
                check_payload_structure=check_payload_structure,
            )
            _validate_manifest_entry_agreement(entry, meta, shard_path)
    else:
        runtime = _runtime_from_suffix(path)
        shard, meta = _read_shard_by_runtime(
            path,
            runtime,
            check_payload_structure=check_payload_structure,
        )

    return {"ok": True}


def _load_panel_root(root: Path, manifest: dict, ref_shards: list) -> list[list[LDShard]]:
    entries = list(manifest["shards"])
    groups = []
    for ref_shard in ref_shards:
        if ref_shard.label == "X":
            selected = [e for e in entries if e.get("chr") == "X"]
            if not selected:
                raise FileNotFoundError("LD manifest has no chrX shard for reference shard X")
        else:
            selected = [
                e
                for e in entries
                if e.get("chr") == ref_shard.label and e.get("sex") is None
            ]
            if len(selected) != 1:
                raise FileNotFoundError(
                    f"LD manifest must contain exactly one autosomal shard for chr {ref_shard.label}"
                )

        group = []
        seen_sex = set()
        for entry in selected:
            shard_path = root / entry["file"]
            shard, meta = _load_npz_shard(shard_path)
            _validate_manifest_entry_agreement(entry, meta, shard_path)
            _validate_reference_compatibility(shard, ref_shard, shard_path)
            key = shard.sex
            if key in seen_sex:
                raise ValueError(f"LD manifest has duplicate shard for chr {shard.chr} sex {key!r}")
            seen_sex.add(key)
            group.append(shard)
        groups.append(group)
    return groups


def _load_single_shard(path: Path, ref_shards: list) -> list[list[LDShard]]:
    shard, _meta = _load_npz_shard(path)
    matching = [s for s in ref_shards if s.label == shard.chr]
    if len(ref_shards) != 1 or len(matching) != 1:
        raise ValueError(
            "single-shard LD loads require a single-shard reference with the same chromosome"
        )
    _validate_reference_compatibility(shard, matching[0], path)
    return [[shard]]


def _read_manifest(path: Path, expected_runtime: str | None) -> dict:
    _require_file(path)
    with open(path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    if manifest.get("object_type") != "ld_panel_manifest":
        raise ValueError(f"{path}: object_type must be 'ld_panel_manifest'")
    if manifest.get("schema_version") != _MANIFEST_SCHEMA:
        raise ValueError(f"{path}: unsupported LD manifest schema_version {manifest.get('schema_version')!r}")
    runtime = manifest.get("runtime_format")
    if runtime not in {_PY_RUNTIME_FORMAT, _MAT_RUNTIME_FORMAT}:
        raise ValueError(f"{path}: unsupported LD runtime_format {runtime!r}")
    if expected_runtime is not None and runtime != expected_runtime:
        raise ValueError(f"{path}: expected runtime_format {expected_runtime!r}, got {runtime!r}")
    shards = manifest.get("shards")
    if not isinstance(shards, list) or not shards:
        raise ValueError(f"{path}: shards must be a non-empty list")

    seen = set()
    for i, entry in enumerate(shards):
        _validate_manifest_entry(entry, f"{path}:shards[{i}]")
        key = (entry["chr"], entry["sex"])
        if key in seen:
            raise ValueError(f"{path}: duplicate manifest entry for chr {key[0]} sex {key[1]!r}")
        seen.add(key)
    return manifest


def _validate_manifest_entry(entry: dict, where: str) -> None:
    required = {"chr", "sex", "file", "file_md5", "num_snp", "nnz", "reference_checksum"}
    missing = sorted(required.difference(entry))
    if missing:
        raise ValueError(f"{where}: missing required fields: {', '.join(missing)}")
    _validate_chr_sex(entry["chr"], entry["sex"], where)
    if not isinstance(entry["file"], str) or entry["file"] == "":
        raise ValueError(f"{where}: file must be a non-empty relative path")
    if Path(entry["file"]).is_absolute() or ".." in Path(entry["file"]).parts:
        raise ValueError(f"{where}: file must be a relative path inside the LD distribution")
    _validate_positive_int(entry["num_snp"], f"{where}: num_snp")
    _validate_nonnegative_int(entry["nnz"], f"{where}: nnz")
    if not isinstance(entry["reference_checksum"], str) or entry["reference_checksum"] == "":
        raise ValueError(f"{where}: reference_checksum must be a non-empty string")
    if not isinstance(entry["file_md5"], str) or len(entry["file_md5"]) != 32:
        raise ValueError(f"{where}: file_md5 must be a lowercase MD5 hex string")


def _load_npz_shard(path: Path) -> tuple[LDShard, dict]:
    shard, meta, _warning = _read_npz_shard(path, check_payload_structure=False)
    return shard, meta


def _read_shard_by_runtime(path: Path, runtime: str, check_payload_structure: bool):
    if runtime == _PY_RUNTIME_FORMAT:
        shard, meta, _unused = _read_npz_shard(path, check_payload_structure)
        return shard, meta
    raise ValueError(f"{path}: unsupported LD runtime {runtime!r}")


def _read_npz_shard(path: Path, check_payload_structure: bool):
    _require_file(path)
    with np.load(path, allow_pickle=False) as data:
        names = set(data.files)
        missing = sorted(_NPZ_REQUIRED.difference(names))
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
        _validate_shard_metadata(meta, path, expected_format=_NPZ_FORMAT)

        num_snp = int(meta["num_snp"])
        nnz = int(meta["nnz"])
        _validate_npz_dimensions(path, arr_data, indices, indptr, shape, a1freq, num_snp, nnz)

        if check_payload_structure:
            _validate_csc_payload(path, arr_data, indices, indptr, num_snp, nnz)

        ld_r = sparse.csc_matrix(
            (arr_data.copy(), indices.copy(), indptr.copy()),
            shape=(num_snp, num_snp),
        )
        if check_payload_structure:
            _validate_sparse_matrix_values(path, ld_r)

        shard = LDShard(
            meta["chr"],
            meta["sex"],
            num_snp,
            ld_r,
            a1freq.copy(),
            meta["reference_checksum"],
        )
        return shard, meta, None


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


def _validate_shard_metadata(meta: dict, path: Path, expected_format: str) -> None:
    required = {
        "object_type",
        "schema_version",
        "format",
        "chr",
        "sex",
        "num_snp",
        "nnz",
        "matrix",
        "diagonal",
        "value",
        "reference_checksum",
    }
    missing = sorted(required.difference(meta))
    if missing:
        raise ValueError(f"{path}: metadata missing required fields: {', '.join(missing)}")
    if meta["object_type"] != "ld_shard":
        raise ValueError(f"{path}: metadata object_type must be 'ld_shard'")
    if meta["schema_version"] != _SHARD_SCHEMA:
        raise ValueError(f"{path}: unsupported metadata schema_version {meta['schema_version']!r}")
    if meta["format"] != expected_format:
        raise ValueError(f"{path}: metadata format must be {expected_format!r}")
    _validate_chr_sex(meta["chr"], meta["sex"], f"{path}: metadata")
    _validate_positive_int(meta["num_snp"], f"{path}: metadata num_snp")
    _validate_nonnegative_int(meta["nnz"], f"{path}: metadata nnz")
    if int(meta["nnz"]) >= 2**31:
        raise ValueError(f"{path}: CSC32 metadata nnz must be < 2^31")
    if meta["matrix"] != "symmetric":
        raise ValueError(f"{path}: metadata matrix must be 'symmetric'")
    if meta["diagonal"] != "explicit_unit":
        raise ValueError(f"{path}: metadata diagonal must be 'explicit_unit'")
    if meta["value"] != "r":
        raise ValueError(f"{path}: metadata value must be 'r'")
    if not isinstance(meta["reference_checksum"], str) or meta["reference_checksum"] == "":
        raise ValueError(f"{path}: metadata reference_checksum must be a non-empty string")
    if expected_format == _NPZ_FORMAT:
        if meta.get("sparse_layout") != "csc":
            raise ValueError(f"{path}: metadata sparse_layout must be 'csc'")
        if meta.get("index_base") != 0:
            raise ValueError(f"{path}: metadata index_base must be 0")


def _validate_chr_sex(chr_label, sex, where: str) -> None:
    if chr_label not in [str(i) for i in range(1, 23)] + ["X"]:
        raise ValueError(f"{where}: chr must be one of 1-22 or X")
    if chr_label == "X":
        _validate_chrx_sex(sex, f"{where}: sex")
    elif sex is not None:
        raise ValueError(f"{where}: autosomal LD shards must have sex null")


def _validate_chrx_sex(sex, where: str) -> str:
    if sex not in _VALID_CHRX_SEX:
        raise ValueError(f"{where}: chrX sex must be one of female, male, combined")
    return sex


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


def _validate_manifest_entry_agreement(entry: dict, meta: dict, path: Path) -> None:
    for key in ("chr", "sex", "num_snp", "nnz", "reference_checksum"):
        if entry.get(key) != meta.get(key):
            raise ValueError(f"{path}: manifest/per-file metadata mismatch for {key}")
    suffix_runtime = _runtime_from_suffix(path)
    if suffix_runtime == _PY_RUNTIME_FORMAT and meta.get("format") != _NPZ_FORMAT:
        raise ValueError(f"{path}: manifest file extension and metadata format disagree")


def _validate_reference_compatibility(shard: LDShard, ref_shard, path: Path) -> None:
    if shard.chr != ref_shard.label:
        raise ValueError(f"{path}: LD chr {shard.chr!r} does not match reference shard {ref_shard.label!r}")
    if shard.num_snp != ref_shard.num_snp:
        raise ValueError(f"{path}: LD num_snp does not match reference shard")
    if shard.reference_checksum != ref_shard.checksum:
        raise ValueError(f"{path}: LD reference_checksum does not match reference shard")


def _runtime_from_suffix(path: Path) -> str:
    if path.suffix == ".npz":
        return _PY_RUNTIME_FORMAT
    raise ValueError(f"{path}: unsupported LD shard extension")


def _require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"LD file not found: {path}")


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _validate_positive_int(value, where: str) -> None:
    if not isinstance(value, (int, np.integer)) or int(value) <= 0:
        raise ValueError(f"{where} must be a positive integer")


def _validate_nonnegative_int(value, where: str) -> None:
    if not isinstance(value, (int, np.integer)) or int(value) < 0:
        raise ValueError(f"{where} must be a non-negative integer")
