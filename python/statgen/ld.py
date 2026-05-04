from pathlib import Path
import warnings

import numpy as np
from scipy import sparse

from ._ld_npz import read_npz_shard, validate_npz_distribution
from ._ld_schema import (
    PY_RUNTIME_FORMAT,
    read_manifest,
    require_file,
    validate_chrx_sex,
    validate_manifest_entry_agreement,
)
from ._utils import validate_requested_shards


class LDShard:
    def __init__(
        self,
        chr_label,
        sex,
        num_snp,
        ld_r,
        a1freq,
        reference_checksum,
        num_monomorphic_snps=0,
    ):
        self._chr = str(chr_label)
        self._sex = sex
        self._num_snp = int(num_snp)
        self._ld_r = ld_r.tocsc(copy=False)
        self._ld_r2 = _build_ld_r2(self._ld_r)
        self._a1freq = np.asarray(a1freq, dtype=np.float32).reshape(-1)
        self._reference_checksum = str(reference_checksum)
        self._num_monomorphic_snps = int(num_monomorphic_snps)

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

    @property
    def num_monomorphic_snps(self) -> int:
        return self._num_monomorphic_snps


class LDPanel:
    def __init__(self, shard_groups, default_chrX_sex="female"):
        self._shard_groups = [list(g) for g in shard_groups]
        self._default_chrX_sex = validate_chrx_sex(default_chrX_sex, "default_chrX_sex")
        self._shard_offsets = []
        pos = 0
        for group in self._shard_groups:
            if not group:
                raise ValueError("LDPanel shard groups must be non-empty")
            first = group[0]
            for shard in group:
                if shard.chr != first.chr:
                    raise ValueError("all LD shards in a group must share chr")
                if shard.num_snp != first.num_snp:
                    raise ValueError("all LD shards in a group must share num_snp")
            self._shard_offsets.append(
                {"shard_label": first.label, "start0": pos, "stop0": pos + first.num_snp}
            )
            pos += first.num_snp
        self._num_snp = pos
        self._validate_default_chrX_sex()

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def shard_offsets(self) -> list[dict]:
        return list(self._shard_offsets)

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

    def a1freq(self, chrX_sex=None) -> np.ndarray:
        shards = [shard for shard, _start, _stop in self.iter_shards_with_offsets(chrX_sex)]
        if not shards:
            return np.array([], dtype=np.float32)
        return np.concatenate([s.a1freq for s in shards]).astype(np.float32, copy=False)

    def select_shards(self, shards) -> "LDPanel":
        available = [group[0].label for group in self._shard_groups]
        selected = validate_requested_shards(shards, available, "LDPanel.select_shards")
        by_label = {group[0].label: group for group in self._shard_groups}
        return LDPanel([by_label[label] for label in selected], self._default_chrX_sex)

    def multiply_r2(self, M, chrX_sex=None):
        if sparse.issparse(M):
            return self._multiply_r2_sparse(M, chrX_sex)
        return self._multiply_r2_dense(M, chrX_sex)

    def iter_shards_with_offsets(self, chrX_sex=None):
        for group, offset in zip(self._shard_groups, self._shard_offsets):
            first = group[0]
            if first.chr != "X":
                yield first, int(offset["start0"]), int(offset["stop0"])
                continue

            sex = (
                self._default_chrX_sex
                if chrX_sex is None
                else validate_chrx_sex(chrX_sex, "chrX_sex")
            )
            by_sex = {s.sex: s for s in group}
            try:
                yield by_sex[sex], int(offset["start0"]), int(offset["stop0"])
            except KeyError as exc:
                raise ValueError(
                    "chrX_sex must name a loaded chrX LD shard; "
                    f"got {sex!r}, present {sorted(by_sex)!r}"
                ) from exc

    def _validate_default_chrX_sex(self) -> None:
        for group in self._shard_groups:
            if group and group[0].chr == "X":
                present = {s.sex for s in group}
                if self._default_chrX_sex not in present:
                    raise ValueError(
                        "default_chrX_sex must name a loaded chrX LD shard; "
                        f"got {self._default_chrX_sex!r}, present {sorted(present)!r}"
                    )

    def _multiply_r2_dense(self, M, chrX_sex=None):
        arr = np.asarray(M)
        if arr.ndim == 1:
            if arr.shape[0] != self._num_snp:
                raise ValueError(
                    f"multiply_r2: vector length mismatch: expected {self._num_snp}, got {arr.shape[0]}"
                )
            work = arr.reshape(-1, 1)
            vector_input = True
        elif arr.ndim == 2:
            if arr.shape[0] != self._num_snp:
                raise ValueError(
                    f"multiply_r2: matrix row count mismatch: expected {self._num_snp}, got {arr.shape[0]}"
                )
            work = arr
            vector_input = False
        else:
            raise ValueError("multiply_r2: M must be a vector or 2D matrix")

        if arr.dtype == np.float32:
            out_dtype = np.float32
        elif arr.dtype == np.float64:
            out_dtype = np.float64
        else:
            out_dtype = np.result_type(arr.dtype, np.float32)

        out = np.empty(work.shape, dtype=out_dtype)
        for shard, start, stop in self.iter_shards_with_offsets(chrX_sex):
            out[start:stop, :] = _multiply_ld_r2_dense(
                shard._ld_r2, work[start:stop, :], out_dtype
            )

        if vector_input:
            return out.reshape(arr.shape)
        return out

    def _multiply_r2_sparse(self, M, chrX_sex=None):
        if len(M.shape) != 2:
            raise ValueError("multiply_r2: sparse M must be a 2D matrix")
        if M.shape[0] != self._num_snp:
            raise ValueError(
                f"multiply_r2: matrix row count mismatch: expected {self._num_snp}, got {M.shape[0]}"
            )

        if M.dtype == np.float32:
            out_dtype = np.float32
        elif M.dtype == np.float64:
            out_dtype = np.float64
        else:
            out_dtype = np.result_type(M.dtype, np.float32)

        parts = []
        for shard, start, stop in self.iter_shards_with_offsets(chrX_sex):
            parts.append((shard._ld_r2 @ M[start:stop, :]).astype(out_dtype))
        return sparse.vstack(parts, format=M.getformat()).astype(out_dtype)


def load_ld(path, reference, default_chrX_sex=None) -> LDPanel:
    default_chrX_sex = validate_chrx_sex(
        "female" if default_chrX_sex is None else default_chrX_sex,
        "default_chrX_sex",
    )
    path = Path(path)
    ref_shards = list(reference.shards)
    if not ref_shards:
        raise ValueError("load_ld: reference must contain at least one shard")

    if path.is_dir():
        manifest = read_manifest(path / "ld_manifest.json", expected_runtime=PY_RUNTIME_FORMAT)
        groups = _load_panel_root(path, manifest, ref_shards)
    else:
        groups = _load_single_shard(path, ref_shards)

    return LDPanel(groups, default_chrX_sex=default_chrX_sex)


def validate_ld_distribution(path, check_payload_structure=False) -> dict:
    return validate_npz_distribution(path, check_payload_structure=check_payload_structure)


def fast_prune(logpvec, ld_panel: LDPanel, r2_threshold=0.2, chrX_sex=None) -> np.ndarray:
    threshold = float(r2_threshold)
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("fast_prune: r2_threshold must be a finite non-negative number")

    arr = np.asarray(logpvec)
    if arr.ndim > 2 or (arr.ndim == 2 and 1 not in arr.shape):
        raise ValueError("fast_prune: logpvec must be a vector")
    was_column = arr.ndim == 2 and arr.shape[1] == 1
    was_row = arr.ndim == 2 and arr.shape[0] == 1
    vec = arr.reshape(-1)
    if vec.size != ld_panel.num_snp:
        raise ValueError(
            f"fast_prune: vector length mismatch: expected {ld_panel.num_snp}, got {vec.size}"
        )

    out = vec.astype(np.result_type(vec.dtype, np.float32), copy=True)
    for shard, start, stop in ld_panel.iter_shards_with_offsets(chrX_sex):
        out[start:stop] = _fast_prune_shard(out[start:stop], shard._ld_r2, threshold)

    if was_column:
        return out.reshape(-1, 1)
    if was_row:
        return out.reshape(1, -1)
    return out


def _build_ld_r2(ld_r):
    data = np.square(ld_r.data).astype(ld_r.data.dtype, copy=False)
    return sparse.csc_matrix(
        (data, ld_r.indices, ld_r.indptr),
        shape=ld_r.shape,
        copy=False,
    )


def _multiply_ld_r2_dense(ld_r2, M, dtype):
    return np.asarray(ld_r2 @ M, dtype=dtype)


def _fast_prune_shard(values: np.ndarray, ld_r2, threshold: float) -> np.ndarray:
    out = values.copy()
    finite_idx = np.flatnonzero(np.isfinite(out))
    if finite_idx.size == 0:
        return out

    order = finite_idx[np.argsort(-np.abs(out[finite_idx]), kind="mergesort")]
    pruned = np.zeros(out.size, dtype=bool)
    retained = np.zeros(out.size, dtype=bool)

    for idx in order:
        idx = int(idx)
        if pruned[idx]:
            continue
        retained[idx] = True
        start = int(ld_r2.indptr[idx])
        stop = int(ld_r2.indptr[idx + 1])
        if start == stop:
            continue
        col_indices = ld_r2.indices[start:stop]
        col_data = ld_r2.data[start:stop]
        neighbors = col_indices[col_data >= threshold]
        neighbors = neighbors[(neighbors != idx) & (~retained[neighbors])]
        if neighbors.size:
            pruned[neighbors] = True
            out[neighbors] = np.nan
    return out


def _write_ld_npz_distribution(*args, **kwargs):
    from ._ld_writer import write_ld_npz_distribution

    return write_ld_npz_distribution(*args, **kwargs)


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
            validate_manifest_entry_agreement(entry, meta, shard_path)
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
    if len(ref_shards) != 1:
        raise ValueError(
            "single-shard LD loads require a single-shard reference; "
            f"reference has {len(ref_shards)} shards"
        )
    if ref_shards[0].label != shard.chr:
        raise ValueError(
            "single-shard LD chromosome does not match single-shard reference: "
            f"LD chr {shard.chr!r}, reference chr {ref_shards[0].label!r}"
        )
    _validate_reference_compatibility(shard, ref_shards[0], path)
    return [[shard]]


def _load_npz_shard(path: Path) -> tuple[LDShard, dict]:
    require_file(path)
    ld_r, a1freq, meta = read_npz_shard(path, check_payload_structure=False)
    shard = LDShard(
        meta["chr"],
        meta["sex"],
        int(meta["num_snp"]),
        ld_r,
        a1freq,
        meta["reference_checksum"],
        int(meta["num_monomorphic_snps"]),
    )
    if shard.num_monomorphic_snps:
        warnings.warn(
            f"{path}: LD shard contains {shard.num_monomorphic_snps} monomorphic SNPs; "
            "LD involving those SNPs is undefined and represented by omitted off-diagonal entries",
            RuntimeWarning,
            stacklevel=2,
        )
    return shard, meta


def _validate_reference_compatibility(shard: LDShard, ref_shard, path: Path) -> None:
    if shard.chr != ref_shard.label:
        raise ValueError(f"{path}: LD chr {shard.chr!r} does not match reference shard {ref_shard.label!r}")
    if shard.num_snp != ref_shard.num_snp:
        raise ValueError(f"{path}: LD num_snp does not match reference shard")
    if shard.reference_checksum != ref_shard.checksum:
        raise ValueError(f"{path}: LD reference_checksum does not match reference shard")
