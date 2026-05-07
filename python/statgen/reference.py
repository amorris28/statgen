import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from ._bfile_utils import _parse_bim
from ._bfile_utils import _validate_source_sort_order as _validate_reference_sort_order
from ._utils import CANONICAL_CHR_ORDER, allele_hash64, validate_requested_shards

logger = logging.getLogger(__name__)

_CACHE_SCHEMA = "reference_cache/0.1"


def _checksum_from_arrays(
    chr_arr: np.ndarray, bp_arr: np.ndarray, a1_arr: np.ndarray, a2_arr: np.ndarray
) -> str:
    if chr_arr.size == 0:
        return hashlib.md5(b"").hexdigest()

    line_series = (
        pd.Series(chr_arr, dtype="string")
        .str.cat(pd.Series(bp_arr, dtype=np.int64).astype("string"), sep=":")
        .str.cat(pd.Series(a1_arr, dtype="string"), sep=":")
        .str.cat(pd.Series(a2_arr, dtype="string"), sep=":")
    )
    text = line_series.str.cat(sep="\n") + "\n"
    return hashlib.md5(text.encode()).hexdigest()


class ReferenceShard:
    def __init__(
        self,
        label: str,
        chr_arr,
        snp_arr,
        bp_arr,
        a1_arr,
        a2_arr,
    ):
        self._label = label
        self._chr = np.asarray(chr_arr, dtype=object)
        self._snp = np.asarray(snp_arr, dtype=object)
        self._bp = np.asarray(bp_arr, dtype=np.int64)
        self._a1 = np.asarray(a1_arr, dtype=object)
        self._a2 = np.asarray(a2_arr, dtype=object)
        if not (
            self._chr.size
            == self._snp.size
            == self._bp.size
            == self._a1.size
            == self._a2.size
        ):
            raise ValueError("ReferenceShard vector lengths must match")
        if self._chr.size and not np.all(self._chr == self._label):
            raise ValueError(f"Reference shard {self._label} contains multiple chr labels")
        self._a1_hash64 = allele_hash64(self._a1)
        self._a2_hash64 = allele_hash64(self._a2)
        self._checksum = _checksum_from_arrays(self._chr, self._bp, self._a1, self._a2)

    @property
    def label(self) -> str:
        return self._label

    @property
    def num_snp(self) -> int:
        return len(self._chr)

    @property
    def chr(self) -> np.ndarray:
        return self._chr

    @property
    def snp(self) -> np.ndarray:
        return self._snp

    @property
    def bp(self) -> np.ndarray:
        return self._bp

    @property
    def a1(self) -> np.ndarray:
        return self._a1

    @property
    def a2(self) -> np.ndarray:
        return self._a2

    @property
    def a1_hash64(self) -> np.ndarray:
        return self._a1_hash64

    @property
    def a2_hash64(self) -> np.ndarray:
        return self._a2_hash64

    @property
    def checksum(self) -> str:
        return self._checksum

    @classmethod
    def _from_arrays(
        cls,
        label,
        chr_arr,
        snp_arr,
        bp_arr,
        a1_arr,
        a2_arr,
        checksum,
        a1_hash64=None,
        a2_hash64=None,
    ):
        obj = cls.__new__(cls)
        obj._label = label
        obj._chr = np.asarray(chr_arr, dtype=object)
        obj._snp = np.asarray(snp_arr, dtype=object)
        obj._bp = np.asarray(bp_arr, dtype=np.int64)
        obj._a1 = np.asarray(a1_arr, dtype=object)
        obj._a2 = np.asarray(a2_arr, dtype=object)
        if not (
            obj._chr.size
            == obj._snp.size
            == obj._bp.size
            == obj._a1.size
            == obj._a2.size
        ):
            raise ValueError("Invalid reference cache: panel-wide vector lengths mismatch")
        if obj._chr.size and not np.all(obj._chr == obj._label):
            raise ValueError(f"Invalid reference cache: shard {obj._label} contains multiple chr labels")
        obj._a1_hash64 = (
            allele_hash64(obj._a1)
            if a1_hash64 is None
            else np.asarray(a1_hash64, dtype=np.uint64).reshape(-1)
        )
        obj._a2_hash64 = (
            allele_hash64(obj._a2)
            if a2_hash64 is None
            else np.asarray(a2_hash64, dtype=np.uint64).reshape(-1)
        )
        if obj._a1_hash64.size != obj._bp.size or obj._a2_hash64.size != obj._bp.size:
            raise ValueError("Invalid reference cache: allele hash vector lengths mismatch")
        obj._checksum = checksum
        return obj


class ReferencePanel:
    def __init__(self, shards: list):
        self._shards = list(shards)
        offsets = []
        pos = 0
        for s in self._shards:
            offsets.append({"shard_label": s.label, "start0": pos, "stop0": pos + s.num_snp})
            pos += s.num_snp
        self._shard_offsets = offsets
        self._num_snp = pos

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def chr(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.chr for s in self._shards])

    @property
    def snp(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.snp for s in self._shards])

    @property
    def bp(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=np.int64)
        return np.concatenate([s.bp for s in self._shards])

    @property
    def a1(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.a1 for s in self._shards])

    @property
    def a2(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.a2 for s in self._shards])

    @property
    def a1_hash64(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=np.uint64)
        return np.concatenate([s.a1_hash64 for s in self._shards])

    @property
    def a2_hash64(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=np.uint64)
        return np.concatenate([s.a2_hash64 for s in self._shards])

    @property
    def shard_offsets(self) -> list:
        return list(self._shard_offsets)

    @property
    def shards(self) -> list:
        return list(self._shards)

    def select_shards(self, shards) -> "ReferencePanel":
        available = [s.label for s in self._shards]
        selected = validate_requested_shards(shards, available, "ReferencePanel.select_shards")
        by_label = {s.label: s for s in self._shards}
        return ReferencePanel([by_label[label] for label in selected])

    def validate_checksums(self) -> bool:
        for shard in self._shards:
            computed = _checksum_from_arrays(shard.chr, shard.bp, shard.a1, shard.a2)
            if computed != shard.checksum:
                raise ValueError(f"Reference checksum mismatch for shard {shard.label}")
        return True

    def is_object_compatible(self, obj) -> bool:
        log_fn = logger.warning

        obj_shards = getattr(obj, "shards", None)
        if obj_shards is None:
            log_fn("statgen: is_object_compatible: object has no shards attribute")
            return False

        obj_shards = list(obj_shards)
        if len(obj_shards) != len(self._shards):
            log_fn(
                "statgen: is_object_compatible: shard count mismatch: "
                f"reference has {len(self._shards)}, object has {len(obj_shards)}"
            )
            return False

        ok = True
        for ref_s, obj_s in zip(self._shards, obj_shards):
            obj_label = getattr(obj_s, "label", None)
            if obj_label is None:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: object shard has no label"
                )
                ok = False
                continue
            if obj_label != ref_s.label:
                log_fn(
                    "statgen: is_object_compatible: shard label mismatch: "
                    f"reference {ref_s.label}, object {obj_label}"
                )
                ok = False
                continue
            obj_num = getattr(obj_s, "num_snp", None)
            if obj_num is None:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: "
                    "object shard has no num_snp"
                )
                ok = False
                continue
            if obj_num != ref_s.num_snp:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: "
                    f"row count mismatch: reference {ref_s.num_snp}, object {obj_num}"
                )
                ok = False
                continue
            obj_chk = getattr(obj_s, "reference_checksum", None)
            if obj_chk is None:
                obj_chk = getattr(obj_s, "checksum", None)
            if obj_chk is not None and obj_chk != ref_s.checksum:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: reference_checksum mismatch"
                )
                ok = False
        return ok


def load_reference(path, shards=None) -> ReferencePanel:
    path_str = str(path)
    path_obj = Path(path_str)

    if "@" in path_str:
        available_labels = []
        available_paths = {}
        for label in CANONICAL_CHR_ORDER:
            candidate = Path(path_str.replace("@", label))
            if candidate.is_file():
                available_labels.append(label)
                available_paths[label] = candidate

        if not available_labels:
            raise FileNotFoundError(
                f"No BIM shards found matching template: {path_str}"
            )

        selected = validate_requested_shards(
            shards, available_labels, "load_reference"
        )
        shards = []
        for label in selected:
            bim_df = _parse_bim(available_paths[label])
            shards.append(
                ReferenceShard(
                    label,
                    bim_df["chr"].to_numpy(dtype=object),
                    bim_df["snp"].to_numpy(dtype=object),
                    bim_df["bp"].to_numpy(dtype=np.int64),
                    bim_df["a1"].to_numpy(dtype=object),
                    bim_df["a2"].to_numpy(dtype=object),
                )
            )
        return ReferencePanel(shards)

    bim_df = _parse_bim(path_obj)

    available_labels = [
        c for c in CANONICAL_CHR_ORDER if (bim_df["chr"] == c).any()
    ]
    selected = validate_requested_shards(shards, available_labels, "load_reference")
    out_shards = []
    for c in selected:
        shard_df = bim_df.loc[bim_df["chr"] == c]
        out_shards.append(
            ReferenceShard(
                c,
                shard_df["chr"].to_numpy(dtype=object),
                shard_df["snp"].to_numpy(dtype=object),
                shard_df["bp"].to_numpy(dtype=np.int64),
                shard_df["a1"].to_numpy(dtype=object),
                shard_df["a2"].to_numpy(dtype=object),
            )
        )
    return ReferencePanel(out_shards)


def save_reference_cache(panel: ReferencePanel, path, mode: str = "full") -> None:
    mode = str(mode).lower()
    if mode not in {"full", "thin"}:
        raise ValueError("save_reference_cache mode must be 'full' or 'thin'")

    # Metadata (schema, labels, checksums) as a compact JSON blob stored in the npz.
    # SNP-axis numeric vectors (bp) and string arrays are stored as native
    # binary numpy arrays — not JSON — per the performance contract.
    # Python keeps reference caches full even when mode='thin' is requested.
    meta = {
        "schema": _CACHE_SCHEMA,
        "mode": "full",
        "shard_labels": [s.label for s in panel.shards],
        "shard_checksums": [s.checksum for s in panel.shards],
    }
    arrays: dict = {
        "_meta": np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    }
    for i, s in enumerate(panel.shards):
        p = f"s{i}_"
        arrays[p + "chr"] = np.asarray(s.chr, dtype=str)
        arrays[p + "snp"] = np.asarray(s.snp, dtype=str)
        arrays[p + "bp"]  = s.bp
        arrays[p + "a1"]  = np.asarray(s.a1, dtype=str)
        arrays[p + "a2"]  = np.asarray(s.a2, dtype=str)
        arrays[p + "a1_hash64"] = s.a1_hash64
        arrays[p + "a2_hash64"] = s.a2_hash64
    np.savez_compressed(path, **arrays)


def load_reference_cache(path, shards=None) -> ReferencePanel:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(bytes(data["_meta"]).decode())
        schema = meta.get("schema")
        if schema != _CACHE_SCHEMA:
            raise ValueError(f"Unsupported reference cache schema: {schema!r}")
        mode = meta.get("mode")
        if mode is None:
            raise ValueError("reference cache missing mode; delete and rebuild old cache")
        if mode != "full":
            raise ValueError(f"Unsupported Python reference cache mode: {mode!r}")

        labels = list(meta.get("shard_labels", []))
        checksums = list(meta.get("shard_checksums", []))
        if len(labels) != len(checksums):
            raise ValueError("Invalid reference cache: shard_labels and shard_checksums length mismatch")
        selected = validate_requested_shards(shards, labels, "load_reference_cache")
        label_to_index = {label: i for i, label in enumerate(labels)}

        shard_objs = []
        for label in selected:
            i = label_to_index[label]
            p = f"s{i}_"
            if p + "a1_hash64" not in data or p + "a2_hash64" not in data:
                raise ValueError("reference cache missing a1_hash64/a2_hash64; rebuild cache")
            shard_objs.append(
                ReferenceShard._from_arrays(
                    label,
                    data[p + "chr"],
                    data[p + "snp"],
                    data[p + "bp"],
                    data[p + "a1"],
                    data[p + "a2"],
                    checksums[i],
                    data[p + "a1_hash64"],
                    data[p + "a2_hash64"],
                )
            )
    return ReferencePanel(shard_objs)
