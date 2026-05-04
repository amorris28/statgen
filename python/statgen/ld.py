from pathlib import Path

import numpy as np

from ._ld_npz import read_npz_shard, validate_npz_distribution
from ._ld_schema import (
    PY_RUNTIME_FORMAT,
    read_manifest,
    require_file,
    validate_chrx_sex,
    validate_manifest_entry_agreement,
)


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
        self._default_chrX_sex = validate_chrx_sex(default_chrX_sex, "default_chrX_sex")
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
    matching = [s for s in ref_shards if s.label == shard.chr]
    if len(ref_shards) != 1 or len(matching) != 1:
        raise ValueError(
            "single-shard LD loads require a single-shard reference with the same chromosome"
        )
    _validate_reference_compatibility(shard, matching[0], path)
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
    )
    return shard, meta


def _validate_reference_compatibility(shard: LDShard, ref_shard, path: Path) -> None:
    if shard.chr != ref_shard.label:
        raise ValueError(f"{path}: LD chr {shard.chr!r} does not match reference shard {ref_shard.label!r}")
    if shard.num_snp != ref_shard.num_snp:
        raise ValueError(f"{path}: LD num_snp does not match reference shard")
    if shard.reference_checksum != ref_shard.checksum:
        raise ValueError(f"{path}: LD reference_checksum does not match reference shard")
