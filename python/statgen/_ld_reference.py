from pathlib import Path

from ._ld_schema import require_file
from .reference import load_reference


def validate_bundled_reference_bim(path, entry: dict, *, target: str):
    path = Path(path)
    require_file(path)
    panel = load_reference(path)
    ref_shards = panel.shards
    if len(ref_shards) != 1:
        raise ValueError(f"{path}: bundled reference_bim must contain exactly one shard")
    ref = ref_shards[0]
    validate_reference_shard_for_ld_entry(ref, entry, path, target=target)
    return ref


def validate_reference_shard_for_ld_entry(ref_shard, entry: dict, path, *, target: str) -> None:
    path = Path(path)
    if ref_shard.label != entry["chr"]:
        raise ValueError(f"{path}: bundled reference_bim chr does not match {target}")
    if ref_shard.num_snp != int(entry["num_snp"]):
        raise ValueError(f"{path}: bundled reference_bim num_snp does not match {target}")
    if ref_shard.checksum != entry["reference_checksum"]:
        raise ValueError(f"{path}: bundled reference_bim reference_checksum does not match {target}")
