from pathlib import Path

from ._ld_schema import require_file
from ._utils import validate_requested_shards
from .reference import ReferencePanel, load_reference


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


def load_bundled_reference_panel(root, manifest: dict, shards, *, where: str) -> ReferencePanel:
    root = Path(root)
    available = _manifest_chr_labels(manifest)
    selected = validate_requested_shards(shards, available, where)
    entries_by_chr = {}
    for entry in manifest["shards"]:
        entries_by_chr.setdefault(entry["chr"], []).append(entry)

    ref_shards = []
    for label in selected:
        entries = entries_by_chr[label]
        reference_path = root / entries[0]["reference_bim"]
        ref_shard = validate_bundled_reference_bim(reference_path, entries[0], target="manifest")
        for entry in entries[1:]:
            validate_reference_shard_for_ld_entry(ref_shard, entry, reference_path, target="manifest")
        ref_shards.append(ref_shard)
    return ReferencePanel(ref_shards)


def _manifest_chr_labels(manifest: dict) -> list[str]:
    labels = []
    seen = set()
    for entry in manifest["shards"]:
        label = entry["chr"]
        if label not in seen:
            labels.append(label)
            seen.add(label)
    return labels
