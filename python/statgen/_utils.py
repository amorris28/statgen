import warnings

import numpy as np

CANONICAL_CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
CHR_RANK = {c: i for i, c in enumerate(CANONICAL_CHR_ORDER)}
_ALLELE_HASH_P = 2_147_483_647
_ALLELE_HASH_BASE1 = 257
_ALLELE_HASH_BASE2 = 263
_ALLELE_HASH_MAX_CHARS = 150


def validate_requested_shards(shards, available_labels, where: str) -> list[str]:
    labels = list(available_labels)
    if shards is None:
        return labels
    if isinstance(shards, str):
        raise ValueError(f"{where}: shards must be a non-empty list of unique canonical contig labels")

    requested = list(shards)
    if not requested:
        raise ValueError(f"{where}: shards must be a non-empty list of unique canonical contig labels")

    seen = set()
    prev_rank = None
    for label in requested:
        if label not in CHR_RANK:
            raise ValueError(f"{where}: unsupported shard label {label!r}; expected canonical labels 1-22 or X")
        if label in seen:
            raise ValueError(f"{where}: duplicate shard label {label!r} in shards")
        rank = CHR_RANK[label]
        if prev_rank is not None and rank <= prev_rank:
            raise ValueError(f"{where}: shards must be in canonical subsequence order")
        seen.add(label)
        prev_rank = rank
        if label not in labels:
            raise ValueError(f"{where}: requested shard {label!r} is not present")
    return requested


def allele_hash64(alleles) -> np.ndarray:
    arr = np.asarray(alleles, dtype=object).reshape(-1)
    n = arr.size
    if n == 0:
        return np.array([], dtype=np.uint64)

    strings = [str(allele) for allele in arr]
    lengths = np.fromiter((len(s) for s in strings), dtype=np.int64, count=n)
    if np.any(lengths > _ALLELE_HASH_MAX_CHARS):
        warnings.warn(
            "Allele length exceeds 150 characters; hashing uses first 150 characters",
            RuntimeWarning,
            stacklevel=2,
        )
        strings = [s[:_ALLELE_HASH_MAX_CHARS] for s in strings]
        lengths = np.minimum(lengths, _ALLELE_HASH_MAX_CHARS)

    encoded = [s.encode("utf-8") for s in strings]
    byte_lengths = np.fromiter((len(b) for b in encoded), dtype=np.int64, count=n)
    max_len = int(byte_lengths.max(initial=0))
    mat = np.zeros((n, max_len), dtype=np.uint8)
    if max_len:
        flat = np.frombuffer(b"".join(encoded), dtype=np.uint8)
        row_idx = np.repeat(np.arange(n, dtype=np.int64), byte_lengths)
        starts = np.cumsum(byte_lengths, dtype=np.int64) - byte_lengths
        col_idx = np.arange(flat.size, dtype=np.int64) - np.repeat(starts, byte_lengths)
        mat[row_idx, col_idx] = flat

    h1 = np.ones(n, dtype=np.int64)
    h2 = np.ones(n, dtype=np.int64)
    for j in range(max_len):
        active = j < byte_lengths
        x = mat[active, j].astype(np.int64) + 1
        h1[active] = (h1[active] * _ALLELE_HASH_BASE1 + x) % _ALLELE_HASH_P
        h2[active] = (h2[active] * _ALLELE_HASH_BASE2 + x) % _ALLELE_HASH_P

    return (h1.astype(np.uint64) << np.uint64(32)) + h2.astype(np.uint64)
