import numpy as np

_KEY_DTYPE = np.dtype([("bp", "<i8"), ("a1_hash64", "<u8"), ("a2_hash64", "<u8")])


def numeric_variant_keys(bp_arr, a1_hash64, a2_hash64) -> np.ndarray:
    bp = np.asarray(bp_arr, dtype=np.int64).reshape(-1)
    a1h = np.asarray(a1_hash64, dtype=np.uint64).reshape(-1)
    a2h = np.asarray(a2_hash64, dtype=np.uint64).reshape(-1)
    if not (bp.size == a1h.size == a2h.size):
        raise ValueError("numeric variant key vector lengths mismatch")
    keys = np.empty(bp.size, dtype=_KEY_DTYPE)
    keys["bp"] = bp
    keys["a1_hash64"] = a1h
    keys["a2_hash64"] = a2h
    return keys


def check_unique_sorted_variant_keys(sorted_keys: np.ndarray, where: str, object_name: str = "source") -> None:
    if sorted_keys.size < 2:
        return
    dup = sorted_keys[1:] == sorted_keys[:-1]
    if np.any(dup):
        pos = int(np.flatnonzero(dup)[0] + 1)
        key = sorted_keys[pos]
        raise ValueError(
            f"Ambiguous duplicate {object_name}/reference matching key in {where}: "
            f"bp={int(key['bp'])}, a1_hash64={int(key['a1_hash64'])}, "
            f"a2_hash64={int(key['a2_hash64'])}"
        )


def match_shard_numeric(
    ref_keys: np.ndarray,
    src_keys: np.ndarray,
    shard_label: str,
    *,
    source_name: str = "source",
) -> np.ndarray:
    ref_order = np.argsort(ref_keys, kind="mergesort", order=("bp", "a1_hash64", "a2_hash64"))
    ref_sorted = ref_keys[ref_order]
    check_unique_sorted_variant_keys(ref_sorted, f"reference shard {shard_label}", source_name)

    if src_keys.size == 0:
        return np.full(ref_keys.size, -1, dtype=np.int64)

    src_order = np.argsort(src_keys, kind="mergesort", order=("bp", "a1_hash64", "a2_hash64"))
    src_sorted = src_keys[src_order]
    check_unique_sorted_variant_keys(src_sorted, f"{source_name} shard {shard_label}", source_name)

    pos = np.searchsorted(src_sorted, ref_keys)
    in_range = pos < src_sorted.size
    matched = np.zeros(ref_keys.size, dtype=bool)
    matched[in_range] = src_sorted[pos[in_range]] == ref_keys[in_range]

    out = np.full(ref_keys.size, -1, dtype=np.int64)
    out[matched] = src_order[pos[matched]]
    return out
