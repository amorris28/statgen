from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

IGNORED_CHR = {"Y", "MT"}


def merge_intervals(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    if starts.size == 0:
        return np.empty((0, 2), dtype=np.int64)

    order = np.lexsort((ends, starts))
    s = starts[order]
    e = ends[order]

    running_max_end = np.maximum.accumulate(e)
    new_group = np.empty(len(s), dtype=bool)
    new_group[0] = True
    new_group[1:] = s[1:] > running_max_end[:-1]

    group_starts_idx = np.where(new_group)[0]
    merged_starts = s[group_starts_idx]
    merged_ends = np.maximum.reduceat(e, group_starts_idx)

    return np.column_stack([merged_starts.astype(np.int64), merged_ends.astype(np.int64)])


def binary_intervals_by_chr_from_arrays(
    chr_values: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
) -> dict[str, np.ndarray]:
    intervals_by_chr: dict[str, np.ndarray] = {}
    grouped = pd.DataFrame({"chr": chr_values, "start": starts, "end": ends}).groupby(
        "chr", sort=False
    )
    for chr_label, grp in grouped:
        if chr_label in IGNORED_CHR:
            continue
        intervals_by_chr[str(chr_label)] = merge_intervals(
            grp["start"].to_numpy(dtype=np.int64),
            grp["end"].to_numpy(dtype=np.int64),
        )
    return intervals_by_chr


def validate_numeric_non_overlapping(
    chr_values: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    path: Path,
    row_base0: int = 0,
) -> None:
    source = pd.DataFrame(
        {
            "chr": chr_values,
            "start": starts,
            "end": ends,
            "row0": np.arange(starts.size, dtype=np.int64),
        }
    )
    for chr_label, grp in source.groupby("chr", sort=False):
        if chr_label in IGNORED_CHR or grp.shape[0] <= 1:
            continue
        order = np.lexsort((grp["end"].to_numpy(), grp["start"].to_numpy()))
        s = grp["start"].to_numpy(dtype=np.int64)[order]
        e = grp["end"].to_numpy(dtype=np.int64)[order]
        rows = grp["row0"].to_numpy(dtype=np.int64)[order]
        prev_max = np.maximum.accumulate(e)[:-1]
        overlap = s[1:] < prev_max
        if overlap.any():
            j = int(np.flatnonzero(overlap)[0]) + 1
            raise ValueError(
                f"{path}: row {row_base0 + int(rows[j]) + 1}: "
                f"numeric annotation intervals overlap on chromosome {chr_label}"
            )


def paint_mask(bp: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    n = bp.size
    if intervals.size == 0:
        return np.zeros(n, dtype=np.float64)

    pos0 = np.asarray(bp, dtype=np.int64) - 1
    starts = intervals[:, 0]
    ends = intervals[:, 1]

    idx = np.searchsorted(starts, pos0, side="right") - 1
    valid = idx >= 0
    out = np.zeros(n, dtype=np.float64)
    if valid.any():
        valid_idx = idx[valid]
        inside = pos0[valid] < ends[valid_idx]
        out[np.flatnonzero(valid)[inside]] = 1.0
    return out


def paint_numeric_values(
    bp: np.ndarray,
    intervals: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    n = bp.size
    k = values.shape[1]
    out = np.zeros((n, k), dtype=np.float64)
    if intervals.size == 0:
        return out

    pos0 = np.asarray(bp, dtype=np.int64) - 1
    starts = intervals[:, 0]
    ends = intervals[:, 1]
    idx = np.searchsorted(starts, pos0, side="right") - 1
    valid = idx >= 0
    if valid.any():
        valid_idx = idx[valid]
        inside = pos0[valid] < ends[valid_idx]
        out[np.flatnonzero(valid)[inside], :] = values[valid_idx[inside], :]
    return out


def paint_binary_column(intervals_by_chr: dict[str, np.ndarray], reference) -> sparse.csr_matrix:
    n = int(reference.num_snp)
    mask = np.zeros(n, dtype=np.float64)
    for ref_shard, off in zip(reference.shards, reference.shard_offsets):
        start = int(off["start0"])
        stop = int(off["stop0"])
        intervals = intervals_by_chr.get(ref_shard.label)
        if intervals is None:
            continue
        mask[start:stop] = paint_mask(ref_shard.bp, intervals)
    return sparse.csr_matrix(mask.reshape(-1, 1))
