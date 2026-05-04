#!/usr/bin/env python3
"""Build statgen Python LD distributions from PLINK bfiles."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from statgen._ld_npz import validate_npz_distribution
from statgen._ld_writer import write_ld_npz_shards
from statgen._utils import CANONICAL_CHR_ORDER
from statgen.reference import load_reference


DEFAULT_LD_WINDOW_KB = 10000
DEFAULT_LD_R2_THRESHOLD = 0.05
VCOR_PARSE_CHUNK_ROWS = 1_000_000
VARIANT_KEY_DTYPE = np.dtype([("chr", object), ("bp", np.int64), ("a1", object), ("a2", object)])
VARIANT_KEY_ORDER = ("chr", "bp", "a1", "a2")


def main(argv=None) -> int:
    args = _parse_args(argv)
    build_ld_distribution(
        bfile=args.bfile,
        out=args.out,
        shard=args.shard,
        plink2=args.plink2,
        ld_window_kb=args.ld_window_kb,
        ld_r2_threshold=args.ld_r2_threshold,
        no_sex_split=args.no_sex_split,
        scratch=args.scratch,
        mind=args.mind[0] if args.mind is not None else None,
        mind_mode=args.mind[1] if args.mind is not None and len(args.mind) == 2 else None,
        keep_samples=args.keep_samples,
        remove_samples=args.remove_samples,
        keep_families=args.keep_families,
        remove_families=args.remove_families,
        threads=args.threads,
        memory=args.memory,
    )
    return 0


def build_ld_distribution(
    *,
    bfile,
    out,
    shard,
    plink2="plink2",
    ld_window_kb=DEFAULT_LD_WINDOW_KB,
    ld_r2_threshold=DEFAULT_LD_R2_THRESHOLD,
    no_sex_split=False,
    scratch=None,
    mind=None,
    mind_mode=None,
    keep_samples=None,
    remove_samples=None,
    keep_families=None,
    remove_families=None,
    threads=None,
    memory=None,
) -> None:
    shards = _discover_bfile_shards(bfile, shard)
    if not shards:
        raise ValueError(f"No bfile shard {shard!r} found for {bfile}")

    plink_options = _normalize_plink_options(
        mind=mind,
        mind_mode=mind_mode,
        threads=threads,
        memory=memory,
    )
    sample_filters = {
        "keep_samples": keep_samples,
        "remove_samples": remove_samples,
        "keep_families": keep_families,
        "remove_families": remove_families,
    }
    has_sample_filters = any(value is not None for value in sample_filters.values())

    plink_version = _plink_version(plink2)
    tmp_root = _prepare_scratch_dir(out, scratch)
    shard_spec = shards[0]
    shard_records = write_ld_npz_shards(
        out,
        _build_shard_specs(
            shard_spec,
            tmp_root,
            plink2=plink2,
            ld_window_kb=int(ld_window_kb),
            ld_r2_threshold=float(ld_r2_threshold),
            no_sex_split=bool(no_sex_split),
            plink_version=plink_version,
            plink_options=plink_options,
            sample_filters=sample_filters,
            has_sample_filters=has_sample_filters,
        ),
    )
    for record in shard_records:
        validate_npz_distribution(Path(out) / record["file"], check_payload_structure=True)

    return None


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile", required=True, help="PLINK bfile prefix, or sharded prefix with @ placeholder")
    parser.add_argument("--out", required=True, help="Output LD distribution directory")
    parser.add_argument("--shard", required=True, choices=CANONICAL_CHR_ORDER, help="Canonical shard label to build, one of 1-22 or X")
    parser.add_argument("--plink2", default="plink2", help="PLINK2 executable path")
    parser.add_argument("--ld-window-kb", type=int, default=DEFAULT_LD_WINDOW_KB)
    parser.add_argument("--ld-r2-threshold", type=float, default=DEFAULT_LD_R2_THRESHOLD)
    parser.add_argument(
        "--scratch",
        default=None,
        help="Directory for PLINK2 intermediate outputs; defaults to <out>/_statgen_build_ld_work",
    )
    parser.add_argument("--no-sex-split", action="store_true", help="Build one combined chrX shard instead of female/male shards")
    parser.add_argument("--mind", nargs="+", metavar="VALUE", help="PLINK2 sample missingness filter: --mind <x> [dosage|hh-missing]")
    parser.add_argument("--keep-samples", default=None, help="Sample ID file to keep, intersected with statgen chrX keep files")
    parser.add_argument("--remove-samples", default=None, help="Sample ID file to remove before building LD")
    parser.add_argument("--keep-families", default=None, help="Family ID file to keep before building LD")
    parser.add_argument("--remove-families", default=None, help="Family ID file to remove before building LD")
    parser.add_argument("--threads", type=int, default=None, help="PLINK2 --threads value")
    parser.add_argument("--memory", type=int, default=None, help="PLINK2 --memory value in MB")
    args = parser.parse_args(argv)
    if args.mind is not None:
        if len(args.mind) not in {1, 2}:
            parser.error("--mind expects <x> optionally followed by dosage or hh-missing")
        if len(args.mind) == 2 and args.mind[1] not in {"dosage", "hh-missing"}:
            parser.error("--mind mode must be dosage or hh-missing")
    return args


def _prepare_scratch_dir(out, scratch) -> Path:
    root = Path(scratch) if scratch is not None else Path(out) / "_statgen_build_ld_work"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _discover_bfile_shards(bfile, shard) -> list[dict]:
    if shard not in CANONICAL_CHR_ORDER:
        raise ValueError("shard must be one of 1-22 or X")
    bfile = str(bfile)
    if "@" in bfile:
        prefix = Path(bfile.replace("@", shard))
        _require_bfile(prefix)
        reference = load_reference(_bfile_component(prefix, ".bim"))
        if len(reference.shards) != 1 or reference.shards[0].label != shard:
            raise ValueError(f"{prefix}.bim must contain exactly canonical shard {shard}")
        return [{"chr": shard, "bfile": prefix, "reference_shard": reference.shards[0]}]

    prefix = Path(bfile)
    _require_bfile(prefix)
    reference = load_reference(_bfile_component(prefix, ".bim"))
    for reference_shard in reference.shards:
        if reference_shard.label == shard:
            return [{"chr": shard, "bfile": prefix, "reference_shard": reference_shard}]
    raise ValueError(f"{prefix}.bim does not contain shard {shard}")


def _require_bfile(prefix: Path) -> None:
    for suffix in (".bed", ".bim", ".fam"):
        path = _bfile_component(prefix, suffix)
        if not path.is_file():
            raise FileNotFoundError(f"missing PLINK bfile component: {path}")


def _bfile_component(prefix: Path, suffix: str) -> Path:
    return Path(str(prefix) + suffix)


def _build_shard_specs(
    shard,
    tmp_root: Path,
    *,
    plink2,
    ld_window_kb,
    ld_r2_threshold,
    no_sex_split,
    plink_version,
    plink_options,
    sample_filters,
    has_sample_filters,
) -> list[dict]:
    chr_label = shard["chr"]
    fam = _read_fam(_bfile_component(shard["bfile"], ".fam"))
    filtered_fam = _apply_user_sample_filters(fam, **sample_filters)

    if chr_label == "X" and not no_sex_split:
        _require_nonmissing_chrx_sex(filtered_fam, _bfile_component(shard["bfile"], ".fam"))
        jobs = [
            ("female", filtered_fam.loc[filtered_fam["sex"] == "2", ["fid", "iid"]]),
            ("male", filtered_fam.loc[filtered_fam["sex"] == "1", ["fid", "iid"]]),
        ]
    elif chr_label == "X":
        keep_samples = filtered_fam.loc[:, ["fid", "iid"]] if has_sample_filters else None
        jobs = [("combined", keep_samples)]
    else:
        keep_samples = filtered_fam.loc[:, ["fid", "iid"]] if has_sample_filters else None
        jobs = [(None, keep_samples)]

    ref_key_index = _reference_key_index(shard["reference_shard"])
    specs = []
    for sex, keep_samples in jobs:
        if keep_samples is not None and keep_samples.empty and sex is not None:
            raise ValueError(f"chrX {sex} build has no samples")
        if keep_samples is not None and keep_samples.empty:
            raise ValueError(f"chr{chr_label} build has no samples")

        tag = f"chr{chr_label}" if sex is None else f"chr{chr_label}_{sex}"
        out_prefix = tmp_root / tag
        keep_file = None
        num_sample = int(len(fam) if keep_samples is None else len(keep_samples))
        if keep_samples is not None:
            keep_file = tmp_root / f"{tag}.keep"
            keep_samples.to_csv(keep_file, sep="\t", header=False, index=False)

        freq_cmd = _plink_freq_command(
            plink2,
            shard["bfile"],
            out_prefix,
            chr_label=chr_label,
            keep_file=keep_file,
            plink_options=plink_options,
        )
        _run(freq_cmd)
        a1freq = _parse_afreq(out_prefix.with_suffix(".afreq"), shard["reference_shard"], ref_key_index)

        ld_cmd = _plink_ld_command(
            plink2,
            shard["bfile"],
            out_prefix,
            chr_label=chr_label,
            keep_file=keep_file,
            ld_window_kb=ld_window_kb,
            ld_r2_threshold=ld_r2_threshold,
            plink_options=plink_options,
        )
        _run(ld_cmd)
        ld_pairs = _parse_vcor(out_prefix.with_suffix(".vcor"), shard["reference_shard"], ref_key_index=ref_key_index)

        specs.append(
            {
                "reference_shard": shard["reference_shard"],
                "sex": sex,
                "a1freq": a1freq,
                "ld_pairs": ld_pairs,
                "build_metadata": {
                    "build_tool": "plink2",
                    "build_command": _format_command(ld_cmd),
                    "plink_version": plink_version,
                    "ld_window_kb": int(ld_window_kb),
                    "ld_r2_threshold": float(ld_r2_threshold),
                    "num_sample": num_sample,
                },
            }
        )
    return specs


def _plink_freq_command(plink2, bfile, out_prefix, *, chr_label, keep_file, plink_options):
    cmd = [
        str(plink2),
        "--bfile",
        str(bfile),
        "--chr",
        str(chr_label),
        "--freq",
        "cols=chrom,pos,ref,alt1,alt1freq,nobs",
        "--out",
        str(out_prefix),
        "--keep-allele-order",
    ]
    _append_plink_filter_and_resource_args(cmd, keep_file=keep_file, plink_options=plink_options)
    return cmd


def _plink_ld_command(
    plink2,
    bfile,
    out_prefix,
    *,
    chr_label,
    keep_file,
    ld_window_kb,
    ld_r2_threshold,
    plink_options,
):
    cmd = [
        str(plink2),
        "--bfile",
        str(bfile),
        "--chr",
        str(chr_label),
        "--r-unphased",
        "ref-based",
        "cols=chrom,pos,ref,alt1",
        "--ld-window-kb",
        str(int(ld_window_kb)),
        "--ld-window-r2",
        str(float(ld_r2_threshold)),
        "--out",
        str(out_prefix),
        "--keep-allele-order",
    ]
    _append_plink_filter_and_resource_args(cmd, keep_file=keep_file, plink_options=plink_options)
    return cmd


def _append_plink_filter_and_resource_args(cmd: list[str], *, keep_file, plink_options: dict) -> None:
    if keep_file is not None:
        cmd.extend(["--keep", str(keep_file)])
    if plink_options["mind"] is not None:
        cmd.extend(["--mind", plink_options["mind"]])
        if plink_options["mind_mode"] is not None:
            cmd.append(plink_options["mind_mode"])
    if plink_options["threads"] is not None:
        cmd.extend(["--threads", str(plink_options["threads"])])
    if plink_options["memory"] is not None:
        cmd.extend(["--memory", str(plink_options["memory"])])


def _normalize_plink_options(*, mind, mind_mode, threads, memory) -> dict:
    if mind is None and mind_mode is not None:
        raise ValueError("mind_mode requires mind")
    if mind is not None:
        mind_value = float(mind)
        if not np.isfinite(mind_value) or mind_value < 0 or mind_value > 1:
            raise ValueError("mind must be a finite value in [0, 1]")
        mind = str(mind)
    if mind_mode is not None and mind_mode not in {"dosage", "hh-missing"}:
        raise ValueError("mind_mode must be 'dosage' or 'hh-missing'")
    if threads is not None:
        threads = int(threads)
        if threads < 1:
            raise ValueError("threads must be a positive integer")
    if memory is not None:
        memory = int(memory)
        if memory < 1:
            raise ValueError("memory must be a positive integer")
    return {
        "mind": mind,
        "mind_mode": mind_mode,
        "threads": threads,
        "memory": memory,
    }


def _run(cmd):
    try:
        subprocess.run(cmd, check=True, text=True, capture_output=True)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"PLINK2 executable not found: {cmd[0]}") from exc
    except subprocess.CalledProcessError as exc:
        message = exc.stderr.strip() or exc.stdout.strip() or f"exit code {exc.returncode}"
        raise RuntimeError(f"PLINK2 command failed: {_format_command(cmd)}\n{message}") from exc


def _plink_version(plink2):
    try:
        result = subprocess.run([str(plink2), "--version"], check=True, text=True, capture_output=True)
    except Exception:
        return None
    text = (result.stdout or result.stderr).strip()
    return text.splitlines()[0] if text else None


def _read_fam(path: Path) -> pd.DataFrame:
    rows = []
    with open(path, newline="") as f:
        for line_no, raw in enumerate(f, start=1):
            fields = raw.rstrip("\n\r").split()
            if len(fields) != 6:
                raise ValueError(f"{path}:{line_no}: FAM must have 6 columns")
            rows.append(fields)
    if not rows:
        raise ValueError(f"{path}: FAM must contain at least one sample")
    return pd.DataFrame(rows, columns=["fid", "iid", "father", "mother", "sex", "phenotype"])


def _apply_user_sample_filters(
    fam: pd.DataFrame,
    *,
    keep_samples,
    remove_samples,
    keep_families,
    remove_families,
) -> pd.DataFrame:
    mask = np.ones(len(fam), dtype=bool)
    if keep_families is not None:
        mask &= fam["fid"].isin(_read_family_filter_file(keep_families)).to_numpy()
    if remove_families is not None:
        mask &= ~fam["fid"].isin(_read_family_filter_file(remove_families)).to_numpy()
    if keep_samples is not None:
        mask &= _sample_filter_mask(fam, keep_samples)
    if remove_samples is not None:
        mask &= ~_sample_filter_mask(fam, remove_samples)
    return fam.loc[mask].reset_index(drop=True)


def _read_family_filter_file(path) -> set[str]:
    path = Path(path)
    families = set()
    saw_header = False
    with open(path, newline="") as f:
        for raw in f:
            fields = raw.rstrip("\n\r").split()
            if not fields:
                continue
            if not saw_header and fields[0].startswith("#FID"):
                saw_header = True
                continue
            saw_header = True
            families.add(fields[0])
    return families


def _sample_filter_mask(fam: pd.DataFrame, path) -> np.ndarray:
    sample_filter = _read_sample_filter_file(path)
    if sample_filter["iid_only"]:
        return fam["iid"].isin(sample_filter["iid"]).to_numpy()
    selected = sample_filter["pairs"]
    return np.fromiter(
        (((fid, iid) in selected) for fid, iid in zip(fam["fid"], fam["iid"])),
        dtype=bool,
        count=len(fam),
    )


def _read_sample_filter_file(path) -> dict:
    path = Path(path)
    rows = []
    header = None
    with open(path, newline="") as f:
        for raw in f:
            fields = raw.rstrip("\n\r").split()
            if not fields:
                continue
            if header is None and fields[0].startswith("#"):
                header = [field.lstrip("#") for field in fields]
                continue
            rows.append(fields)

    if header is None:
        iid_only = all(len(row) == 1 for row in rows)
        if iid_only:
            return {"iid_only": True, "iid": {row[0] for row in rows}, "pairs": set()}
        bad = [i for i, row in enumerate(rows, start=1) if len(row) < 2]
        if bad:
            raise ValueError(f"{path}:{bad[0]}: sample filter rows must have one IID column or at least FID IID columns")
        return {"iid_only": False, "iid": set(), "pairs": {(row[0], row[1]) for row in rows if len(row) >= 2}}

    header_index = {name: i for i, name in enumerate(header)}
    if "IID" not in header_index:
        raise ValueError(f"{path}: sample filter header must contain IID")
    if "FID" not in header_index:
        iid_i = header_index["IID"]
        _require_sample_filter_width(path, rows, iid_i)
        return {"iid_only": True, "iid": {row[iid_i] for row in rows}, "pairs": set()}

    fid_i = header_index["FID"]
    iid_i = header_index["IID"]
    _require_sample_filter_width(path, rows, max(fid_i, iid_i))
    pairs = {(row[fid_i], row[iid_i]) for row in rows}
    return {"iid_only": False, "iid": set(), "pairs": pairs}


def _require_sample_filter_width(path: Path, rows: list[list[str]], required_index: int) -> None:
    for i, row in enumerate(rows, start=2):
        if len(row) <= required_index:
            raise ValueError(f"{path}:{i}: sample filter row is shorter than its header")


def _require_nonmissing_chrx_sex(fam: pd.DataFrame, path: Path) -> None:
    bad = ~fam["sex"].isin(["1", "2"])
    if bad.any():
        idx = int(np.flatnonzero(bad.to_numpy())[0])
        raise ValueError(f"{path}:{idx + 1}: chrX sex-specific builds require FAM sex codes 1 or 2")


def _parse_afreq(path: Path, reference_shard, ref_key_index=None) -> np.ndarray:
    header, chunks = _iter_table_chunks_with_hash_header(path, chunk_rows=None)
    required = {"CHROM", "POS", "REF", "ALT1", "ALT1_FREQ"}
    if not required.issubset(header):
        missing = ", ".join(sorted(required.difference(header)))
        raise ValueError(f"{path}: missing required .afreq columns: {missing}")

    ref_keys, ref_idx = ref_key_index if ref_key_index is not None else _reference_key_index(reference_shard)
    out = np.full(reference_shard.num_snp, np.nan, dtype=np.float32)
    seen = np.zeros(reference_shard.num_snp, dtype=bool)
    for chunk in chunks:
        keys = _variant_key_array(chunk, "CHROM", "POS", "ALT1", "REF")
        idx, missing = _map_variant_keys_to_reference(keys, ref_keys, ref_idx)
        if missing.any():
            raise ValueError(f"{path}: frequency row references variant absent from reference shard")
        if np.unique(idx).size != idx.size or seen[idx].any():
            raise ValueError(f"{path}: duplicate frequency row")
        freq = pd.to_numeric(chunk["ALT1_FREQ"], errors="coerce").to_numpy(dtype=np.float32)
        out[idx] = freq
        seen[idx] = True
    if not seen.all():
        raise ValueError(f"{path}: missing frequency row for reference variant")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{path}: frequency values must be finite")
    return out


def _parse_vcor(path: Path, reference_shard, chunk_rows=VCOR_PARSE_CHUNK_ROWS, ref_key_index=None) -> pd.DataFrame:
    header, chunks = _iter_table_chunks_with_hash_header(path, chunk_rows)
    required = {"CHROM_A", "POS_A", "REF_A", "ALT1_A", "CHROM_B", "POS_B", "REF_B", "ALT1_B"}
    if not required.issubset(header):
        missing = ", ".join(sorted(required.difference(header)))
        raise ValueError(f"{path}: missing required .vcor columns: {missing}")
    r_col = _vcor_r_column(header, path)

    ref_keys, ref_idx = ref_key_index if ref_key_index is not None else _reference_key_index(reference_shard)
    idx1_parts = []
    idx2_parts = []
    r_parts = []
    pair_key_parts = []
    for chunk in chunks:
        key_a = _variant_key_array(chunk, "CHROM_A", "POS_A", "ALT1_A", "REF_A")
        key_b = _variant_key_array(chunk, "CHROM_B", "POS_B", "ALT1_B", "REF_B")
        arr_a, missing_a = _map_variant_keys_to_reference(key_a, ref_keys, ref_idx)
        arr_b, missing_b = _map_variant_keys_to_reference(key_b, ref_keys, ref_idx)
        missing = missing_a | missing_b
        if missing.any():
            raise ValueError(f"{path}: LD row references variant absent from reference shard")

        keep = arr_a != arr_b
        if not keep.any():
            continue
        arr_a = arr_a[keep]
        arr_b = arr_b[keep]
        r = pd.to_numeric(chunk.loc[keep, r_col], errors="coerce").to_numpy(dtype=np.float32)
        if not np.all(np.isfinite(r)):
            raise ValueError(f"{path}: signed LD r values must be finite")

        lo = np.minimum(arr_a, arr_b)
        hi = np.maximum(arr_a, arr_b)
        pair_keys = lo.astype(np.int64) * np.int64(reference_shard.num_snp) + hi.astype(np.int64)

        idx1_parts.append(arr_a)
        idx2_parts.append(arr_b)
        r_parts.append(r)
        pair_key_parts.append(pair_keys)

    if not idx1_parts:
        return pd.DataFrame(
            {
                "idx1": np.array([], dtype=np.int32),
                "idx2": np.array([], dtype=np.int32),
                "r": np.array([], dtype=np.float32),
            }
        )
    pair_keys = np.concatenate(pair_key_parts)
    if np.unique(pair_keys).size != pair_keys.size:
        raise ValueError(f"{path}: duplicate unordered LD pair")
    return pd.DataFrame(
        {
            "idx1": np.concatenate(idx1_parts),
            "idx2": np.concatenate(idx2_parts),
            "r": np.concatenate(r_parts),
        }
    )


def _vcor_r_column(header, path: Path) -> str:
    candidates = ["UNPHASED_R", "PHASED_R", "R"]
    for name in candidates:
        if name in header:
            return name
    raise ValueError(f"{path}: missing signed r column in .vcor")


def _read_table_with_hash_header(path: Path) -> dict:
    header, chunks = _iter_table_chunks_with_hash_header(path, chunk_rows=None)
    rows = []
    for chunk in chunks:
        rows.extend(chunk.to_dict("records"))
    return {"header": header, "rows": rows}


def _iter_table_chunks_with_hash_header(path: Path, chunk_rows):
    if not path.is_file():
        raise FileNotFoundError(f"expected PLINK output file not found: {path}")
    with open(path, newline="") as f:
        first = f.readline()
        if first == "":
            raise ValueError(f"{path}: empty table")
        header = first.rstrip("\n\r").split()
        if header and header[0].startswith("#"):
            header[0] = header[0][1:]
    extra_col = "__statgen_extra_column__"
    while extra_col in header:
        extra_col = f"_{extra_col}"
    read_names = [*header, extra_col]

    def chunks():
        kwargs = {
            "sep": r"\s+",
            "header": None,
            "names": read_names,
            "skiprows": 1,
            "dtype": str,
            "keep_default_na": False,
            "na_filter": False,
            "index_col": False,
            "engine": "c",
        }
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("error", category=pd.errors.ParserWarning)
                if chunk_rows is None:
                    yield _validate_table_chunk(path, pd.read_csv(path, **kwargs), header, extra_col)
                    return
                reader = pd.read_csv(path, chunksize=chunk_rows, **kwargs)
                for chunk in reader:
                    yield _validate_table_chunk(path, chunk, header, extra_col)
        except pd.errors.EmptyDataError:
            return
        except pd.errors.ParserWarning as exc:
            raise ValueError(f"{path}: malformed table row") from exc

    return set(header), chunks()


def _validate_table_chunk(path: Path, chunk: pd.DataFrame, header: list[str], extra_col: str) -> pd.DataFrame:
    if (chunk[extra_col] != "").to_numpy().any() or chunk[header].eq("").to_numpy().any():
        raise ValueError(f"{path}: malformed table row")
    return chunk.loc[:, header]


def _reference_key_index(reference_shard) -> tuple[np.ndarray, np.ndarray]:
    keys = np.empty(reference_shard.num_snp, dtype=VARIANT_KEY_DTYPE)
    keys["chr"] = reference_shard.chr
    keys["bp"] = reference_shard.bp
    keys["a1"] = reference_shard.a1
    keys["a2"] = reference_shard.a2
    order = np.argsort(keys, order=VARIANT_KEY_ORDER)
    sorted_keys = keys[order]
    if sorted_keys.size > 1 and np.any(sorted_keys[:-1] == sorted_keys[1:]):
        raise ValueError("reference shard contains duplicate chr:bp:a1:a2 keys")
    return sorted_keys, order.astype(np.int32, copy=False)


def _variant_key_array(df: pd.DataFrame, chr_col: str, bp_col: str, a1_col: str, a2_col: str) -> np.ndarray:
    keys = np.empty(len(df), dtype=VARIANT_KEY_DTYPE)
    keys["chr"] = df[chr_col].to_numpy(dtype=object, copy=False)
    keys["bp"] = pd.to_numeric(df[bp_col], errors="raise").to_numpy(dtype=np.int64)
    keys["a1"] = df[a1_col].to_numpy(dtype=object, copy=False)
    keys["a2"] = df[a2_col].to_numpy(dtype=object, copy=False)
    return keys


def _map_variant_keys_to_reference(
    keys: np.ndarray,
    sorted_ref_keys: np.ndarray,
    sorted_ref_idx: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    pos = np.searchsorted(sorted_ref_keys, keys)
    valid = pos < sorted_ref_keys.size
    matched = valid.copy()
    matched[valid] = sorted_ref_keys[pos[valid]] == keys[valid]
    idx = np.empty(keys.size, dtype=np.int32)
    idx[matched] = sorted_ref_idx[pos[matched]]
    return idx, ~matched


def _format_command(cmd) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


if __name__ == "__main__":
    raise SystemExit(main())
