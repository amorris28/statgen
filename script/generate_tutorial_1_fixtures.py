#!/usr/bin/env python3
"""Generate minimal Tutorial 1 source fixtures.

Writes all source files that Tutorial 1 reads as inputs:

    source/ld_reference/chr{21,22,X}.{bed,bim,fam}
    source/genotypes/chr{21,22,X}.{bed,bim,fam}
    source/genotypes/chrX.ploidy
    source/sumstats/trait_a.tsv.gz
    source/sumstats/trait_b.tsv.gz
    source/annotations/{coding_exon,exon,intron,utr3,utr5,whole_gene}.bed
    source/annotations/functional_groups.annot
    source/annotations/conservation.annot
    source/annotations/conservation.meta

Design highlights
-----------------
Different variant counts per chromosome: chr21=40, chr22=26, chrX=20.

Partial overlap between reference, genotype, and sumstats:
  Each bfile has N SNPs, the reference has N SNPs, but only N/2 are shared.
  The first N/2 reference variants are matched in the genotype bfile
  (is_present=True); the last N/2 are absent (is_present=False).
  The last N/2 reference variants are matched in the sumstats file; the
  first N/2 are absent (logpvec=NaN). Both files also carry N/2 variants
  that are not in the reference at all (no match → silently ignored).
  This exercises the reference-alignment feature of load_genotype and
  load_sumstats.

chrX genotype bfile uses only 50 of the 60 panel subjects.  Every 6th
  subject is absent from genotypes/chrX.fam.  load_genotype detects the
  missing 10 via FAM matching and sets subject_present=False for those rows.

LD reference and annotation are always reference-complete (full N SNPs,
  all 6 annotation BEDs based on reference positions).

Each annotation has a different number of genomic block intervals and a
  different fraction of reference variants covered.
The functional_groups annotation is a grouped BED-like file where one selected
  column expands to multiple binary annotations. The conservation annotation is
  a headered BED-like file with two continuous numeric value columns and a
  column-aligned metadata sidecar.

PLINK bfiles are generated directly in Python; no external tools are needed
for this step. plink2 is required later for the Tutorial 1 LD build step.

Usage (run from the repository root):
    python script/generate_tutorial_1_fixtures.py [--out OUTDIR]

Set the Tutorial 1 working directory to OUTDIR after running.
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SHARDS = ["21", "22", "X"]
N_VARIANTS = {"21": 40, "22": 26, "X": 20}  # different per shard, all even
N_SAMPLES = 60          # 30 male (sex=1) + 30 female (sex=2)
N_MALE = N_SAMPLES // 2
N_SAMPLES_CHRX_GENO = 50  # chrX genotype bfile uses a subset of individuals
RNG_SEED = 42

# Realistic-ish GRCh38 starting positions per chromosome
CHR_START_BP = {"21": 5_010_000, "22": 10_510_000, "X": 2_700_000}
BP_STEP = 100_000  # 100 kbp between variants

_ALLELE_PAIRS = [("A", "G"), ("T", "C"), ("A", "C"), ("T", "G"), ("G", "C"), ("A", "T")]

# Indices of the 50 subjects included in the chrX genotype bfile.
# Every 6th subject (1-based positions 6, 12, ..., 60) is absent, leaving 50.
_CHRX_GENO_SAMPLE_IDX = np.array([i for i in range(N_SAMPLES) if (i + 1) % 6 != 0])

# Each annotation is defined by (name, n_blocks, block_size_in_variants, start_variant_offset).
# Intervals are 0-based half-open BED coordinates spanning block_size consecutive variants.
# Different n_blocks and block_size give visibly different BED line counts and coverage fractions.
_ANNOTATION_SPEC = [
    #  name           n_blocks  block_size  start_offset
    ("coding_exon",   5,        1,          2),  # 5 single-SNP intervals per chr
    ("exon",          8,        1,          0),  # 8 single-SNP intervals per chr
    ("intron",        3,        4,          5),  # 3 × 4-SNP blocks per chr
    ("utr3",          4,        1,          1),  # 4 single-SNP intervals per chr
    ("utr5",          3,        1,          4),  # 3 single-SNP intervals per chr
    ("whole_gene",    2,        6,          3),  # 2 × 6-SNP blocks per chr
]

ANNOTATION_NAMES = [spec[0] for spec in _ANNOTATION_SPEC]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--out",
        default="docs/tutorial_1_fixtures",
        help="Output directory (default: %(default)s)",
    )
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Variant helpers
# ---------------------------------------------------------------------------

def _allele_pair(bp: int) -> tuple[str, str]:
    return _ALLELE_PAIRS[(bp // BP_STEP) % len(_ALLELE_PAIRS)]


def _make_variant(chr_label: str, bp: int) -> dict:
    a1, a2 = _allele_pair(bp)
    return {"chr": chr_label, "bp": bp, "snp": f"{chr_label}:{bp}:{a1}:{a2}", "a1": a1, "a2": a2}


def _make_variants_at(chr_label: str, start_bp: int, count: int) -> list[dict]:
    return [_make_variant(chr_label, start_bp + i * BP_STEP) for i in range(count)]


# ---------------------------------------------------------------------------
# PLINK bfile writers
# ---------------------------------------------------------------------------

def _write_bim(path: Path, variants: list[dict]) -> None:
    with open(path, "w") as f:
        for v in variants:
            f.write(f"{v['chr']}\t{v['snp']}\t0\t{v['bp']}\t{v['a1']}\t{v['a2']}\n")


def _write_fam(path: Path, sample_indices=None) -> None:
    """Write FAM for the given sample indices (default: all N_SAMPLES subjects)."""
    indices = range(N_SAMPLES) if sample_indices is None else sample_indices
    with open(path, "w") as f:
        for i in indices:
            iid = f"SAMPLE{i + 1:04d}"
            sex = 1 if i < N_MALE else 2
            f.write(f"{iid}\t{iid}\t0\t0\t{sex}\t-9\n")


def _pack_bed_row(genos: np.ndarray) -> bytes:
    """Pack one SNP-major BED row.

    genos: int8 array of length n_samples with values 0, 1, 2 (copies of a1)
    or -1 (missing).

    BED 2-bit encoding per sample:
        00 → 2 copies of a1 (hom effect allele)
        01 → missing
        10 → 1 copy of a1 (het)
        11 → 0 copies of a1 (hom reference allele)
    Four samples packed per byte, LSB first.
    """
    n = len(genos)
    two_bit = np.full(n, 3, dtype=np.uint8)  # 11 = hom a2 (0 copies of a1)
    two_bit[genos == 2] = 0   # 00 = hom a1
    two_bit[genos == -1] = 1  # 01 = missing
    two_bit[genos == 1] = 2   # 10 = het
    padded_len = ((n + 3) // 4) * 4
    padded = np.zeros(padded_len, dtype=np.uint8)
    padded[:n] = two_bit
    padded = padded.reshape(-1, 4)
    byte_vals = (
        padded[:, 0]
        | (padded[:, 1] << 2)
        | (padded[:, 2] << 4)
        | (padded[:, 3] << 6)
    ).astype(np.uint8)
    return byte_vals.tobytes()


def _write_bed(path: Path, geno_matrix: np.ndarray) -> None:
    """Write PLINK BED file in SNP-major mode.

    geno_matrix: int8 array of shape (n_samples, n_snp).
    """
    with open(path, "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))  # magic bytes + SNP-major flag
        n_snp = geno_matrix.shape[1]
        for snp_idx in range(n_snp):
            f.write(_pack_bed_row(geno_matrix[:, snp_idx]))


# ---------------------------------------------------------------------------
# chrX ploidy sidecar
# ---------------------------------------------------------------------------

def _write_ploidy(path: Path, n_snp: int) -> None:
    """Write chrX ploidy sidecar: male_ploidy=1, female_ploidy=2 for all variants."""
    with open(path, "w") as f:
        for _ in range(n_snp):
            f.write("1\t2\n")


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def _write_sumstats(path: Path, variants: list[dict], rng: np.random.Generator) -> None:
    """Write gzip-compressed sumstats TSV with chr, bp, a1, a2, p, z, n."""
    n = len(variants)
    pvec = rng.uniform(1e-8, 0.99, size=n)
    zvec = rng.standard_normal(size=n)
    with gzip.open(path, "wt") as f:
        f.write("chr\tbp\ta1\ta2\tp\tz\tn\n")
        for v, p, z in zip(variants, pvec, zvec):
            f.write(
                f"{v['chr']}\t{v['bp']}\t{v['a1']}\t{v['a2']}"
                f"\t{p:.6g}\t{z:.6f}\t50000\n"
            )


# ---------------------------------------------------------------------------
# Annotation BED files (always reference-complete)
# ---------------------------------------------------------------------------

def _annotation_intervals_for_chr(
    variants: list[dict],
    n_blocks: int,
    block_size: int,
    start_offset: int,
) -> list[tuple[int, int]]:
    """Return 0-based half-open [start, end) intervals for one chromosome.

    Each interval spans exactly block_size consecutive reference variants
    beginning at start_offset, then at start_offset + step, etc.
    A variant at bp X is covered when start <= (X - 1) < end.
    """
    n = len(variants)
    remaining = n - start_offset
    if remaining <= 0:
        return []
    step = max(block_size, remaining // n_blocks)
    intervals = []
    for i in range(n_blocks):
        first_idx = start_offset + i * step
        if first_idx >= n:
            break
        last_idx = min(first_idx + block_size, n)
        bp_start = variants[first_idx]["bp"] - 1     # 0-based, includes bp[first]
        bp_end = variants[last_idx - 1]["bp"]        # exclusive; pos0 = bp-1 < bp = end
        intervals.append((bp_start, bp_end))
    return intervals


def _write_annotation_beds(annot_dir: Path, ref_variants_by_chr: dict) -> None:
    annot_dir.mkdir(parents=True, exist_ok=True)
    total_ref = sum(len(v) for v in ref_variants_by_chr.values())

    for name, n_blocks, block_size, start_offset in _ANNOTATION_SPEC:
        total_intervals = 0
        total_covered = 0
        with open(annot_dir / f"{name}.bed", "w") as f:
            for chr_label, variants in ref_variants_by_chr.items():
                ivs = _annotation_intervals_for_chr(
                    variants, n_blocks, block_size, start_offset
                )
                for bp_start, bp_end in ivs:
                    f.write(f"{chr_label}\t{bp_start}\t{bp_end}\n")
                    total_intervals += 1
                    total_covered += sum(
                        1 for v in variants if bp_start <= v["bp"] - 1 < bp_end
                    )
        pct = 100.0 * total_covered / total_ref if total_ref else 0
        print(
            f"  {name:<14} {total_intervals:3d} intervals,"
            f" {total_covered}/{total_ref} variants ({pct:.0f}%)"
        )


def _write_continuous_annotation(annot_dir: Path) -> None:
    """Write a small headered BED-like continuous annotation plus sidecar."""
    rows = [
        ("21", 5_000_000, 5_400_000, 0.82, 1.5),
        ("21", 5_600_000, 6_100_000, 0.35, 0.4),
        ("22", 10_500_000, 11_000_000, 0.67, 2.2),
        ("22", 12_100_000, 12_650_000, 0.18, 0.7),
        ("X", 2_600_000, 3_300_000, 0.91, 1.9),
        ("X", 3_900_000, 4_300_000, 0.44, 0.2),
    ]
    with open(annot_dir / "conservation.annot", "w") as f:
        f.write("chrom\tstart0\tend0\tconservation\tpromoter_activity\n")
        for chrom, start0, end0, conservation, promoter_activity in rows:
            f.write(
                f"{chrom}\t{start0}\t{end0}\t"
                f"{conservation:g}\t{promoter_activity:g}\n"
            )

    with open(annot_dir / "conservation.meta", "w") as f:
        f.write('{"column":"chrom","role":"coordinate"}\n')
        f.write('{"column":"start0","role":"coordinate"}\n')
        f.write('{"column":"end0","role":"coordinate"}\n')
        f.write(
            '{"name":"conservation","type":"continuous",'
            '"description":"Synthetic conservation-like interval score"}\n'
        )
        f.write(
            '{"name":"promoter_activity","type":"continuous",'
            '"description":"Synthetic promoter-activity interval score"}\n'
        )


def _write_grouped_annotation(annot_dir: Path) -> None:
    rows = [
        ("21", 5_000_000, 5_500_000, "coding"),
        ("21", 5_700_000, 6_300_000, "regulatory"),
        ("22", 10_500_000, 11_100_000, "coding"),
        ("X", 2_600_000, 3_400_000, "regulatory"),
    ]
    with open(annot_dir / "functional_groups.annot", "w") as f:
        f.write("chrom\tstart0\tend0\tgroup\n")
        for chrom, start0, end0, group in rows:
            f.write(f"{chrom}\t{start0}\t{end0}\t{group}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = parse_args(argv)
    out_root = Path(args.out)
    rng = np.random.default_rng(RNG_SEED)

    src = out_root / "source"
    ld_ref_dir = src / "ld_reference"
    geno_dir = src / "genotypes"
    sumstats_dir = src / "sumstats"
    annot_dir = src / "annotations"

    for d in [ld_ref_dir, geno_dir, sumstats_dir, annot_dir]:
        d.mkdir(parents=True, exist_ok=True)

    ref_variants_by_chr: dict[str, list[dict]] = {}
    ss_variants_all: list[dict] = []

    print("Generating bfiles and sumstats variants:")
    for chr_label in SHARDS:
        n_ref = N_VARIANTS[chr_label]
        n_half = n_ref // 2
        start_bp = CHR_START_BP[chr_label]

        # ── Reference variants (defines the panel) ────────────────────────────
        ref_variants = _make_variants_at(chr_label, start_bp, n_ref)
        ref_variants_by_chr[chr_label] = ref_variants

        # ── Genotype bfile: first n_half of ref + n_half novel after ref range ─
        # Variants ref[0:n_half] → is_present=True in genotype shard.
        # Variants ref[n_half:]  → is_present=False (no match in genotype BIM).
        # Novel variants are in the genotype BIM but outside the reference → ignored.
        geno_variants = (
            ref_variants[:n_half]
            + _make_variants_at(chr_label, start_bp + n_ref * BP_STEP, n_half)
        )  # sorted by bp ✓

        # ── Sumstats: last n_half of ref + n_half novel before ref range ───────
        # Variants ref[n_half:] → logpvec=finite in sumstats shard.
        # Variants ref[0:n_half] → logpvec=NaN (absent from sumstats file).
        # Novel variants in the TSV are outside the reference → ignored.
        ss_novel_start_bp = start_bp - n_half * BP_STEP
        ss_variants_all.extend(
            ref_variants[n_half:]
            + _make_variants_at(chr_label, ss_novel_start_bp, n_half)
        )

        # ── Genotype matrices ──────────────────────────────────────────────────
        def _geno_matrix(n_snp, n_samp):
            freq = rng.uniform(0.1, 0.9, size=n_snp)
            return np.column_stack(
                [rng.binomial(2, f, size=n_samp).astype(np.int8) for f in freq]
            )

        # LD reference: all 60 samples, all n_ref reference variants
        ldref_mat = _geno_matrix(n_ref, N_SAMPLES)
        _write_bim(ld_ref_dir / f"chr{chr_label}.bim", ref_variants)
        _write_fam(ld_ref_dir / f"chr{chr_label}.fam")
        _write_bed(ld_ref_dir / f"chr{chr_label}.bed", ldref_mat)

        # Genotypes: partially overlapping BIM, subset of samples for chrX
        geno_mat_full = _geno_matrix(n_ref, N_SAMPLES)  # n_ref = n_half + n_half
        if chr_label == "X":
            geno_mat = geno_mat_full[_CHRX_GENO_SAMPLE_IDX, :]
            n_geno_samples = N_SAMPLES_CHRX_GENO
            _write_fam(geno_dir / f"chr{chr_label}.fam", _CHRX_GENO_SAMPLE_IDX)
            _write_ploidy(geno_dir / "chrX.ploidy", n_ref)
        else:
            geno_mat = geno_mat_full
            n_geno_samples = N_SAMPLES
            _write_fam(geno_dir / f"chr{chr_label}.fam")
        _write_bim(geno_dir / f"chr{chr_label}.bim", geno_variants)
        _write_bed(geno_dir / f"chr{chr_label}.bed", geno_mat)

        print(
            f"  chr{chr_label}: ref={n_ref} SNPs ({n_half} matched by geno, {n_half} matched by sumstats),"
            f" geno BIM={n_ref} SNPs ({n_half} novel), ld_ref={N_SAMPLES} samples,"
            f" geno={n_geno_samples} samples"
        )

    print()
    _write_sumstats(sumstats_dir / "trait_a.tsv.gz", ss_variants_all, rng)
    _write_sumstats(sumstats_dir / "trait_b.tsv.gz", ss_variants_all, rng)
    n_total_ref = sum(N_VARIANTS[c] for c in SHARDS)
    n_total_ss_shared = n_total_ref // 2
    print(
        f"  sumstats: {len(ss_variants_all)} rows"
        f" ({n_total_ss_shared} match reference, {len(ss_variants_all) - n_total_ss_shared} novel)"
    )

    print()
    _write_annotation_beds(annot_dir, ref_variants_by_chr)
    _write_grouped_annotation(annot_dir)
    _write_continuous_annotation(annot_dir)

    print(f"\nFixtures written to: {src}")
    print(f"Use '{out_root}' as the Tutorial 1 working directory.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
