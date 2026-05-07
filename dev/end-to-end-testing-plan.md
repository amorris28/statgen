# End-to-end integration testing plan

This plan describes high-value end-to-end tests for global alignment
invariants across reference, genotype, LD, Python, and MATLAB/Octave runtimes.
The tests are intentionally synthetic and generated on the fly because their
purpose is to stress re-indexing, sharding, and subject mapping logic rather
than to define portable fixture bytes.

The main invariant is that `a1freq` computed by PLINK2 during LD construction
and embedded in an `LDPanel` must match the `a1` dosage frequency recomputed
from the same genotype payload through `GenotypePanel.fetch_genotypes`. This
turns allele frequency into a compact, variant-specific checksum of the
genotype vector and sample subset used for that LD shard.

## Fixture shape

Generate one deterministic master VCF/PSAM pair and derive all bfiles from it
with PLINK2.

Recommended master data:

- `200` total variants across shards `1`, `2`, and `X`;
- approximately `65` variants on chr1, `63` on chr2, and `72` on chrX;
- `37` to `43` subjects, for example `20` male, `19` female, and `2` unknown
  sex subjects;
- unique marker names for every variant, using a stable pattern such as
  `e2e_chr1_0001_A_G`, so failures are easy to trace in test code and PLINK
  output;
- deterministic seeded genotypes with post-processing to avoid monomorphic
  variants in the ordinary build path;
- missing calls in selected variants;
- chrX genotypes where male and female allele frequencies differ visibly;
- sentinel correlated SNP pairs for limited LD `r` checks.

Derive two partially overlapping analysis panels from the master variant set:

- a genotype panel subset with `150` variants;
- an LD panel subset with `150` variants;
- exactly `100` variants shared by both derived panels;
- genotype-only and LD/reference-only variants distributed across chr1, chr2,
  and chrX.

With a `200`-variant master universe, both derived panels should be `150`
variants. A larger LD panel, such as `1500` variants, is useful only if the
master universe is expanded accordingly. The important invariant is the partial
overlap, not the exact absolute size.

Use PLINK2 `--extract` to create the derived bfiles. This stresses projection
from a larger source and ensures genotype and LD panels have independent source
row indices while retaining a controlled overlap.

## PLINK conversion

Prefer generating VCF/PSAM and using PLINK2 to create `.bed/.bim/.fam` bfiles
instead of writing PLINK BED bytes directly. This reduces fixture-writer risk
and exercises the same representation expected from real users.

The generator must verify after `--make-bed` that PLINK2 did not reorder or
reinterpret alleles relative to the intended contract:

- BIM `a1` must match VCF `ALT`;
- BIM `a2` must match VCF `REF`;
- marker names and positions must match the generated manifest;
- `--keep-allele-order` must be used in PLINK2 frequency and LD commands;
- if any PLINK2 version changes this behavior, the test should fail clearly
  before running statgen assertions.

Write `.ploidy` sidecars directly as deterministic text. For chrX, use male
ploidy `1` and female ploidy `2` for the ordinary test variants.

## Same-position allele clusters

Each shard should contain at least one cluster with many variants at the same
`chr:bp` and different `(a1, a2)` pairs. Use about `10` to `12` variants in
each cluster.

The cluster should be designed so that:

- allele pairs require matching by `(bp, a1_hash64, a2_hash64)`, not by
  marker name or physical row position;
- hash-sort order differs from original VCF/BIM row order for at least one
  cluster;
- every clustered variant has a distinctive genotype vector and allele
  frequency;
- sorting or matching keys cannot disconnect genotype payload row indices from
  their variants without changing the recomputed `a1freq`.

The generator should assert the intended hash-order perturbation before the
test proceeds. This makes the test explicit rather than relying on accidental
allele string ordering.

## Core test cases

### 0. Pre-flight one-shard baseline

Start with a deliberately simple chr1-only test before the full re-indexing
workflow. This test should prove that PLINK2 conversion, LD building, LD
loading, genotype loading, frequency recomputation, and selected LD `r`
comparisons agree in the simplest possible case.

Recommended shape:

- one shard: chr1 only;
- one non-sharded bfile;
- reference, genotype, and LD variants are identical;
- exactly `31` SNPs;
- exactly `29` subjects, for example `14` male, `14` female, and `1` unknown
  sex subject;
- several missing calls;
- two or three sentinel LD pairs covering positive and negative correlation,
  and, if easy to control, one below-threshold omitted pair.

Assertions:

- PLINK2 `--make-bed` output preserves the allele contract: BIM `a1 == ALT`
  and BIM `a2 == REF`;
- all variants have `GenotypePanel.is_present == true`;
- `LDPanel.a1freq` matches frequency recomputed from
  `GenotypePanel.fetch_genotypes` for every SNP;
- selected `ld_r` entries match direct Pearson `r` from fetched genotypes;
- `LDPanel.multiply_r2(M)` matches an explicit sparse calculation for small
  deterministic vector and matrix inputs.

This case should contain no chrX, no sharded input, no variant subset mismatch,
and no same-position allele clusters. If this test fails, the likely problem is
basic PLINK2 conversion/build/load or numerical comparison logic. The
full-case tests below should run only after the pre-flight path is reliable.

### 1. Non-sharded genotype source, LD built from extracted source

Create one non-sharded genotype bfile from the genotype extract and one
non-sharded LD-source bfile from the LD extract. Build LD independently for
shards `1`, `2`, and `X`.

Assertions:

- `load_ld` loads the LD-derived reference in shard order `1`, `2`, `X`;
- `load_genotype(genotype_bfile, ld.reference)` projects genotype rows onto the
  LD reference;
- for the `100` shared variants, `LDPanel.a1freq` matches frequency recomputed
  from `GenotypePanel.fetch_genotypes`;
- LD/reference-only variants have `GenotypePanel.is_present == false` and are
  excluded from frequency recomputation;
- same-position allele clusters retain their expected marker-to-frequency
  mapping.

This case primarily stresses non-sharded BIM splitting and source-to-reference
projection.

### 2. Sharded genotype source with chrX FAM subset

Create sharded genotype bfiles with `@` substitution. Autosomal shards use the
same full FAM in the same order. The chrX shard uses a non-contiguous subset of
subjects in a deliberately different FAM order.

Assertions:

- `GenotypePanel.is_subject_present("X")` matches the intended panel-axis
  subject mask;
- chrX `fetch_genotypes` returns rows in panel sample order and fills subjects
  absent from chrX with missing calls;
- LD chrX `female` and `male` `a1freq` values match recomputation using both
  sex masks and the chrX subject-present mask;
- autosomal frequencies match the non-sharded source case for shared variants.

This is the highest-value test for subject-index bugs.

### 3. Mixed panel-global SNP requests

Use the same loaded genotype and LD objects. Fetch a deliberately awkward list
of present SNP indices:

- last present SNP on chr1;
- first present SNP on chrX;
- middle present SNP on chr2;
- repeated chr1 SNP;
- repeated chrX SNP;
- several same-position cluster members.

Assertions:

- returned genotype columns follow the requested order exactly;
- repeated SNP requests return repeated columns;
- frequencies recomputed from the selected columns map back to the same LD
  `a1freq` entries.

This catches panel-global to shard-local splitting and reassembly errors.

### 4. Python versus MATLAB/Octave consistency

Read the same generated source bytes in Python and Octave. LD is runtime
native, so Python reads `.npz` and Octave reads converted `.mat` files derived
from the same `.npz` output.

Compare:

- reference shard labels, offsets, positions, alleles, and checksums;
- genotype `is_present`, ploidy vectors, and `subject_present` masks;
- selected `fetch_genotypes_int8` slices;
- LD `a1freq` for autosomes and chrX `female`/`male`;
- `multiply_r2` on deterministic vector and matrix inputs.

The Octave side should print compact arrays or save a small temporary output
file for pytest to compare. The test should skip cleanly when Octave is not
available.

## Frequency recomputation

Recompute frequencies from `fetch_genotypes`, not from internal BED helpers.
This keeps the test on public APIs.

For autosomes:

```text
a1freq = sum(a1_count over finite calls) / (2 * number of finite calls)
```

For chrX female shards:

```text
a1freq = sum(a1_count over present finite female calls)
         / (2 * number of present finite female calls)
```

For chrX male shards:

```text
a1freq = sum(a1_count over present finite male calls)
         / (1 * number of present finite male calls)
```

For optional chrX combined shards:

```text
a1freq = sum(a1_count over present finite male/female calls)
         / sum(ploidy over the same observed calls)
```

Unknown-sex subjects should participate in autosomal frequency checks but not
in sex-specific chrX frequency checks.

## LD r checks

Do not recompute the full LD matrix. Add a small number of sentinel checks:

- one positive autosomal correlation;
- one negative autosomal correlation;
- one below-threshold pair expected to be absent from the sparse matrix and
  read as zero;
- one chrX pair where `female` and `male` `r` differ.

Recompute Pearson `r` directly from fetched genotype vectors using the same
sample subset as the corresponding LD shard. Compare only those selected sparse
entries. Also compare `LDPanel.multiply_r2(M)` with an explicit per-shard
`elementwise(ld_r .^ 2) * M` computation for a small deterministic vector and
matrix.

## Tolerances

Use tolerance based on PLINK2 text output precision and runtime dtype:

- `a1freq`: absolute tolerance `1e-5`;
- selected LD `r`: absolute tolerance `1e-3`;
- `multiply_r2`: absolute tolerance `1e-5` for `float64` comparisons and
  `1e-4` for `float32` comparisons;
- Python versus Octave LD matrix and `a1freq` comparisons should use the same
  numeric tolerances, not dtype identity.

PLINK2 versions may print `ALT1_FREQ` with more precision than `1e-5`; the
fixed tolerance is intentionally conservative for frequency checks while still
small enough to detect one-subject or one-row indexing errors at the proposed
sample sizes. PLINK2 LD text commonly prints fewer decimals for `r`, so `1e-3`
is the safer comparison tolerance for selected `r` checks.

## Cost control

Keep this as one generated fixture family with multiple assertions rather than
several independent large fixtures. Tests requiring real PLINK2 should be
optional and skipped when PLINK2 is unavailable. Negative validation cases are
not part of this plan because unit tests already cover malformed inputs and
failure paths.
