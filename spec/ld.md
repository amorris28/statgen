# LDShard and LDPanel

## Represents

LD describes the pairwise signed correlation structure between variants in the
same processed shard. Values are signed correlation `r`, not `r²`; consumers
that need `r²` compute it after loading. Cross-chromosome LD is not stored.
Signs are oriented to `a1` dosages, consistent with the allele contract in
[SPEC.md](SPEC.md).

LD is paired with a reference: matrix row/column order, SNP count, allele
frequencies, and `reference_checksum` correspond to a specific
`ReferenceShard`.

LD panels are a deliberate exception to the package-wide portable-source plus
optional-cache pattern. LD distribution artifacts are runtime-native sparse
matrix files because the expensive sparse-matrix construction is part of
reference-panel preparation, not routine user loading.

## Logical LD shard contract

Each LD shard represents:

- `chr`: shard label, one of `1`-`22` or `X`;
- `sex`: `null` for autosomes, one of `"female"`, `"male"`, or `"combined"`
  for chrX;
- `num_snp`: number of SNPs in the paired reference shard;
- `ld_r`: symmetric sparse signed-correlation matrix, shape
  `num_snp x num_snp`;
- diagonal entries of `ld_r`: explicit `1.0`;
- off-diagonal entries: signed `r` oriented to `a1` dosage;
- missing off-diagonal entries: interpreted as zero;
- `a1freq`: float vector of `a1` allele frequencies aligned to the paired
  `ReferenceShard`;
- `reference_checksum`: MD5 checksum of the paired reference shard;
- build metadata: PLINK command, PLINK version if available, window size,
  `r²` storage threshold, sample count, and schema version.

The paired `ReferencePanel` and `LDPanel` must originate from the same BIM row
order and allele orientation. The LD build additionally depends on the paired
BED/FAM genotype source used to estimate LD and allele frequencies.

## Distribution formats

`load_ld` reads the distribution artifact native to the runtime:

- Python reads NumPy/SciPy `.npz` files with CSC sparse components.
- MATLAB/Octave reads `.mat` files containing native sparse matrices.

The Python `.npz` format is also the documented build handoff format from
`statgen_build_ld.py` into the MATLAB/Octave converter. MATLAB/Octave user
loading reads `.mat`, not `.npz`.

### Python `.npz` shard format

Each `.npz` file represents one LD shard. Required arrays:

```text
data       float32, length nnz
indices    int32, length nnz, zero-based row indices
indptr     int32, length num_snp + 1, zero-based CSC column pointers
shape      int64, length 2, equals [num_snp, num_snp]
a1freq     float32, length num_snp
metadata   uint8, UTF-8 JSON bytes
```

`data`, `indices`, and `indptr` encode a symmetric sparse signed-`r` matrix in
CSC layout with explicit unit diagonal. Sparse indices are zero-based because
the `.npz` format is Python/SciPy-oriented. The archive must not contain object
arrays or pickled payloads.

CSC32 shards must have `nnz < 2^31`. Builders must reject larger shards rather
than overflowing sparse indices; a future CSC64 schema would be a separate
format/version.

The metadata JSON must include:

```json
{
  "object_type": "ld_shard",
  "schema_version": "1.0",
  "format": "statgen_ld_npz_csc32",
  "chr": "1",
  "sex": null,
  "num_snp": 123,
  "nnz": 4567,
  "sparse_layout": "csc",
  "index_base": 0,
  "matrix": "symmetric",
  "diagonal": "explicit_unit",
  "value": "r",
  "reference_checksum": "...",
  "build_tool": "plink2",
  "build_command": "plink2 ...",
  "plink_version": "...",
  "ld_window_kb": 10000,
  "ld_r2_threshold": 0.05,
  "num_sample": 1000
}
```

`plink_version` may be `null` if unavailable.

### MATLAB/Octave `.mat` shard format

Each `.mat` file represents one LD shard. Required variables:

```text
ld_r       sparse double, num_snp x num_snp
a1freq     numeric vector, length num_snp
metadata   struct
```

The `.mat` file must be MAT-file version 7.3 (HDF5-based). The legacy v5
MAT-file format is not acceptable because realistic LD shards can exceed its
variable-size limits.

`ld_r` must be native sparse, symmetric, signed `r`, explicitly unit-diagonal,
and aligned to the paired reference shard. The loader must load `ld_r`
directly; it must not reconstruct the sparse matrix from COO or CSC arrays on
the default user path.

`metadata` must contain the same logical fields as the `.npz` metadata, with
MATLAB-valid struct field names and `format: "statgen_ld_mat_sparse_double"`.

## Panel layout and manifest

A panel root is a directory containing one or more LD shard files plus
`ld_manifest.json`. A loader may also accept a single shard file.

Recommended shard filenames:

```text
ld_chr1.npz
ld_chr2.npz
...
ld_chr22.npz
ld_chrX_female.npz
ld_chrX_male.npz
ld_chrX_combined.npz
```

MATLAB/Octave distributions use the same names with `.mat`. Autosomal shards
use no sex suffix. chrX shards always use an explicit sex suffix. Not all chrX
sex labels need to be present. `female` and `male` are the default
sex-specific build outputs; `combined` is optional and assumption-dependent.

For panel-root loads, `ld_manifest.json` is authoritative for file discovery.
Loaders must not infer panels by globbing arbitrary filenames. Files present in
the directory but absent from the manifest are ignored.

Required manifest fields:

```json
{
  "object_type": "ld_panel_manifest",
  "schema_version": "1.0",
  "runtime_format": "python_npz_csc32",
  "shards": [
    {
      "chr": "1",
      "sex": null,
      "file": "ld_chr1.npz",
      "file_md5": "...",
      "num_snp": 123,
      "nnz": 4567,
      "reference_checksum": "..."
    }
  ]
}
```

`runtime_format` is `"python_npz_csc32"` for Python distributions and
`"matlab_mat_sparse_double"` for MATLAB/Octave distributions. Manifest entries
and per-file metadata must agree on `chr`, `sex`, `num_snp`, `nnz`,
`reference_checksum`, and runtime/storage format. `file` and `file_md5` are
manifest-only fields. Runtime loaders are not required to compute `file_md5` on
the default load path.

If the supplied `ReferencePanel` is a shard subset, `load_ld` loads only the LD
files required for those reference shards. LD files for other reference shards
may be present in the directory and are ignored.

## Building and conversion

`statgen_build_ld.py` builds the Python `.npz` LD distribution and handoff
files from a PLINK bfile. It accepts both input layouts:

- **Sharded bfile input** (`@` in path, e.g. `chr@`): one bfile per
  chromosome; discovery is by canonical substitution (`1`-`22`, `X`) and each
  discovered shard is processed independently.
- **Non-sharded bfile input** (single path, no `@`): the script splits by
  chromosome.

For each processed unit, the builder:

1. Computes `a1` allele frequencies using PLINK2 `--freq`, aligned to BIM row
   order.
2. Computes pairwise signed LD using PLINK2 `--r` and `--keep-allele-order`
   with default window 10,000 kb and default `r²` storage threshold `0.05`.
3. Writes one validated `.npz` shard file and updates `ld_manifest.json`.

Before publishing an `.npz` shard, the builder must validate at least:

- all required array names are present;
- `shape == [num_snp, num_snp]`;
- `len(indptr) == num_snp + 1`;
- `len(data) == len(indices) == indptr[-1]`;
- `len(a1freq) == num_snp`;
- sparse index arrays use int32 and `nnz < 2^31`;
- sparse index bounds are valid;
- metadata schema version is supported;
- metadata and array dimensions agree.

The r² threshold, window size, PLINK2 command, PLINK version if available,
sample count, and reference checksum are build metadata and must be recorded.
The sharding of the resulting LD panel must match the sharding of the reference
used with it.

By default, chrX builds produce `female` and `male` shards using FAM column 5
(1 = male, 2 = female). The script requires non-missing sex for chrX
sex-specific builds. `combined` chrX output is not the default; it is an
explicit user decision because it depends on modeling and encoding
assumptions. If users suppress sex splitting and write only a single combined
chrX shard, the manifest records `sex: "combined"` and the rationale must be
recorded in build metadata.

The MATLAB/Octave converter reads `.npz` shard files and writes `.mat` shard
files plus a MATLAB/Octave manifest. It must work under Octave for testability
and build reproducibility, but MATLAB remains the normative runtime where
MATLAB and Octave differ. The converter may assume `.npz` shards were produced
and validated by `statgen_build_ld.py`; it must validate metadata consistency
before writing but is not required to repeat expensive O(nnz) structure checks
such as full symmetry or diagonal scans.

## In-memory objects

An `LDShard` holds one sparse signed-`r` matrix and its aligned `a1freq`. Each
`LDShard` is unambiguously identified by its `(chr, sex)` pair.

An `LDPanel` is an ordered collection of shard lists matching the paired
reference panel. For autosomal chromosome shards, each chromosome has exactly
one `LDShard` (`sex = null`). For chrX, the panel may hold any non-empty subset
of `"female"`, `"male"`, and `"combined"`. Loading always loads the full set
of shards selected from disk for each requested reference shard.

`LDPanel` carries a `default_chrX_sex` field (`"female"`, `"male"`, or
`"combined"`) that selects which chrX shard downstream operations use when they
need a single LD matrix for chrX. If chrX is present in the panel,
`default_chrX_sex` must name a shard actually present in that chrX shard group.
If chrX is absent, `default_chrX_sex` is retained as metadata and ignored by
chrX-unrelated operations. For default sex-specific panels built by
`statgen_build_ld.py`, the default is `"female"`.

`LDShard` fields:

- `chr`: contig/shard label (`1`-`22` or `X`);
- `sex`: `None`/`null` for autosomal/sex-agnostic shards; `"female"`,
  `"male"`, or `"combined"` for chrX shards;
- `num_snp`: number of SNPs in the shard;
- `ld_r`: language-native sparse signed-correlation matrix;
- `a1freq`: vector of allele frequencies aligned to the paired
  `ReferenceShard`;
- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`.

Python stores `ld_r` as a SciPy CSC matrix from the `.npz` CSC components.
MATLAB/Octave stores `ld_r` as native sparse double loaded from `.mat`.
Python retains `a1freq` as `float32` from the `.npz` payload; MATLAB/Octave
returns `a1freq` as double.

## API

```text
load_ld(path, reference, optional default_chrX_sex) -> LDPanel
validate_ld_distribution(path, optional check_payload_structure) -> report

LDPanel.a1freq(optional chrX_sex) -> num_snp float vector
LDPanel.default_chrX_sex -> "female" | "male" | "combined"
LDPanel.select_shards(shards) -> LDPanel
LDPanel.multiply_r2(M, optional chrX_sex) -> vector or matrix with same shape as M
fast_prune(logpvec, ld_panel, optional r2_threshold, optional chrX_sex) -> logpvec
```

`path` may identify a panel root containing `ld_manifest.json` or a single LD
shard file. `reference` is required. `load_ld` must validate that every loaded
LD shard is compatible with the corresponding `ReferenceShard`: shard label,
`num_snp`, and `reference_checksum` must match. A mismatch is an error.

For a panel-root path, `load_ld` derives expected LD files from the supplied
`ReferencePanel` and resolves them through the manifest. Missing required LD
files are errors.

`validate_ld_distribution` is an explicit distribution-QA path, not part of
default loading. It validates manifest file MD5 checksums and manifest/per-file
metadata agreement. When `check_payload_structure` is true, it may also perform
expensive payload checks such as sparse index bounds, explicit diagonal, and
symmetry validation.

LD-specific public cache APIs are not part of the LD contract. There is no
`save_ld_cache` or `load_ld_cache`; the runtime distribution artifacts are the
load-efficient representation.

Expected behavior:

- Panel accessors are read-only and concatenate shard data in reference panel
  order.
- `default_chrX_sex` defaults to `"female"` when omitted.
- `r2_threshold` defaults to `0.2` when omitted.
- `LDPanel.a1freq(optional chrX_sex)` uses `default_chrX_sex` for chrX by
  default. `chrX_sex` applies only to chrX; chr1-22 always use their
  sex-agnostic shard. If chrX is present, an override must name a chrX shard
  present in the panel. If chrX is absent, `chrX_sex` is ignored.
- `LDPanel.multiply_r2(M, optional chrX_sex)` computes `LD_r² * M`
  independently per reference shard, never materializes a dense genome-wide LD
  matrix, accepts both 1-D and 2-D inputs, and returns the same shape as `M`.
- Python should return `float32` for `float32` input and `float64` for
  `float64` input when multiplying by LD. MATLAB/Octave may return double
  because the stored MATLAB LD matrix is sparse double. Cross-runtime tests
  compare numerical values with dtype-aware tolerances, not dtype identity.
- `fast_prune` applies greedy significance-based pruning independently per
  reference shard using the LD selected by `chrX_sex`.
- Shard-level matrix multiplication is an internal implementation detail, not
  part of the public API.

## r² Matrix Multiply

The core LD operation multiplies the per-shard squared LD matrix by a
genome-wide matrix or vector. This is the only numerical primitive exposed by
`LDPanel`; all downstream computations (LD scores, annotation LD weighting,
per-SNP variance) are expressed in terms of it.

```text
LDPanel.multiply_r2(M, optional chrX_sex) -> vector or matrix with same shape as M
```

`M` is a `num_snp x k` matrix or a `num_snp` vector aligned to the reference
panel row order. The result has the same shape as `M`.

`LDPanel.multiply_r2(M, optional chrX_sex)` splits `M` by reference shard,
computes `LD_r² * M_shard` for each shard, and concatenates the results in
reference panel order. Omitting `chrX_sex` uses `default_chrX_sex` for the chrX
shard. Pass `"female"`, `"male"`, or `"combined"` to override in
language-specific syntax. The `chrX_sex` selector affects only chrX; autosomes
are always sex-agnostic. If chrX is absent from the panel, `chrX_sex` is
ignored.

The operation must not materialize a dense genome-wide LD matrix. Callers
pre-exclude SNPs by setting the corresponding rows of `M` to zero.

## Pruning

```text
fast_prune(logpvec, ld_panel, optional r2_threshold, optional chrX_sex) -> logpvec
```

Greedy significance-based LD pruning. Input and output are genome-wide vectors
aligned to the reference panel.

Algorithm (applied independently per reference shard):

1. Stable-sort SNPs by `|logpvec|` descending; ties keep reference order. Skip
   `NaN` entries (treated as pre-excluded).
2. For each SNP in sorted order: if not already pruned, retain the current SNP
   and mark all other SNPs with `r² >= r2_threshold` as pruned (output set to
   `NaN`). The current SNP is never pruned by its own diagonal.
3. Retained positions keep their original `logpvec` values.

Omitting `chrX_sex` uses `default_chrX_sex` for the chrX shard; if chrX is
absent, the selector is ignored. Omitting `r2_threshold` uses `0.2`.
