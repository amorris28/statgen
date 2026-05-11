# LD Distribution Format: Runtime-Native Proposal v2

## Summary

This document proposes revising the LD specification around a runtime-native
distribution model:

- Python users receive NumPy/SciPy `.npz` LD files.
- MATLAB/Octave users receive native `.mat` LD files containing preconstructed
  sparse matrices.
- The Python `.npz` format is also the documented build handoff format from
  `statgen_build_ld.py` into the MATLAB/Octave converter.
- LD-specific cache APIs are removed from the public contract.

This is a deliberate exception to the current package-wide "portable source
plus optional cache" model. LD panels are derived, huge, sparse numerical
objects whose practical usefulness depends on load and matrix-operation
performance. For LD, the distribution artifact should be the load-efficient
runtime representation, not a slower portable source that every user must
convert locally.

Backward compatibility is not a constraint for this change.

---

## 1. Design Decisions

The proposal is based on the following decisions:

1. **LD stores signed correlations `r`, not `r²`.**
   Consumers that need `r²` compute it from signed `r`. This preserves the
   existing semantic contract.

2. **Python distribution uses `.npz` with a CSC32 storage schema.**
   Each `.npz` file stores one symmetric sparse signed-`r` matrix in
   SciPy-compatible CSC components, with explicit unit diagonal and int32
   sparse indices. The Python in-memory representation is CSC.

3. **MATLAB/Octave distribution uses `.mat`.**
   Each `.mat` file stores one native sparse signed-`r` matrix, already
   constructed, symmetric, with explicit unit diagonal.

4. **Python `.npz` is the build handoff format.**
   `statgen_build_ld.py` writes `.npz` files. A MATLAB/Octave converter reads
   those `.npz` files and writes native sparse `.mat` files.

5. **The converter must work under Octave.**
   The `.mat` files produced by the Octave converter must be loadable by
   downstream MATLAB/Octave `statgen` code. MATLAB remains the normative runtime
   where MATLAB and Octave differ, but Octave compatibility is required for this
   build conversion path.

6. **chrX sex-specific files use explicit suffixes.**
   For example: `ld_chrX_female.npz`, `ld_chrX_male.npz`,
   `ld_chrX_combined.npz`, and corresponding `.mat` files.

7. **Panel roots include `ld_manifest.json`.**
   The manifest is authoritative for panel-root discovery and distribution-file
   integrity checks. Per-file metadata remains authoritative for shard-level
   reference compatibility.

8. **LD cache APIs are not part of the public LD API.**
   `load_ld` reads the runtime distribution artifact directly.

---

## 2. Motivation

The current spec defines a language-agnostic raw directory format:

```text
ld/chrN/
  metadata.json
  ld_idx1.i32
  ld_idx2.i32
  ld_r.f32
  mafvec.f32       ← historical name; semantically a1 frequency
```

This format is simple and portable, but it is not load-optimal. In particular,
MATLAB sparse matrices are efficient once constructed and efficient to load
from `.mat`, but constructing huge sparse matrices from COO-like triplets can
consume excessive CPU and peak RAM. For million-SNP LD matrices, relying on
MATLAB-only sparse construction at user load time is not acceptable.

A user-side cache conversion model puts the expensive conversion on each user:

```text
portable LD triplets -> load_ld -> LDPanel -> save_ld_cache
```

For LD, this is the wrong boundary. The expensive construction should happen
once during reference-panel distribution, not on every user machine.

---

## 3. Proposed Architecture

```text
PLINK bfile
    |
    v
statgen_build_ld.py
    |
    v
Python LD distribution / build handoff
ld_chr1.npz
ld_chr2.npz
...
ld_chrX_female.npz
ld_chrX_male.npz
    |
    v
MATLAB/Octave converter
    |
    v
MATLAB LD distribution
ld_chr1.mat
ld_chr2.mat
...
ld_chrX_female.mat
ld_chrX_male.mat
```

End users load the artifact for their runtime:

```text
Python:       load_ld(ld_dir_or_file, reference, ...) -> LDPanel
MATLAB:       load_ld(ld_dir_or_file, reference, ...) -> LDPanel
MATLAB/Octave load_ld reads .mat, not .npz, on the user-facing path.
```

The Python `.npz` files serve two roles:

- Python runtime distribution format.
- Documented handoff format into the MATLAB/Octave converter.

The MATLAB `.mat` files are not caches in the old sense. They are the canonical
MATLAB/Octave LD distribution artifacts.

---

## 4. Logical LD Contract

Both runtime formats must represent the same logical LD shard:

- `chr`: shard label, one of `1`-`22` or `X`;
- `sex`: `null` for autosomes, one of `"female"`, `"male"`, `"combined"` for
  chrX;
- `num_snp`: number of SNPs in the paired reference shard;
- `ld_r`: symmetric sparse signed correlation matrix, shape
  `num_snp x num_snp`;
- diagonal entries of `ld_r`: explicit `1.0`;
- off-diagonal entries: signed `r` oriented to `a1` dosage;
- missing off-diagonal entries: interpreted as zero;
- `a1freq`: float vector of `a1` allele frequencies aligned to the paired
  `ReferenceShard`;
- `reference_checksum`: MD5 checksum of the paired reference shard;
- build metadata: PLINK command, PLINK version if available, window size,
  `r²` threshold, sample count, and schema version.

The paired `ReferencePanel` and `LDPanel` must originate from the same BIM row
order and allele orientation. The LD build additionally depends on the paired
BED/FAM genotype source used to estimate LD and allele frequencies.

The logical contract is cross-runtime. The physical byte layout is
runtime-specific.

`LDPanel.multiply_r2(M, optional chrX_sex)` computes with `ld_r .^ 2`, not
with signed `ld_r` directly.

`multiply_r2` should return the same floating dtype/class as `M` where the
runtime can do so without changing the stored LD representation. Python should
return `float32` for `float32` input and `float64` for `float64` input.
MATLAB/Octave may return double because the stored MATLAB LD matrix is sparse
double. Cross-runtime tests compare numerical values with dtype-aware
tolerances, not dtype identity.

`fast_prune(logpvec, ld_panel, optional r2_threshold, optional chrX_sex)`
returns a floating vector with the same shape as `logpvec`. Retained SNPs keep
their original `logpvec` values; pruned SNPs are set to `NaN`. Input `NaN`
entries are treated as pre-excluded and remain `NaN` in output. Python preserves
`float32` and `float64` inputs and promotes non-floating inputs to `float64` so
pruned values can be represented as `NaN`. MATLAB/Octave returns double.

---

## 5. Python `.npz` Distribution and Handoff Format

Each `.npz` file represents one LD shard. The sparse matrix is stored in CSC
form because this maps naturally to MATLAB's column-oriented sparse format
during conversion.

The public object remains `LDShard`; the `32` qualifier belongs to the storage
format, not to the logical object name.

Required arrays:

```text
data       float32, length nnz
indices    int32, length nnz, zero-based row indices
indptr     int32, length num_snp + 1, zero-based CSC column pointers
shape      int64,   length 2, equals [num_snp, num_snp]
a1freq     float32, length num_snp
metadata   uint8,   UTF-8 JSON bytes
```

`data`, `indices`, and `indptr` encode a symmetric sparse signed-`r` matrix
with explicit diagonal. `indices` are zero-based because `.npz` is a Python /
SciPy-oriented format and because zero-based sparse coordinates are already the
portable coordinate convention in the current spec.

The `.npz` archive must not contain object arrays or pickled payloads.

`metadata` is a UTF-8 JSON document encoded as a `uint8` vector, not as a NumPy
string/object array. This keeps the converter simple and avoids NumPy object
loading.

Example metadata:

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
  "dtype": {
    "data": "float32",
    "indices": "int32",
    "indptr": "int32",
    "a1freq": "float32"
  },
  "reference_checksum": "...",
  "build_tool": "plink2",
  "build_command": "plink2 ...",
  "ld_window_kb": 10000,
  "ld_r2_threshold": 0.05,
  "num_sample": 1000
}
```

Artifact builders and the MATLAB/Octave converter must validate before writing
published artifacts:

- required array names are present;
- `shape == [num_snp, num_snp]`;
- `len(indptr) == num_snp + 1`;
- `len(data) == len(indices) == indptr[-1]`;
- `len(a1freq) == num_snp`;
- `indices` and `indptr` are int32;
- `nnz < 2^31`;
- index bounds are valid;
- metadata schema version is supported;
- metadata and array dimensions agree.

Runtime loaders perform lightweight compatibility gates before loading large
payloads: schema/format checks, manifest/per-file metadata agreement,
compatibility with the supplied `ReferencePanel`, and constant-time shape/count
checks such as `shape`, `num_snp`, `nnz`, `len(indptr)`, and `indptr[-1]`.
Runtime loaders may assume the sparse payload was validated at
artifact-build time. They are not required to perform full sparse structural
validation, diagonal scans, symmetry scans, or file-hash validation on the
default load path.

The Python loader may construct a SciPy `csc_matrix` directly from these arrays.
It should avoid dense materialization.

Python `LDShard` stores the matrix in CSC form in memory. Operations that would
normally prefer CSR should use the symmetry of `ld_r` to work efficiently with
CSC column access rather than converting the matrix to CSR as part of the
object representation.

---

## 6. MATLAB `.mat` Distribution Format

Each `.mat` file represents one LD shard and stores a native sparse matrix.
The MATLAB distribution format is MAT-file v7.3. The legacy MATLAB v5 MAT-file
format is not acceptable for LD distribution artifacts because realistic LD
shards exceed its large-variable/file-size limits. The Octave converter must
produce files loadable by downstream MATLAB/Octave code with equivalent
large-sparse-variable behavior; if Octave cannot write compatible files for a
target panel, MATLAB must be used for the distribution conversion step.

Required variables:

```text
ld_r       sparse double, num_snp x num_snp
a1freq     numeric vector, length num_snp
metadata   struct
```

`ld_r` must be:

- native MATLAB sparse;
- symmetric;
- signed `r`;
- explicit unit diagonal;
- aligned to the paired reference shard.

`metadata` must contain the same logical fields as the `.npz` metadata:

```text
object_type
schema_version
format
chr
sex
num_snp
nnz
matrix
diagonal
value
reference_checksum
build_tool
build_command
ld_window_kb
ld_r2_threshold
num_sample
```

The exact MATLAB struct field names must be valid MATLAB identifiers.

The MATLAB/Octave loader should load `ld_r` directly and must not reconstruct a
sparse matrix from COO/CSC arrays on the default user path.

---

## 7. MATLAB/Octave Converter

The converter reads Python `.npz` files and writes MATLAB `.mat` files:

```text
convert_ld_npz_to_mat(npz_path, mat_path)
```

or equivalent batch API:

```text
convert_ld_panel_npz_to_mat(npz_dir, mat_dir)
```

The converter must work in Octave for testability and build reproducibility.

Because `.npz` is a ZIP archive of `.npy` arrays, the converter needs only a
limited reader for the exact `.npy` subset written by `statgen_build_ld.py`:

- little-endian `float32`;
- little-endian `int32`;
- `uint8`;
- C-order 1-D arrays.

The converter does not need to be a general NumPy reader.

Conversion steps:

1. Unzip the `.npz` archive into a temporary directory.
2. Read `data.npy`, `indices.npy`, `indptr.npy`, `shape.npy`, `a1freq.npy`,
   and `metadata.npy`.
3. Decode `metadata.npy` as UTF-8 JSON bytes.
4. Expand CSC column pointers into column indices.
5. Convert row and column indices from zero-based to one-based.
6. Validate metadata consistency before sparse construction.
7. Construct native sparse `ld_r = sparse(row, col, data, n, n)`.
8. Copy `reference_checksum` and `reference_bim` from the `.npz` metadata
   into the `.mat` metadata struct unchanged.
9. Save `ld_r`, `a1freq`, and `metadata` to `.mat`.

The sparse construction cost is paid once by the reference-panel maintainer,
not by end users.

The converter assumes the `.npz` payload was produced and validated by
`statgen_build_ld.py`. It is not required to perform a full O(nnz) symmetry
scan or diagonal scan during conversion.

---

## 8. File Naming and Panel Layout

Recommended filenames:

```text
ld_chr1.npz
ld_chr2.npz
...
ld_chr22.npz
ld_chrX_female.npz
ld_chrX_male.npz
ld_chrX_combined.npz
```

and:

```text
ld_chr1.mat
ld_chr2.mat
...
ld_chr22.mat
ld_chrX_female.mat
ld_chrX_male.mat
ld_chrX_combined.mat
```

Autosomal shards use no sex suffix. chrX shards always use an explicit sex
suffix.

Not all chrX sex labels need to be present. `female` and `male` remain the
default sex-specific build outputs. `combined` is optional and
assumption-dependent.

The Python and MATLAB/Octave distributions are separate directory trees, each
with its own `ld_manifest.json`. A Python distribution root contains `.npz`
shard files and a manifest with `runtime_format: "python_npz_csc32"`. A
MATLAB/Octave distribution root contains `.mat` shard files and a manifest with
`runtime_format: "matlab_mat_sparse_double"`. Bundled reference `.bim` files
are present in both roots and are byte-identical.

A panel root is a directory containing one or more LD shard files, bundled
`.bim` files, and `ld_manifest.json`. A loader may also accept a single shard
file.

For panel-root loads, the directory must contain `ld_manifest.json`. The
manifest is authoritative for file discovery; loaders must not infer the panel
by globbing arbitrary filenames.

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
    },
    {
      "chr": "X",
      "sex": "female",
      "file": "ld_chrX_female.npz",
      "file_md5": "...",
      "num_snp": 123,
      "nnz": 4567,
      "reference_checksum": "..."
    }
  ]
}
```

`runtime_format` is `"python_npz_csc32"` for Python distributions and
`"matlab_mat_sparse_double"` for MATLAB/Octave distributions. The manifest
records distribution-file MD5 checksums (`file_md5`) so users or distribution
tools can check for missing, truncated, or wrong files. Runtime loaders are not
required to compute file MD5 checksums on the default load path. Each shard
file still carries its own metadata, including `reference_checksum`; loaders
must validate that manifest entries and per-file metadata agree.

Manifest entries and per-file metadata must agree on `chr`, `sex`, `num_snp`,
`nnz`, `reference_checksum`, and runtime/storage format. `file` and `file_md5`
are manifest-only fields.

For panel-root loads:

- autosomal expected files are derived from the supplied `ReferencePanel`
  shard labels and resolved through the manifest;
- chrX discovery checks the explicit allowed suffixes
  `female`, `male`, `combined`, as listed in the manifest;
- files present in the directory but absent from the manifest are ignored.

If the supplied `ReferencePanel` is a shard subset, `load_ld` loads only the LD
files required for those reference shards. LD files for other reference shards
may be present in the directory and are ignored.

---

## 9. Public API Revision

LD should expose one user-facing load path:

```text
load_ld(path, reference, optional default_chrX_sex) -> LDPanel
validate_ld_distribution(path, optional check_payload_structure=false) -> report
statgen_create_ld_manifest.py --ld <root>               creates Python ld_manifest.json
create_ld_mat_manifest(npz_root, mat_root, shards)       creates MATLAB/Octave ld_manifest.json

LDPanel.a1freq(optional chrX_sex) -> num_snp float vector
LDPanel.default_chrX_sex -> "female" | "male" | "combined"
LDPanel.select_shards(shards) -> LDPanel
LDPanel.multiply_r2(M, optional chrX_sex) -> vector or matrix with same shape as M
fast_prune(logpvec, ld_panel, optional r2_threshold, optional chrX_sex)
    -> logpvec with same shape; pruned SNPs set to NaN, retained SNPs unchanged
```

`load_ld` reads the runtime distribution artifact:

- Python reads `.npz`.
- MATLAB/Octave reads `.mat`.

`validate_ld_distribution` is an explicit distribution-QA path, not part of
default loading. It validates manifest file MD5 checksums and manifest/per-file
metadata agreement. When `check_payload_structure=true`, it may also perform
expensive payload checks such as sparse index monotonicity, index bounds,
explicit diagonal, and symmetry validation.

`reference` is required. `load_ld` must validate that every loaded LD shard is
compatible with the corresponding `ReferenceShard`: shard label, `num_snp`, and
`reference_checksum` must match. A mismatch is an error at load time.

For a panel-root path, `load_ld` derives expected LD files from the supplied
`ReferencePanel`; missing required LD files are errors. If `reference` contains
a subset of reference shards, only those shards are loaded.

The contract is that the `ReferencePanel` and `LDPanel` originate from the same
BIM file rows. The LD panel additionally records the BED/FAM-derived LD build
metadata.

LD-specific public cache APIs are removed from the LD contract:

```text
save_ld_cache(panel, path)     removed
load_ld_cache(path, ...)       removed
```

`default_chrX_sex` defaults to `"female"` when omitted. `load_ld` loads all
chrX LD shards listed in the manifest when the supplied `ReferencePanel`
contains chrX. If the supplied `ReferencePanel` has no chrX shard, `load_ld`
does not load or retain chrX LD shards even if chrX files are present in the LD
directory.

If chrX is loaded, the effective `default_chrX_sex` must name a loaded chrX
shard. A combined-only panel must therefore be loaded with
`default_chrX_sex="combined"`. If chrX is not loaded because the supplied
`ReferencePanel` has no chrX, `default_chrX_sex` is retained as metadata and
ignored by chrX-unrelated operations.

---

## 10. Empirical Size and Sparse Index Width

The following empirical panels use symmetric LD matrices with explicit
diagonal. The reported `nnz` is already the final stored sparse count, including
upper triangle, lower triangle, and diagonal.

Memory estimates are sparse payload estimates only. They exclude object
overhead, decompression buffers, matrix-operation temporaries, input/output
vectors, and converter peak-memory overhead.

Python estimates assume CSC with `float32 data`:

- int32 index layout: `float32 data + int32 indices + int32 indptr`;
- int64 index layout: `float32 data + int64 indices + int64 indptr`.

MATLAB estimates assume native sparse double:

- `double values + MATLAB sparse index storage`;
- MATLAB sparse single is out of scope.

Using int32 sparse indices in Python reduces the Python sparse payload from
about 12 bytes per stored entry to about 8 bytes per stored entry. This means
the int32 representation uses roughly two thirds of the int64 representation,
or saves roughly one third of the sparse payload RAM, when `nnz < 2^31`.

### HRC EUR QC

| Chr | # variants | nnz | Python int32 GiB | Python int64 GiB | MATLAB double GiB |
|---:|---:|---:|---:|---:|---:|
| 1 | 947,323 | 283,872,190 | 2.12 | 3.18 | 4.24 |
| 2 | 1,013,059 | 345,401,163 | 2.58 | 3.87 | 5.15 |
| 3 | 851,931 | 319,942,501 | 2.39 | 3.58 | 4.77 |
| 4 | 868,251 | 320,512,054 | 2.39 | 3.59 | 4.78 |
| 5 | 774,184 | 288,002,405 | 2.15 | 3.22 | 4.30 |
| 6 | 778,335 | 461,261,998 | 3.44 | 5.16 | 6.88 |
| 7 | 698,208 | 249,942,342 | 1.86 | 2.80 | 3.73 |
| 8 | 673,056 | 284,817,738 | 2.12 | 3.19 | 4.25 |
| 9 | 520,727 | 149,278,817 | 1.11 | 1.67 | 2.23 |
| 10 | 608,173 | 214,328,926 | 1.60 | 2.40 | 3.20 |
| 11 | 592,848 | 330,855,843 | 2.47 | 3.70 | 4.93 |
| 12 | 571,483 | 218,096,827 | 1.63 | 2.44 | 3.25 |
| 13 | 437,449 | 149,004,313 | 1.11 | 1.67 | 2.22 |
| 14 | 393,607 | 125,384,632 | 0.94 | 1.40 | 1.87 |
| 15 | 344,991 | 95,631,865 | 0.71 | 1.07 | 1.43 |
| 16 | 383,058 | 91,721,098 | 0.68 | 1.03 | 1.37 |
| 17 | 321,325 | 77,314,477 | 0.58 | 0.87 | 1.15 |
| 18 | 342,710 | 90,774,590 | 0.68 | 1.02 | 1.36 |
| 19 | 263,086 | 74,304,301 | 0.55 | 0.83 | 1.11 |
| 20 | 270,782 | 64,408,350 | 0.48 | 0.72 | 0.96 |
| 21 | 163,979 | 36,756,460 | 0.27 | 0.41 | 0.55 |
| 22 | 161,946 | 36,371,718 | 0.27 | 0.41 | 0.54 |
| **Total** | **11,980,511** | **4,307,984,608** | **32.14** | **48.23** | **64.28** |

### UKB Imp V3 QC

| Chr | # variants | nnz | Python int32 GiB | Python int64 GiB | MATLAB double GiB |
|---:|---:|---:|---:|---:|---:|
| 1 | 1,008,419 | 385,817,100 | 2.88 | 4.32 | 5.76 |
| 2 | 1,102,083 | 475,112,369 | 3.54 | 5.32 | 7.09 |
| 3 | 932,470 | 447,064,342 | 3.33 | 5.00 | 6.67 |
| 4 | 944,938 | 443,339,396 | 3.31 | 4.96 | 6.61 |
| 5 | 851,028 | 408,517,599 | 3.05 | 4.57 | 6.09 |
| 6 | 864,841 | 825,900,174 | 6.16 | 9.24 | 12.31 |
| 7 | 759,403 | 354,703,715 | 2.65 | 3.97 | 5.29 |
| 8 | 725,143 | 384,044,264 | 2.86 | 4.30 | 5.73 |
| 9 | 562,514 | 201,616,696 | 1.50 | 2.26 | 3.01 |
| 10 | 656,697 | 306,436,181 | 2.29 | 3.43 | 4.57 |
| 11 | 642,215 | 467,757,577 | 3.49 | 5.23 | 6.97 |
| 12 | 610,598 | 300,219,342 | 2.24 | 3.36 | 4.48 |
| 13 | 467,348 | 202,077,822 | 1.51 | 2.26 | 3.01 |
| 14 | 415,551 | 173,674,399 | 1.30 | 1.94 | 2.59 |
| 15 | 363,625 | 131,685,827 | 0.98 | 1.47 | 1.96 |
| 16 | 399,314 | 120,109,519 | 0.90 | 1.35 | 1.79 |
| 17 | 342,283 | 114,459,758 | 0.85 | 1.28 | 1.71 |
| 18 | 364,669 | 123,293,388 | 0.92 | 1.38 | 1.84 |
| 19 | 285,517 | 109,132,239 | 0.81 | 1.22 | 1.63 |
| 20 | 285,440 | 90,350,343 | 0.67 | 1.01 | 1.35 |
| 21 | 172,405 | 50,072,876 | 0.37 | 0.56 | 0.75 |
| 22 | 170,168 | 50,863,956 | 0.38 | 0.57 | 0.76 |
| **Total** | **12,926,669** | **6,166,248,882** | **45.99** | **69.01** | **91.98** |

The empirical per-chromosome `nnz` values are all below signed int32 capacity,
even though genome-wide totals exceed 4 billion. Since LD is distributed
per chromosome, Python `.npz` builders should use int32 sparse indices for
these panels.

---

## 11. Large `nnz` and Future CSC64

A chromosome-level LD matrix can have more than 2 billion stored entries after
symmetrization and explicit diagonal insertion, and may plausibly exceed
`2^32` entries. The current Python distribution format is intentionally CSC32:
it requires int32 sparse indices and `nnz < 2^31` for each LD file.

If a future panel has a chromosome shard whose `nnz` does not fit in int32,
that should be introduced as an explicit new storage schema, for example
`statgen_ld_npz_csc64`, rather than silently widening the v1 CSC32 format. The
logical object can remain `LDShard`; the bit width belongs to the storage
schema.

MATLAB sparse matrices are column-compressed and can address dimensions beyond
32-bit limits on 64-bit platforms, but practical limits are dominated by RAM,
temporary construction memory, and `.mat` load/save behavior. A single native
sparse chromosome matrix with billions of entries may be operationally too
large even if it is theoretically representable.

---

## 12. Required Spec Changes

### `spec/SPEC.md`

Revise the portable storage policy. The current rule says portable object files
should be language-agnostic and explicitly lists raw LD triplets. Replace this
with an LD-specific exception:

- most canonical source objects remain language-agnostic;
- LD is a derived, large, performance-bound object;
- LD distribution artifacts are runtime-native and schema-versioned;
- cross-runtime byte identity is not required;
- cross-runtime logical equivalence is required.

### `spec/ld.md`

Replace the current raw directory "Disk representation" section with:

- logical LD contract;
- Python `.npz` distribution / handoff format;
- MATLAB `.mat` distribution format;
- file naming and chrX suffix rules;
- required `ld_manifest.json` for panel-root discovery and file MD5 sums;
- required `ReferencePanel` compatibility checks in `load_ld`;
- CSC32 index rules and future CSC64 extension boundary;
- converter requirements;
- revised API without LD cache functions.

Keep signed `r` as the semantic LD value.

### `spec/performance-contract.md`

Clarify that LD distribution artifacts are not optional caches. They are
runtime-native distribution formats by design.

Cache rules continue to apply to object families that still define cache APIs.

### `spec/testing.md`

Update LD fixture expectations:

- checked-in Python `.npz` LD fixtures are canonical build-handoff fixtures;
- checked-in or generated `.mat` fixtures must be produced by the
  MATLAB/Octave converter;
- Python and MATLAB/Octave tests compare logical behavior, not byte identity;
- tests must cover `multiply_r2`, `fast_prune`, `a1freq`, chrX sex selection,
  manifest validation, schema validation, and reference checksum compatibility.

### `dev/implementation_plan.md`

Update LD phases:

- implement Python `.npz` LD writer/loader;
- implement limited MATLAB/Octave `.npz`/`.npy` reader for conversion only;
- implement `.npz` to `.mat` converter;
- implement MATLAB/Octave `.mat` LD loader;
- remove LD cache API work.

---

## 13. Testing Requirements

Tests should verify:

- Python can load `.npz` LD fixtures directly.
- Octave converter can read `.npz` and write `.mat`.
- MATLAB/Octave `load_ld` can load converter-produced `.mat`.
- Python and MATLAB/Octave agree on:
  - shard labels;
  - `a1freq`;
  - reference checksums;
  - load-time compatibility against the supplied `ReferencePanel`;
  - chrX default-sex behavior;
  - `multiply_r2` results within tolerance;
  - `fast_prune` results.
- Invalid `.npz` files fail clearly:
  - missing arrays;
  - wrong dtype;
  - inconsistent `indptr`;
  - out-of-range indices;
  - non-int32 sparse index arrays;
  - `nnz >= 2^31`;
  - malformed metadata;
  - unsupported schema version.
- Invalid manifests fail clearly:
  - missing required shard entries;
  - wrong runtime format;
  - missing files;
  - manifest metadata disagreement with per-file metadata.
- Invalid `.mat` files fail clearly:
  - missing `ld_r`;
  - dense `ld_r`;
  - wrong dimensions;
  - malformed metadata;
  - unsupported schema version.

The tests should not require Python and MATLAB/Octave cache files to match
byte-for-byte.

---

## 14. Remaining Open Questions

The major design decisions are settled. Remaining details are implementation
level:

1. If a future reference panel requires `nnz >= 2^31` in one LD file, should the
   project add a separate `statgen_ld_npz_csc64` storage schema?

This is a future schema-extension question outside the v1 CSC32 contract.
