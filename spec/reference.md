# ReferenceShard and ReferencePanel

## Represents

A reference panel is a genome-wide table of SNP positions and alleles that
defines the coordinate system for all other objects. Each SNP is described by
chromosome, base-pair position, identifier, and an ordered allele pair under
the input contract in [SPEC.md](SPEC.md):

- `chr`: upstream contig label, expected to use the selected `genomatch`
  contig naming mode, normally NCBI labels such as `1`-`22` and `X`;
- `snp`: SNP identifier;
- `bp`: base-pair position;
- `a1`, `a2`: ordered alleles under the allele contract in
  [SPEC.md](SPEC.md); `a1` is the non-reference allele, `a2` is
  the reference allele.

Optional columns include `source` and `variant_id`.

## Disk representation

Canonical disk representation is PLINK-style `.bim`:

```text
chr  snp  cm  bp  a1  a2
```

The `cm` column is a legacy placeholder required only for 6-column schema
compatibility. Loaders MUST validate that the `cm` column is present; `cm`
MUST NOT be stored in in-memory objects and MUST NOT participate in checksums,
compatibility checks, caches, or public APIs.

A reference panel may be sharded or non-sharded on disk:

- **Sharded**: one `.bim` file per chromosome via an `@` placeholder in the
  path (e.g. `chr@.bim`).
- **Non-sharded**: a single `.bim` file spanning all chromosomes.

Shard row order defines shard-local SNP indices. Panel row order defines global
SNP indices.

## In-memory objects

A `ReferenceShard` is the in-memory representation of one `.bim` shard.
A `ReferencePanel` is an ordered collection of `ReferenceShard` objects and
defines the multi-chromosome SNP order used for all aligned in-memory objects.
Each shard carries the MD5 reference checksum defined below. A
`ReferenceShard` contains exactly one chromosome label: all loaded BIM rows in
the shard must have `chr == ReferenceShard.label`. Shard constructors/loaders
must validate this invariant.

## Reference identity keys

Each `ReferenceShard` has a shard-local reference checksum. The checksum input
is the UTF-8 text sequence:

```text
chr:bp:a1:a2\n
```

one line per variant in shard row order. The `chr` value is the contig label as
loaded; no normalization is applied. `bp` is the decimal integer base-pair
coordinate; `a1`/`a2` are the allele strings as loaded. The algorithm is MD5
over that byte sequence, encoded as a lowercase hexadecimal string. The checksum
identifies row order, contig naming, position, and allele orientation for
alignment checks; it is not a security primitive.

Each variant also has fixed-width allele identity hashes for
sumstats-to-reference matching, which is the fast object-alignment path that
needs exact variant identity. Chromosome identity is supplied by the shard
label, and base-pair position is already stored as a numeric vector, so only the
variable-length allele strings need compact hashed representations.

`a1_hash64` and `a2_hash64` are `uint64` vectors, length `num_snp`, computed
from allele strings as loaded. If any allele exceeds 150 characters,
implementations must warn and compute the hash from only the first 150
characters. The hash algorithm is a deterministic two-lane polynomial hash over
the UTF-8 bytes of the allele string prefix, chosen so all intermediate
arithmetic is exactly representable in MATLAB/Octave `double` without relying
on unsigned integer overflow:

```text
p = 2147483647  # 2^31 - 1, prime
base1 = 257
base2 = 263
h1 = 1
h2 = 1
for byte in allele_utf8_bytes(first_150_characters(allele)):
    x = double(byte) + 1
    h1 = mod(h1 * base1 + x, p)
    h2 = mod(h2 * base2 + x, p)
allele_hash64 = bitshift(uint64(h1), 32) + uint64(h2)
```

Because `h1 < p`, `h2 < p`, and both bases are small, all products and sums in
the recurrence are below `2^53` and can be computed exactly in IEEE-754 double
precision. Implementations must not use saturating `uint64` multiplication to
compute allele hashes.

Sumstats-to-reference matching uses the tuple `(shard label, bp, a1_hash64,
a2_hash64)`. The allele hashes are not used for reference compatibility checks;
compatibility continues to use the stored per-shard `reference_checksum`. The
allele hashes are not security primitives. Collision risk is negligible for
non-adversarial allele strings in this use; implementations should still fail
clearly on duplicate matching keys within a shard or source input when such
duplicates would make alignment ambiguous.

MATLAB/Octave sumstats loaders should match within each shard using a native
numeric sort/merge over `(bp, a1_hash64, a2_hash64)`. Equivalent native numeric
implementations are allowed, but implementations must not reconstruct string
join keys and must not cast `uint64` hashes to `double` for matching because
values above `2^53` would lose precision. Duplicate matching keys within a
shard or source input must fail clearly when they would make alignment
ambiguous.

Regardless of how the panel is sharded, callers access genome-wide vectors and
matrices through read-only panel-level accessors that concatenate across shards
transparently.

## Representation

Implementations may use language-native containers. They must preserve shard
order, expose the required vectors, retain per-shard checksums, and provide
zero-based shard offsets for compatibility with portable metadata.

## Cache layout

MATLAB/Octave reference caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
chr
snp
bp
a1
a2
a1_hash64
a2_hash64
```

`metadata` is a struct with:

```text
schema = "reference_cache/0.1"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
```

`chr`, `snp`, `a1`, and `a2` are panel-wide cell arrays of strings. `bp` is a
panel-wide numeric vector. `a1_hash64` and `a2_hash64` are panel-wide `uint64`
vectors. The BIM `cm` field is not cached. Shard offsets are zero-based
half-open intervals into the panel-wide variables and are sufficient to
reconstruct `ReferenceShard` objects. Cache metadata validation should be cheap,
depending on shard count and array dimensions rather than scanning all SNP
values. `load_reference_cache` trusts stored shard checksums, matching the
cache behavior of reference-aligned objects that cannot recompute those
checksums themselves. Callers who want to verify reference cache integrity may
explicitly call `ReferencePanel.validate_checksums()`.

## MATLAB/Octave thin cache loading

The `a1_hash64`/`a2_hash64` cache fields and sumstats matching semantics are
cross-runtime object contracts. Thin reference-cache loading is
MATLAB/Octave-specific implementation behavior, motivated by the high cost of
deserializing large `.mat` cell arrays of strings.

MATLAB/Octave `load_reference_cache` defaults to a thin load. The thin load
reads only `metadata`, `bp`, `a1_hash64`, and `a2_hash64`; it does not deserialize
the cached `chr`, `snp`, `a1`, or `a2` string cell arrays. A full load is
requested explicitly:

```text
load_reference_cache(path, optional shards, full=true) -> ReferencePanel
```

or, in MATLAB/Octave name-value syntax:

```matlab
reference = statgen.load_reference_cache(path, 'full', true)
reference = statgen.load_reference_cache(path, shards, 'full', true)
```

Thin-loaded references synthesize `ReferencePanel.chr` lazily on access,
without caching the synthesized vector. For each internal shard, synthesis
repeats that shard's label `num_snp` times; the panel accessor concatenates
those shard-level vectors in panel order. Thin-loaded `snp`, `a1`, and `a2`
accessors must fail clearly and instruct the caller to reload with `full=true`.
`bp`, `a1_hash64`, and `a2_hash64` are available in both thin and full loads.
Full loads may validate that saved `chr` agrees with shard labels.

## API

```text
load_reference(path, optional shards) -> ReferencePanel
save_reference_cache(panel, path, optional format)
load_reference_cache(path, optional shards, optional full=false) -> ReferencePanel

ReferencePanel.num_snp -> int
ReferencePanel.chr -> num_snp chromosome-label vector
ReferencePanel.snp -> num_snp string vector
ReferencePanel.bp  -> num_snp integer vector
ReferencePanel.a1  -> num_snp string vector
ReferencePanel.a2  -> num_snp string vector
ReferencePanel.a1_hash64 -> num_snp uint64 vector
ReferencePanel.a2_hash64 -> num_snp uint64 vector
ReferencePanel.shard_offsets -> table with shard_label, start0, stop0
ReferencePanel.select_shards(shards) -> ReferencePanel
ReferencePanel.is_object_compatible(object) -> bool
ReferencePanel.validate_checksums() -> bool
```

Expected behavior:

- Reference panels are external inputs; reference APIs never write `.bim`
  files. LD distribution builders (`statgen_build_ld.py` and the MATLAB/Octave
  converter) are the only tools that write `.bim` files, and only as exact
  copies of input BIM rows bundled alongside LD shards.
- Cache is a single file (non-sharded). The MATLAB/Octave cache layout is
  specified in the "Cache layout" section above.
- Shard discovery, contig validation, row-order validation, and shard subsetting
  follow [contigs-and-shards.md](contigs-and-shards.md).
- Each `ReferenceShard` must contain exactly one chromosome label. This is a
  global invariant of the object model, not just a cache-loader assumption.
- Compute and retain each shard reference checksum and each per-row
  `a1_hash64` and `a2_hash64` value.
- Accessors are read-only, concatenate shard columns in reference panel order,
  and return plain language-native vectors or tables.
- MATLAB/Octave thin-loaded cache references may omit `snp`, `a1`, and `a2`
  payloads in memory; those accessors must fail clearly unless the reference
  was loaded with `full=true`. Thin-loaded `chr` is synthesized lazily by
  repeating each shard label `num_snp` times and concatenating shard vectors in
  panel order, without caching. Thin-loaded `bp`, `a1_hash64`, and `a2_hash64`
  are available. Other runtimes may load all cache fields by default while
  preserving the same object contract.
- `shard_offsets` uses zero-based half-open intervals into genome-wide arrays.
- `shard_offsets.start0`/`stop0` are cross-language coordinate metadata, not
  direct language indices. MATLAB/Octave callers convert at use-site
  (`start0 + 1 : stop0`).
- `load_reference_cache` skips source-style row validation; `shards` subsetting
  applies against cached shard labels per
  [contigs-and-shards.md](contigs-and-shards.md).
- `load_reference_cache` validates shard-offset metadata and panel-wide vector
  lengths cheaply before reconstructing shard objects. MATLAB/Octave's default
  thin load validates only the loaded panel-wide payloads (`bp`, `a1_hash64`,
  and `a2_hash64`) plus metadata. A full load also validates full-field vector
  lengths. Cache loaders trust stored shard checksums and do not recompute them.
- `ReferencePanel.validate_checksums()` recomputes each shard checksum from the
  current `chr`, `bp`, `a1`, and `a2` values, compares it to the stored
  `ReferenceShard.checksum`, fails on mismatch with the shard label, and
  returns `true` on success. It requires full `a1` and `a2` fields and, in
  runtimes with thin loading, must fail clearly on thin-loaded references.
- `is_object_compatible` checks whether a loaded statgen object is aligned to
  this reference panel. Compatibility requires the same ordered shard labels,
  matching shard row counts, and matching shard checksums where available.
  Returns `true` when compatible, `false` otherwise. Implementations may log
  informational messages or warnings. Compatibility mismatches do not raise
  exceptions; callers branch on the returned boolean.
