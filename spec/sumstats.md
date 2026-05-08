# SumstatsShard and Sumstats

## Represents

GWAS summary statistics for one trait or source: p-values, optional Z scores,
sample sizes, effect estimates, and quality fields, all oriented to `a1` under
the allele contract in [SPEC.md](SPEC.md). Multiple traits are independent
objects that may share the same reference panel. A future multi-trait container
should use a different name, such as `SumstatsCollection`.

## Disk representation

Portable disk representation is one gzip-compressed TSV file per trait/source:

```text
TRAIT.tsv.gz
```

Canonical sumstats TSV files are produced by `genomatch` or equivalent
post-harmonization pipelines, not by `statgen`; `statgen` consumes these files
for loading, alignment, and caching.

Required columns:

- `chr` or `CHR`
- `bp` or `POS`
- `a1` or `EffectAllele`
- `a2` or `OtherAllele`
- `p` or `P`

Allele-column contract follows [SPEC.md](SPEC.md).

`p` is required and must be present for every source row. `p == 0` is allowed;
`p < 0`, `p > 1`, non-numeric, `NaN`, infinite, and missing `p` values are
validation errors.

Optional columns:

- `z`;
- `n`;
- `beta`;
- `se`;
- `eaf`;
- `info`.

Column recognition is case-insensitive. Apart from case, the only non-internal
column names accepted by `statgen` are the genomatch cleaned-sumstats names
`POS`, `EffectAllele`, and `OtherAllele`, which map to internal fields `bp`,
`a1`, and `a2`. Genomatch vmap-style `bp`, `a1`, and `a2` are accepted
directly. `SNP` may be present but is not used for matching.

The TSV does not have to contain every variant in the reference panel.
`load_sumstats` is a reference-driven projection from a possibly larger source.
Loading against a `ReferencePanel` projects rows into reference order by taking
validated source rows whose `chr` label is in the supplied reference shard set,
then matching within that shard on the tuple `(bp, a1_hash64, a2_hash64)`
defined in [reference.md](reference.md). Validated source rows outside the
supplied reference shard set are ignored for alignment. The allele hashes are
computed from each source row's exact `a1` and `a2` strings. This is
semantically an exact `chr:bp:a1:a2` join; the hashes are only fixed-width
implementation keys for sumstats-to-reference matching. The result is split
into `SumstatsShard`s matching the reference shards, and variants absent from
the TSV are represented as missing values. If a supplied reference shard has no
matching source rows at all, the loader warns and represents that shard as all
missing. The loader does not normalize or alias contig labels; matched sumstats
`chr` values must already match the reference labels.

## In-memory objects

A `SumstatsShard` is aligned to one `ReferenceShard`; a `Sumstats` object is an
ordered collection of `SumstatsShard` objects aligned to a `ReferencePanel`.
Each shard retains the paired reference checksum as `reference_checksum`.

Each `SumstatsShard` contains:

- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`;
- `logpvec`: float vector, length `num_snp`, containing `-log10(p)` with the
  zero-p convention described below;
- `zvec`: optional float vector, length `num_snp`, containing signed Z scores;
- `nvec`: optional float vector, length `num_snp`, containing effective sample
  size;
- `is_present`: derived logical vector, length `num_snp`, equal to
  `~isnan(logpvec)`;
- `beta_vec`: optional float vector, length `num_snp`, containing effect sizes
  oriented to `a1`;
- `se_vec`: optional float vector, length `num_snp`, containing standard errors;
- `eaf_vec`: optional float vector, length `num_snp`, containing effect allele
  frequency for `a1` where available;
- `info_vec`: optional float vector, length `num_snp`, containing imputation or
  variant-quality information where available.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: each shard's vector length and row order correspond to the
paired `ReferenceShard`, and any aligned cache is valid only for that
reference.

## Cache layout

MATLAB/Octave sumstats caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
logpvec
zvec
nvec
beta_vec
se_vec
eaf_vec
info_vec
```

`metadata` is a struct with:

```text
schema = "sumstats_cache/0.1"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
has_z
has_n
has_beta
has_se
has_eaf
has_info
```

`logpvec` is a required panel-wide numeric vector. Optional vectors are also
panel-wide numeric vectors when present; absent optional fields are saved as
`[]` and indicated by the corresponding `has_*` flag. There is no `has_p`
metadata field because `p` is represented by required `logpvec`. `is_present`
is derived from `logpvec` and is not stored. Shard offsets are zero-based
half-open intervals into the panel-wide vectors and are sufficient to
reconstruct `SumstatsShard` objects. `shard_checksums` stores the per-shard
reference checksums used to restore `SumstatsShard.reference_checksum`. Cache
metadata validation should be cheap, depending on shard count and vector
dimensions rather than scanning all SNP values. After structural validation
passes, cache loaders must perform the warning checks for optional
`zvec`/`nvec` completeness; those warning checks inspect vector values as
needed and are not part of the structural compatibility gate.

## Panel-level accessors

`Sumstats` exposes read-only genome-wide accessors that concatenate across
shards in reference panel order. The result is a plain array; downstream code
works with it natively without modifying the `Sumstats` object.

```text
Sumstats.logpvec  -> num_snp float vector
Sumstats.zvec     -> num_snp float vector, or missing optional field sentinel
Sumstats.nvec     -> num_snp float vector, or missing optional field sentinel
Sumstats.is_present -> num_snp logical vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec   -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec  -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
```

SNPs absent from the source file have `NaN` in `logpvec` and any present
optional vectors. Missing optional fields use the language-specific sentinel
defined in [SPEC.md](SPEC.md).

## API

```text
load_sumstats(path, reference) -> Sumstats
save_sumstats_cache(sumstats, path, optional format)
load_sumstats_cache(path, optional shards) -> Sumstats
create_sumstats(reference, pvec, optional zvec, optional nvec, optional beta_vec,
                optional se_vec, optional eaf_vec, optional info_vec) -> Sumstats

Sumstats.num_snp -> int
Sumstats.logpvec -> num_snp float vector
Sumstats.zvec -> num_snp float vector, or missing optional field sentinel
Sumstats.nvec -> num_snp float vector, or missing optional field sentinel
Sumstats.is_present -> num_snp logical vector
Sumstats.beta_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.se_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.eaf_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.info_vec -> num_snp float vector, or missing optional field sentinel
Sumstats.save_cache(path, optional format) -> void
Sumstats.select_shards(shards) -> Sumstats
```

Expected behavior:

- `reference` is required; rows are projected into reference order on load.
- `create_sumstats(...)` creates an in-memory `Sumstats` aligned to `reference`
  from full-panel vectors. Canonical sumstats TSV files remain external inputs.
- sumstats TSV column recognition is case-insensitive and only accepts the
  internal field names plus `POS`, `EffectAllele`, and `OtherAllele` from the
  genomatch cleaned-sumstats schema.
- sumstats-to-reference matching uses shard label plus `(bp, a1_hash64,
  a2_hash64)` from [reference.md](reference.md). Joins are exact after basic
  field parsing; the loader does not normalize chromosome labels, swap alleles,
  or perform strand handling. MATLAB/Octave implementations should use
  a native numeric sort/merge within each shard; equivalent native numeric
  approaches are allowed, but implementations must not reconstruct string join
  keys or cast `uint64` allele hashes to `double` for matching.
- required numeric field `p` must parse as a finite numeric value in the closed
  interval `[0, 1]` for every source row; non-numeric, `NaN`, infinite, missing,
  negative, or greater-than-one values are validation errors and must fail load
  with a clear message.
- optional numeric fields `z` and `n` may be absent. `load_sumstats` warns when
  either column is absent. When present, each returns a full aligned vector with
  `NaN` for unmatched rows and missing source values, and `load_sumstats` warns
  if any value is missing among aligned variants with `is_present == true`.
- cache is a single file (non-sharded). The MATLAB/Octave cache layout is
  specified in the "Cache layout" section above; absent optional fields are
  serialized as empty arrays (`[]`) so all field variables are always present
  in the `.mat` file, distinguished by the `has_*` flags in metadata.
- `load_sumstats_cache` performs cache-internal validation only and supports
  optional `shards` subsetting; per-shard checksums are trusted from cache
  metadata. It warns when cached `zvec` or `nvec` is absent, and warns when a
  present `zvec` or `nvec` has missing values among cached variants with
  `is_present == true`.
- missing variants are represented as `NaN` or masks; row order matches the
  reference panel.
- `logpvec` is derived from the required `p` column as `-log10(p)`. By
  convention, `p == 0` yields `Inf`.
- `is_present` is a derived accessor defined as `~isnan(logpvec)`. It is `true`
  for matched source rows with valid `p`, including `p == 0`, and `false` for
  reference variants absent from the source after alignment.
- for each optional field (`zvec`, `nvec`, `beta_vec`, `se_vec`, `eaf_vec`,
  `info_vec`):
  when the source column is present, the accessor returns a full aligned vector
  with `NaN` for missing or unmatched rows; when the source column is absent,
  the accessor returns the language-specific missing optional-field sentinel
  from [SPEC.md](SPEC.md) (`None` in Python, `[]` in MATLAB/Octave).
- `create_sumstats(...)` validates vector lengths against `reference.num_snp`.
  Unknown shapes fail clearly; required `pvec` values must satisfy the same
  finite numeric range contract as loaded objects, with `p == 0` allowed.
  Optional vectors follow the same absent/present sentinel semantics.
- cache save/load must preserve the same optional-field semantics (field absent
  remains absent; field present remains a vector).
- `Sumstats.save_cache(...)` is a thin convenience method equivalent to
  `save_sumstats_cache(sumstats, ...)`.
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native vectors.
- `Sumstats.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md).
- cache payloads store per-shard reference checksums for compatibility checks
  under the general reference-compatibility contract in
  [reference.md](reference.md).
