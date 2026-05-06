# AnnotationShard and AnnotationPanel

## Represents

Annotations are SNP-level feature masks. A feature may represent functional
annotation, gene membership, geneset membership, or other binary SNP category
definitions.

## Disk representation

Canonical annotation inputs are BED interval files, one file per annotation.
Only BED columns 1–3 are used: chromosome, 0-based start, 0-based exclusive end.
Additional BED columns are ignored. The set of annotations for an analysis is
one BED file per annotation:

```text
annotations/
  baseline.bed
  coding.bed
  conserved.bed
  ...
```

Painted matrices are not canonical—they are derived from BED files and a
specific reference and must be reproducible from them.

BED comment lines are skipped before parsing columns 1–3:

- empty lines;
- lines starting with `#`.

Comment and blank lines are only valid before the first data row. A `#`-prefixed
or blank line after a data row has begun is not supported and will be treated as
a parse error.

`track` and `browser` metadata lines are not supported. If present in a file,
they must be prefixed with `#` before loading.

Field separation for BED data rows is a single tab character, per the BED
specification.

## In-memory objects

`AnnotationShard` and `AnnotationPanel` are reference-specific in-memory
objects. An `AnnotationShard` is produced by painting one or more BED interval
files onto a `ReferenceShard`; an `AnnotationPanel` is an ordered collection of
`AnnotationShard` objects aligned to a `ReferencePanel`. The per-shard
representation keeps annotation data partitioned by chromosome, matching the
shard structure of the reference. Each shard retains the paired reference
checksum.

Annotations loaded from BED files are binary. Non-binary annotation values are
out of scope for this phase. LD-weighted or otherwise numeric annotation
matrices are user-space arrays derived from `annomat`, not
`AnnotationPanel` fields.

After painting, each annotation occupies one column in `annomat`. The full set
of annotations is:

- `annomat`: binary matrix, shape `num_snp x num_annot`, one column per
  annotation, rows aligned to the reference shard;
- `annonames`: string vector, length `num_annot`.

## Painting logic

Painting maps a BED interval file onto a `ReferenceShard` or `ReferencePanel`
to produce a 0/1 column in `annomat`. The procedure follows
`annot/paint_bed_to_bim.py`, except that `statgen` does not normalize
chromosome labels. BED chromosome labels and reference `chr` labels must already
use the same upstream contig naming mode, normally `genomatch` NCBI naming.

1. **Exact chromosome matching**: compare the BED chromosome field to the BIM
   `chr` column exactly as loaded. Do not strip `chr`, map `23` to `X`, or
   apply any other alias conversion. Inputs with inconsistent contig naming must
   be fixed upstream before annotation loading.

2. **Interval merging**: within each chromosome, sort intervals by start
   position and merge overlapping or adjacent intervals into a non-overlapping
   sorted array of `(start, end)` pairs.

3. **Position conversion**: BIM `BP` values are 1-based. Convert to 0-based by
   subtracting 1 before testing interval membership.

4. **Mask assignment**: for each SNP, use binary search on the merged start
   array to find the candidate interval. A SNP at 0-based position `p` is
   inside interval `[start, end)` when `start <= p < end`. SNPs on chromosomes
   absent from the BED file receive `0`.

Interval membership semantics are inherited from BED (`[start, end)`, 0-based
start, 0-based exclusive end); `statgen` does not redefine BED coordinates.

The result is a boolean or 0/1 vector in BIM row order, one entry per SNP.
Complement masks are a user-space operation on the returned `annomat`; the
annotation loader does not provide a negation option.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: row order and SNP count correspond to the paired reference,
and any cache is valid only for that reference.
Internal representation should use sparse storage. `annomat` is exposed as a
sparse binary matrix in both runtimes (`scipy.sparse.csr_matrix` in Python;
MATLAB/Octave sparse matrix). Dense materialization is caller-driven and
explicit (for example `toarray()`/`full(...)`).

## Cache layout

MATLAB/Octave annotation caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
annomat
annonames
```

`metadata` is a struct with:

```text
schema = "annotations_cache/0.1"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
```

`annomat` is a panel-wide sparse binary matrix with rows aligned to the
reference panel. `annonames` is a column cell array of annotation names. Shard
offsets are zero-based half-open intervals into `annomat` rows and are
sufficient to reconstruct `AnnotationShard` objects. Cache metadata validation
should be cheap, depending on shard count, matrix dimensions, and annotation
name count rather than scanning all SNP rows.

## Panel-level accessors

`AnnotationPanel` is immutable after loading. It exposes read-only genome-wide
accessors that concatenate across shards in reference panel order. Downstream
code works with the resulting plain arrays natively.

```text
AnnotationPanel.annomat    -> num_snp × num_annot sparse binary matrix
AnnotationPanel.annonames  -> num_annot string vector
```

`annonames` is identical across shards and returned once. The `annomat` type
is a sparse binary matrix.

## API

```text
load_annotations(bed_paths, reference) -> AnnotationPanel
save_annotations_cache(panel, path, optional format)
load_annotations_cache(path, optional shards) -> AnnotationPanel
create_annotations(reference, annomat, annonames) -> AnnotationPanel
create_annotation(reference, annovec, annoname) -> AnnotationPanel

AnnotationPanel.annomat -> num_snp × num_annot sparse binary matrix
AnnotationPanel.annonames -> num_annot string vector
AnnotationPanel.select_shards(shards) -> AnnotationPanel
AnnotationPanel.select_annotations(names) -> AnnotationPanel
AnnotationPanel.union_annotations(other, optional mode) -> AnnotationPanel
```

Expected behavior:

- `bed_paths` accepts either a single BED file path (string/path scalar) or a
  list of BED files, one per annotation. A single path is treated as a
  one-element list. Annotation names are derived deterministically from BED
  basenames (without extension).
- Empty BED files are invalid input and must fail with a clear error.
  A file that becomes empty after skipping BED comment lines is also invalid
  and must fail with a clear error.
- `reference` is required; the output shard structure matches the reference
  panel; `annomat` rows are aligned to the reference.
- `annonames` must be unique. Duplicate names from BED basenames or
  programmatic inputs are an error.
- caching saves and restores the painted `annomat`; the cache is a single file
  (non-sharded). `load_annotations_cache` performs cache-internal validation
  only and supports optional `shards` subsetting.
- the MATLAB/Octave cache layout is specified in the "Cache layout" section
  above; per-shard checksums are trusted from cache metadata (the painted matrix
  does not contain the raw SNP data needed to recompute them).
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native matrices or vectors.
- `annonames` must be identical across shards and is returned once.
- `AnnotationPanel.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md).
- `create_annotations` and `create_annotation` validate strict alignment to the
  supplied `reference` (`num_snp` shape match, binary values only) and return
  immutable `AnnotationPanel` objects.
- `AnnotationPanel.select_annotations(names)` preserves requested name order and
  fails on unknown names (no silent drops).
- `AnnotationPanel.union_annotations(other, optional mode)` requires
  `ReferencePanel.is_object_compatible(other) == true`. `mode` defaults to
  `by_name`; name collisions are errors.
- cache payloads store per-shard reference checksums so compatibility with a
  `ReferencePanel` can be checked after load via
  `ReferencePanel.is_object_compatible`.
