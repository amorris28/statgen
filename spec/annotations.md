# AnnotationShard and AnnotationPanel

## Represents

Annotations are SNP-level feature values. A feature may represent functional
annotation, gene membership, geneset membership, other binary SNP category
definitions, or continuous numeric values aligned to SNPs.

## Disk representation

Canonical binary annotation inputs are BED interval files, one file per
annotation. Only BED columns 1–3 are used: chromosome, 0-based start, 0-based
exclusive end. Additional BED columns are ignored by the binary batch loader.
The set of binary annotations for an analysis may be one BED file per
annotation:

```text
annotations/
  baseline.bed
  coding.bed
  conserved.bed
  ...
```

Painted matrices are not canonical—they are derived from BED files and a
specific reference and must be reproducible from them.

Canonical numeric annotation inputs use the same interval-coordinate contract
with one source file per load call:

```text
chrom<TAB>start0<TAB>end0<TAB>value1<TAB>value2<TAB>...
```

The first three columns are BED coordinates: chromosome, 0-based start, and
0-based exclusive end. Columns 4 and later are finite numeric annotation values.
Headerless files are the canonical default because common baselineLD-style
continuous annotations are distributed as headerless BED-like files. Headered
files are supported for wide external tables when the loader is told that the
first row is a header and which value columns to use.

The binary 3-column BED format is a special case of interval annotation input:
inside each interval the value is implicitly `1`, and outside all intervals the
value is `0`.

Optional annotation metadata sidecars must be supplied explicitly; `statgen`
does not auto-discover sidecars from source filenames. For column-valued
numeric annotation files, a metadata sidecar is a UTF-8 text file with exactly
one physical non-empty line per physical input column, in source-column order.
Comment syntax is not special in metadata sidecars. The sidecar line count must
equal the input file column count exactly, and the sidecar must include lines
for coordinate columns 1–3 even though those lines are not emitted as output
annotation metadata. The last line may omit the trailing newline. For 3-column
binary BED inputs, metadata sidecars are file-level metadata and the entire
sidecar content is used as the annotation metadata string, preserving the text
content exactly.

Comment and blank-line handling follows each runtime's native table reader
(`pandas.read_csv`, MATLAB/Octave `textscan`, or R `data.table::fread`).

`track` and `browser` metadata lines are not supported. If present in a file,
they must be prefixed with `#` before loading.

Field separation for BED and numeric annotation data rows is a single tab
character, per the BED specification.

## In-memory objects

`AnnotationShard` and `AnnotationPanel` are reference-specific in-memory
objects. An `AnnotationShard` is produced by painting one or more interval
annotation columns onto a `ReferenceShard`; an `AnnotationPanel` is an ordered
collection of
`AnnotationShard` objects aligned to a `ReferencePanel`. The per-shard
representation keeps annotation data partitioned by chromosome, matching the
shard structure of the reference. Each shard retains the paired reference
checksum as `reference_checksum`.

After painting, each annotation occupies one column in `annomat`. The full set
of annotations is:

- `reference_checksum`: MD5 reference checksum for the paired
  `ReferenceShard`;
- `annomat`: numeric matrix, shape `num_snp x num_annot`, one column per
  annotation, rows aligned to the reference shard;
- `annonames`: string vector, length `num_annot`;
- `is_binary`: logical vector, length `num_annot`, indicating whether each
  annotation column should be treated as a binary membership annotation;
- `annotation_metadata`: string vector, length `num_annot`, containing
  user metadata for each annotation.

`annotation_metadata` entries are ordinary strings preserved by `statgen`
without parsing or validation. The recommended content is one JSON object
string per annotation. If omitted by the caller, factory functions use empty
strings and source loaders must generate stable provenance strings.
Columns with `is_binary = true` must contain only `0` and `1`. Columns with
`is_binary = false` are numeric annotation columns and may happen to contain
only `0` and `1`; loaders do not reclassify them as binary from observed
values.

## Painting Logic

Painting maps an interval file onto a `ReferenceShard` or `ReferencePanel` to
produce one or more columns in `annomat`. For 3-column BED input, the produced
column is a 0/1 membership column. For numeric input, each selected value column
produces one numeric annotation column. The binary procedure follows
`annot/paint_bed_to_bim.py`, except that `statgen` does not normalize
chromosome labels. Annotation chromosome labels and reference `chr` labels must
already use the same upstream contig naming mode, normally `genomatch` NCBI
naming.

1. **Exact chromosome matching**: compare the annotation chromosome field to the
   BIM `chr` column exactly as loaded. Do not strip `chr`, map `23` to `X`, or
   apply any other alias conversion. Inputs with inconsistent contig naming must
   be fixed upstream before annotation loading. Annotation loaders follow
   [contigs-and-shards.md](contigs-and-shards.md): canonical labels outside
   the supplied reference shard set may paint all zeros, recognized
   non-supported labels such as `Y` and `MT` are ignored, and ambiguous or
   non-canonical labels such as `chr1`/`chrX` are errors.

2. **Interval validation**: within each annotation column and chromosome, sort
   intervals by start position. For binary 3-column BED input, overlapping and
   adjacent intervals are allowed and are interpreted as a union of membership
   intervals, matching the historical BED-painting behavior. For numeric input,
   overlapping intervals are invalid and must fail with a clear error because
   they would assign multiple numeric values to the same SNP. Adjacent numeric
   intervals may be treated independently.

3. **Position conversion**: BIM `BP` values are 1-based. Convert to 0-based by
   subtracting 1 before testing interval membership.

4. **Value assignment**: for each SNP, use binary search on the sorted start
   array to find the candidate interval. A SNP at 0-based position `p` is
   inside interval `[start, end)` when `start <= p < end`. SNPs on chromosomes
   absent from the annotation file receive `0`. For 3-column BED input, SNPs
   inside an interval receive `1`. For numeric input, SNPs inside an interval
   receive the selected numeric value from that interval row.

Interval membership semantics are inherited from BED (`[start, end)`, 0-based
start, 0-based exclusive end); `statgen` does not redefine BED coordinates.

The result is a numeric vector in BIM row order, one entry per SNP.
Complement masks are a user-space operation on the returned `annomat`; the
annotation loader does not provide a negation option.

Numeric annotation values must be finite. Missing values, non-numeric values,
`NaN`, `Inf`, and `-Inf` are invalid. Categorical interval fields are out of
scope; callers should convert categorical inputs into numeric columns before
loading.

## Representation

Implementations may use language-native containers. The object is tied to one
reference panel: row order and SNP count correspond to the paired reference,
and any cache is valid only for that reference.
Internal representation should use sparse storage. `annomat` is exposed as a
sparse numeric matrix in every runtime (`scipy.sparse.csr_matrix` in Python,
MATLAB/Octave sparse matrix, and an R `Matrix::dgCMatrix`). Dense
materialization is caller-driven and explicit (for example
`toarray()`/`full(...)`).

## Cache layout

MATLAB/Octave annotation caches are `.mat` files with user-inspectable
panel-wide variables at top level:

```text
metadata
annomat
annonames
is_binary
annotation_metadata
```

`metadata` is a struct with:

```text
schema = "annotations_cache/0.2"
n_shards
shard_labels
shard_checksums
shard_start0
shard_stop0
```

`annomat` is a panel-wide sparse numeric matrix with rows aligned to the
reference panel. `annonames` is a column cell array of annotation names.
`is_binary` is a logical column vector and `annotation_metadata` is a column
cell array of metadata strings. Shard offsets are zero-based half-open
intervals into `annomat` rows and are sufficient to reconstruct
`AnnotationShard` objects. `shard_checksums` stores the per-shard reference
checksums used to restore
`AnnotationShard.reference_checksum`. Cache metadata validation should be cheap,
depending on shard count, matrix dimensions, and annotation name count rather
than scanning all SNP rows.

Implementations must support the previous `annotations_cache/0.1` binary-cache
schema and load it as a binary annotation cache by synthesizing
`is_binary = true` for every annotation and `annotation_metadata = ""` for every
annotation. The loaded `annomat` must be promoted to the current sparse numeric
runtime representation; in particular, R loaders must coerce old `lgCMatrix`
payloads to `Matrix::dgCMatrix` (for example with `as(mat, "dgCMatrix")`).
New annotation caches must be written as `annotations_cache/0.2` so older
binary-only loaders reject numeric-capable caches instead of silently coercing
continuous values.

R annotation caches are RDS files containing one named list with the same
logical top-level fields and metadata fields as the MATLAB/Octave cache layout
above. R `metadata` is a named list, `annomat` is a panel-wide
`Matrix::dgCMatrix`, `annonames` is a character vector, `is_binary` is a
logical vector, and `annotation_metadata` is a character vector.

## Panel-level accessors

`AnnotationPanel` is immutable after loading. It exposes read-only genome-wide
accessors that concatenate across shards in reference panel order. Downstream
code works with the resulting plain arrays natively.

```text
AnnotationPanel.annomat    -> num_snp × num_annot sparse numeric matrix
AnnotationPanel.annonames  -> num_annot string vector
AnnotationPanel.is_binary  -> num_annot logical vector
AnnotationPanel.annotation_metadata -> num_annot string vector
```

`annonames`, `is_binary`, and `annotation_metadata` are identical across shards
and returned once. The `annomat` type is a sparse numeric matrix.

## API

```text
load_annotations(bed_paths, reference, optional annotation_metadata,
                 optional annotation_metadata_paths) -> AnnotationPanel
load_annotation(path, reference, optional header, optional value_columns,
                optional annotation_names, optional annotation_metadata,
                optional annotation_metadata_path)
    -> AnnotationPanel
save_annotations_cache(panel, path)
load_annotations_cache(path, optional shards) -> AnnotationPanel
create_annotations(reference, annotation_matrix, annotation_names,
                   optional is_binary,
                   optional annotation_metadata) -> AnnotationPanel
create_annotation(reference, annovec, annotation_name,
                  optional is_binary,
                  optional annotation_metadata) -> AnnotationPanel

AnnotationPanel.num_snp -> int
AnnotationPanel.annomat -> num_snp × num_annot sparse numeric matrix
AnnotationPanel.annonames -> num_annot string vector
AnnotationPanel.is_binary -> num_annot logical vector
AnnotationPanel.annotation_metadata -> num_annot string vector
AnnotationPanel.select_shards(shards) -> AnnotationPanel
AnnotationPanel.select_annotations(names) -> AnnotationPanel
AnnotationPanel.union_annotations(other, optional mode) -> AnnotationPanel
AnnotationPanel.save_cache(path) -> void
```

Expected behavior:

- `load_annotations(bed_paths, reference, ...)` is the batch binary BED loader.
  `bed_paths` accepts either a single BED file path (string/path scalar) or a
  list of BED files, one per annotation. A single path is treated as a
  one-element list. Annotation names are derived deterministically from BED
  basenames (without extension). BED columns 4 and later are ignored by
  `load_annotations`; use `load_annotation` to load numeric values from those
  columns. Every annotation produced by `load_annotations` has
  `is_binary = true`.
- `load_annotation(path, reference, ...)` loads exactly one annotation source
  file. The file may produce one or more annotation columns. It supports
  3-column binary BED input and BED-like numeric interval TSV input with one or
  more selected value columns.
- `header` defaults to false. When `header = false`, `value_columns` refers to
  physical file column indices using the host language's ordinary indexing
  convention: Python uses 0-based indices, while MATLAB/Octave and R use
  1-based indices. Selected value columns must identify physical columns 4 or
  later. Named `value_columns` are invalid without a header. When
  `header = true`, `value_columns` may contain source column names or
  host-language-indexed physical file column indices. Header names are the
  recommended cross-runtime selector for headered files.
- For `load_annotation`, 3-column input with `value_columns` omitted produces
  one binary annotation named from the file stem unless `annotation_names` is
  supplied. Four-column input with `value_columns` omitted uses column 4 and
  names the output from the file stem when `header = false`, or from the source
  column name when `header = true`. Input with five or more columns requires
  explicit `value_columns`; if `header = false`, it also requires
  `annotation_names`. Three-column input produces `is_binary = true`.
  Annotation columns loaded from columns 4 and later produce
  `is_binary = false` by default, even if the observed values happen to be only
  `0` and `1`.
- If `annotation_names` is supplied to `load_annotation`, it overrides output
  names inferred from source column headers or file stem. Its length must equal
  the number of output annotation columns.
- `load_annotation` accepts at most one of `annotation_metadata` and
  `annotation_metadata_path`. If `annotation_metadata` is supplied, it must be a
  string vector with one entry per output annotation; scalar string metadata is
  not accepted for `load_annotation`, even when the source produces one output
  annotation. Pass a length-1 string vector for one-column loads. If
  `annotation_metadata_path` is supplied for column-valued numeric input, it is
  a column-aligned sidecar as defined above; selected output metadata is taken
  from the sidecar lines matching the selected physical file columns after
  converting any host-language-indexed integer selectors to physical column
  numbers. If `annotation_metadata_path` is supplied for 3-column binary input,
  the full sidecar content is used as the single output metadata string. If
  neither argument is supplied, the loader must generate stable provenance
  strings, recommended as compact JSON object strings containing at
  least `source_file`, `source_column0`, and `source_column_name`.
- `load_annotations` accepts at most one of `annotation_metadata` and
  `annotation_metadata_paths`. If `annotation_metadata` is supplied, it must be
  a string vector with one entry per BED path. If `annotation_metadata_paths` is
  supplied, it must have one entry per BED path. Each entry is either a sidecar
  path or the runtime's natural missing/empty path sentinel (`None` in Python,
  `[]` or empty string in MATLAB/Octave, and `NA_character_` in R). A sidecar
  entry is file-level metadata and the full sidecar content is used as the
  metadata string for the corresponding binary BED annotation. A missing/empty
  entry requests generated provenance for that BED file. If neither argument is
  supplied, the loader must generate stable provenance strings, recommended as
  compact JSON object strings containing at least `source_file`,
  `source_column0`, and `source_column_name`.
- Empty annotation source files are invalid input and must fail with a clear
  error. A file that becomes empty after skipping BED comment lines is also
  invalid and must fail with a clear error.
- `reference` is required; the output shard structure matches the reference
  panel; `annomat` rows are aligned to the reference.
- `annonames` must be unique. Duplicate names from BED basenames, selected
  header names, explicit `annotation_names`, or programmatic inputs are an
  error.
- caching saves and restores the painted `annomat`; the cache is a single file
  (non-sharded). `load_annotations_cache` performs cache-internal validation
  only and supports optional `shards` subsetting.
- the MATLAB/Octave cache layout is specified in the "Cache layout" section
  above; per-shard checksums are trusted from cache metadata (the painted matrix
  does not contain the raw SNP data needed to recompute them).
- Accessors are read-only, concatenate shards in reference panel order, and
  return plain language-native matrices or vectors.
- `annonames`, `is_binary`, and `annotation_metadata` must be identical across
  shards and are returned once.
- `AnnotationPanel.select_shards` shard subsetting follows
  [contigs-and-shards.md](contigs-and-shards.md).
- `create_annotations` and `create_annotation` validate strict alignment to the
  supplied `reference` (`num_snp` shape match, finite numeric values only) and
  return immutable `AnnotationPanel` objects. If `is_binary` is omitted, it is
  inferred per annotation column from exact `0`/`1` matrix values for
  convenience. If supplied, it must be a logical vector of length `num_annot`
  for `create_annotations` or a scalar logical for `create_annotation`.
  Declared `is_binary = true` columns must contain only `0` and `1`; declared
  `is_binary = false` columns are not required to contain non-binary values.
- `create_annotations` accepts one or more annotation columns. If
  `annotation_metadata` is supplied, it must be a string vector with length
  `num_annot`; scalar recycling is not allowed. If omitted, metadata defaults to
  empty strings.
- `create_annotation` is the one-column convenience wrapper. If supplied,
  `annotation_metadata` is a single scalar string; if omitted, it defaults to an
  empty string. This scalar shortcut is specific to `create_annotation`.
- `AnnotationPanel.select_annotations(names)` preserves requested name order and
  fails on unknown names or duplicate requested names (no silent drops).
- `AnnotationPanel.union_annotations(other, optional mode)` requires strict
  reference-alignment compatibility between the two annotation panels: same
  shard count, same shard labels in the same order, same per-shard row counts,
  and present, equal per-shard reference checksums. Missing or mismatched
  reference checksum metadata is an error. `mode` defaults to `by_name`; name
  collisions are errors. It column-binds `annomat` and concatenates
  `annonames`, `is_binary`, and `annotation_metadata` without coercing numeric
  values to binary.
- `AnnotationPanel.save_cache(...)` is a thin convenience method equivalent to
  `save_annotations_cache(panel, ...)`.
- cache payloads store per-shard reference checksums so compatibility with a
  `ReferencePanel` can be checked after load via
  `ReferencePanel.is_object_compatible`.
