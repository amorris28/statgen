# MATLAB/Octave implementation conventions

## Compatibility target

The MATLAB implementation is the primary runtime contract. Octave is used as a
unit-test harness and compatibility proxy, not as the normative behavior
source. When MATLAB and Octave differ, implementation behavior is defined by
MATLAB semantics unless explicitly documented otherwise.

## Package layout

Internal implementation uses the `+statgen` package directory, which Octave
resolves as the `statgen` namespace. Flat `statgen_*` wrapper scripts are added
only for the public loader entry points and top-level primitives — not for
internal classes or helpers.

```
matlab/
  +statgen/
    load_reference.m
    load_ld.m
    load_ld_reference.m
    load_sumstats.m
    load_annotations.m
    load_genotype.m
    fast_prune.m
    ReferenceShard.m
    ReferencePanel.m
    LDShard.m
    LDPanel.m
    AnnotationShard.m
    AnnotationPanel.m
    SumstatsShard.m
    Sumstats.m
    GenotypeShard.m
    GenotypePanel.m
    +internal/          % shared implementation helpers, not public API
  statgen_load_reference.m
  statgen_load_ld.m
  statgen_load_ld_reference.m
  statgen_load_sumstats.m
  statgen_load_annotations.m
  statgen_load_genotype.m
  statgen_fast_prune.m
```

The `matlab/` directory itself must be on the Octave/MATLAB path. The
`+statgen` subdirectory is discovered automatically from there; it must not
be added to the path directly.

## Wrapper pattern

Each flat entry point is a thin forwarding wrapper using the standard
`varargout`/`varargin` pattern:

```matlab
function varargout = statgen_load_reference(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_reference(varargin{:});
end
```

The `max(nargout, 1)` allocation avoids an empty index expression when the
caller discards all outputs. The internal `statgen.*` function is responsible
for all argument validation and logic; wrappers contain no other code.

Wrappers exist only for: `statgen_load_reference`, `statgen_load_ld`,
`statgen_load_ld_reference`, `statgen_load_sumstats`,
`statgen_load_annotations`, `statgen_load_genotype`, `statgen_fast_prune`.
No wrappers are created for classes or internal helpers.

Interactive `help` text conventions are specified in
[docstrings.md](docstrings.md).

## Naming conventions

- Functions and variables: `snake_case`.
- Classes: `PascalCase` matching the spec object names.
- Shared implementation helpers that are used by multiple package entry
  points live under the nested `+statgen/+internal/` namespace and are called
  with fully qualified names such as `statgen.internal.read_ld_manifest(...)`.
  The `internal` namespace is not public API and may change without
  compatibility guarantees.
- Truly file-local helpers should be subfunctions in the same `.m` file.
- `+statgen/private/` may be used only for compatibility wrappers or helpers
  whose visibility is verified under both MATLAB and Octave; shared helpers
  should prefer `+statgen/+internal/` because it uses documented namespace
  resolution and is covered by Octave CI.

## Types and return values

- Per-SNP vectors are column vectors (shape `n × 1`).
- Matrices have SNPs along rows (shape `n × k`), except genotype fetch matrices:
  `GenotypePanel.fetch_genotypes_int8` and `GenotypePanel.fetch_genotypes`
  deliberately return samples as rows and requested SNPs as columns
  (`num_sample × length(snp_indices)`) to match the natural hardcall matrix
  layout and the Python API.
- Missing numeric values use `NaN`.
- Optional arguments use `[]` as the absent sentinel; implementations test
  with `isempty(arg)`.
- String arguments are character arrays or string scalars; convert to
  `char` at the loader boundary for compatibility with Octave string handling.

Cache format requirements are defined normatively in
[performance-contract.md](performance-contract.md). In particular, MATLAB/Octave
caches must be `.mat` files with SNP-axis vectors/matrices stored as native
numeric/logical/sparse arrays.

String-typed SNP-axis cache fields may use object-specific `.mat` encodings
when native MATLAB/Octave cell-string arrays are too expensive to load. The
standard pattern for reference string fields is one per-shard character vector
containing newline-delimited values with no trailing newline, as defined in
[reference.md](reference.md). Cache-loaded objects should keep such payloads
encoded and decode them lazily only when the corresponding string accessor is
used.

## Shard representation

Panels hold an ordered cell array or struct array of shard objects. Genome-wide
accessors return a single matrix or vector in panel order.

For panel-like objects (`ReferencePanel`, `Sumstats`, `AnnotationPanel`),
genome-wide accessors should be implemented as **lazy, non-cached**
concatenations from shard payloads:

- constructors should not eagerly materialize/store full panel-wide SNP-axis
  vectors or matrices in addition to per-shard payloads;
- accessors should build the concatenated output on demand from `shards`;
- "non-cached" means the panel object should not persist a memoized full
  panel-wide SNP-axis payload (for example a hidden cached `annomat`/`zvec`)
  across accessor calls; each accessor call may recompute concatenation;
- callers that need repeated use should store one local copy, e.g.
  `A = panel.annomat`.

Panel and shard classes may define `display(obj)`/`disp(obj)` for interactive
summaries. Display output should be metadata-only and must not materialize
large genome-wide dependent accessors or lazily decoded cache string payloads.

This keeps MATLAB/Octave behavior aligned with the memory model in the object
specs and with Python panel-accessor semantics.

## Error handling

Use `error('statgen:errorID', 'message ...')` with a namespaced error ID for
all user-facing validation failures. Do not use `assert` for input validation.

Compatibility checks log via `warning('statgen:compat', ...)` for mismatches.
They return a scalar logical and do not throw.

## Runtime compatibility guardrails

- For BIM/TSV-like source inputs, native table readers must use explicit text
  schema for mixed-label columns (for example `chr`) rather than auto-inference.
  This avoids label coercion (for example `X` becoming `NaN`) and preserves
  input-contract semantics.
- Tabular file reading uses `textscan` throughout. `readtable` is not used;
  `textscan` runs identically under MATLAB and Octave and keeps a single code
  path on both runtimes. Two patterns cover all current formats:
  - **Fixed schema, no filtering** (e.g. BIM): `textscan` with a fixed
    6-column `'%s%s%f%f%s%s'` format accepting tabs or runs of spaces as
    delimiters.
    See `load_reference.m:read_bim_tabular_`.
  - **Dynamic schema or leading comment/blank filtering** (e.g. BED, sumstats
    TSV): probe-and-seek or header-read `textscan`. For BED, a `fgetl` loop
    skips leading blank and `'#'`-prefixed lines to find the first data line,
    then `fseek` repositions there; for TSV, a single `fgetl` reads the header
    line and the file position advances naturally. In both cases the format
    string is built dynamically as `repmat('%s', 1, n_cols)` from the probed
    column count. `textscan` requires an exact column count with tab delimiter;
    a fixed under-count confuses the parser. For BED, `#`-prefixed lines are
    not permitted after the first data row so `CommentStyle` is not needed.
    See `load_annotations.m:read_bed_tabular_`, `load_sumstats.m`.
- Internal struct metadata keys must be valid MATLAB identifiers (for example
  `statgen_var_names__`), not names that rely on permissive dynamic-field
  behavior (for example leading-underscore keys such as `_var_names`), because
  this can fail at runtime under MATLAB even when other environments appear to
  tolerate it.
- Use `fprintf` for formatted output in shared runtime/test snippets intended
  to execute under MATLAB. Octave-only output helpers (for example `printf`)
  must be wrapped or normalized at the harness boundary.
- Cross-runtime differences must be isolated in boundary helper functions
  (I/O parsing, hashing, output formatting, engine/session integration), not
  duplicated across core loader/object logic.
- Test execution may use one persistent MATLAB engine session per pytest run.
  MATLAB-facing code and test harness paths must not assume process-per-call
  isolation.
