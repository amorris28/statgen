# statgen phased implementation plan

This plan turns the object-centered specification in `spec/` into a
working Python and MATLAB/Octave package. The guiding constraints are:

- loaders create immutable, load-only objects from documented disk formats;
- all per-SNP vectors and matrix rows align to `ReferencePanel` order;
- portable storage stays language-agnostic except for LD, whose distribution
  artifacts are runtime-native sparse matrix files by spec;
- MATLAB/Octave and Python behavior is validated from the same fixtures;
- genotype matrix access remains deferred until its contract is specified.

## Phase 0: package skeleton and shared test fixtures

Goal: establish the repository layout, dependency assumptions, and fixture data
used by every later phase.

Implementation tasks:

- create `python/pyproject.toml` and `python/statgen/` with empty module stubs
  for `reference`, `sumstats`, `annotations`, `ld`, `genotype`, and `_utils`;
  package is installable via `pip install -e python/` from the repo root;
- create `matlab/+statgen/` package with a minimal `version.m` function as the
  first callable entry point, plus the corresponding `statgen_version.m` flat
  wrapper using the `max(nargout, 1)` varargout pattern; no other stubs needed
  yet;
- create `script/` with placeholder scripts for LD conversion/build tooling;
- create `tests/conftest.py` with:
  - a session-scoped `octave_available` fixture that probes `shutil.which("octave")`;
  - a `pytest.mark.octave` marker (registered in `pyproject.toml`);
  - a `skipif_no_octave` skip decorator;
- create `tests/fixtures/generate.py`: standalone script that writes all
  committed fixture files deterministically using explicit little-endian dtypes;
  run with `python tests/fixtures/generate.py` to regenerate;
- generate and commit baseline fixtures under `tests/fixtures/` (this is the
  authoritative concrete baseline fixture composition):
  - one autosome and chrX;
  - sharded `.bim` files and one non-sharded `.bim`;
  - inputs using genomatch NCBI contig naming;
  - small BED annotation files with overlapping, adjacent, and boundary
    intervals;
  - one `.tsv.gz` sumstats file with finite, missing, `p == 0`, invalid, and
    absent-reference variants;
  - minimal LD `.npz` and converted `.mat` distribution files with manifests,
    including chrX `male`, `female`, and optional `combined` shards;
  - tiny PLINK bfile metadata fixtures for genotype metadata loading.

Tests and acceptance criteria:

- pytest runs without Octave and skips `@pytest.mark.octave` tests cleanly;
- at least one `@pytest.mark.octave` smoke test: invoke `statgen.version()` via
  subprocess, verify the returned value matches the expected string;
- fixture-generation script is deterministic and its output is committed;
- `Makefile` with three thin targets: `install`, `fixtures`, `test`;
- CI/local test commands are documented in `README.md`.

## Phase 1: shared conventions and reference panels

Goal: implement the base coordinate system that every other object depends on,
while preserving the post-genomatch variant contract exactly.

Implementation tasks:

- implement BIM parsing and validation for required columns:
  `chr`, `snp`, `cm`, `bp`, `a1`, `a2`;
- implement contig label validation per `spec/contigs-and-shards.md`:
  silently ignore `Y`/`MT` rows on load; reject `chr1`/`chrX`-style labels
  with a clear error;
- implement allele syntax validation: both `a1` and `a2` must be non-empty
  uppercase DNA strings; either allele may be multi-base; loaders preserve
  values as-is without normalization;
- loaders must never swap or reinterpret `a1`/`a2` roles and must not attempt
  reference-genome truth validation for allele identity;
- implement row-order validation per `spec/contigs-and-shards.md`: sort key is
  `(chr_rank, bp_numeric, a1_lexicographic, a2_lexicographic)`; duplicate
  `(chr, bp, a1, a2)` tuples must fail; violations must fail clearly without
  reordering or deduplication;
- implement reference checksum computation from exact loaded `chr:bp:a1:a2\n`
  rows with lowercase MD5 output;
- implement Python `ReferenceShard` and `ReferencePanel` classes with
  read-only accessors and shard offsets;
- implement MATLAB/Octave reference loading, checksum, and panel accessors;
- implement `load_reference(path, shards=None)` for:
  - sharded path templates containing `@` (canonical-label substitution, no glob);
  - single BIM split by chromosome;
  - optional `shards` subset per `spec/contigs-and-shards.md`;
- implement `ReferencePanel.select_shards(shards)` returning a subset panel;
- implement reference cache save/load in both languages; `load_reference_cache`
  accepts optional `shards` for subsetting cached shard labels;
- implement `ReferencePanel.is_object_compatible(object)` comparing ordered
  shard labels, shard row counts, and shard checksums where available; may log
  warnings but does not raise on mismatches.

Tests and acceptance criteria:

- Python and Octave return the same shard labels, row order, offsets, and
  checksums for all fixture reference layouts;
- single-file loading splits by chromosome; sharded `@` loading uses canonical
  substitution and excludes non-existent shards without error;
- `Y`/`MT` rows are silently ignored; `chr1`-style labels raise a clear error;
- sort-order violations and duplicate `(chr, bp, a1, a2)` tuples fail clearly;
- `select_shards` and cache `shards` subsetting produce correct sub-panels;
- `is_object_compatible` compares ordered shard labels, not just counts;
- compatibility checks return `false` without raising for label, row-count, or
  checksum mismatches;
- chrX is supported without chromosome 1-22 hard-coding.

## Fixture sourcing policy for Phase 2+

Follow `spec/testing.md` as the source of truth for fixture growth policy.
From Phase 2 onward, tests should use:

- Phase 2 (`Sumstats`): checked-in canonical TSV fixtures for cross-language
  alignment and `logpvec` contract; on-the-fly malformed/missing-column files
  for failure-path tests.
- Phase 3 (`Annotations`): checked-in canonical BED fixtures for overlap,
  adjacency, boundary, and exact chromosome-label behavior; on-the-fly
  permutation/stress boundary inputs.
- Phase 4a (`LD` artifact loading/validation): checked-in canonical LD `.npz`
  and `.mat` distribution fixtures with manifests for parity and compatibility;
  on-the-fly metadata/payload corruption cases.
- Phase 5a (LD artifact writing/conversion): on-the-fly synthetic
  PLINK-like tables for unit tests and on-the-fly generated `.npz`/`.mat`
  outputs for validation against Phase 4a.
- Phase 4b (`LD` panel operations): checked-in canonical fixtures plus
  Phase 5a-generated artifacts for numerical parity and selector behavior.
- Phase 5b (PLINK LD build hardening): on-the-fly generated outputs for
  integration validation; PLINK2-dependent tests are optional/skipped when
  unavailable.
- Phase 6 (`Genotype` metadata): checked-in canonical bfile metadata fixtures;
  on-the-fly inconsistent/missing shard-file cases.
- Phase 7 (workflow validation): default to checked-in canonical workflow
  inputs; only add new checked-in fixtures when a new workflow contract cannot
  be represented by temporary derived inputs.

## Phase 2: summary statistics alignment

Goal: load one portable GWAS TSV per trait/source and project it into reference
order.

Implementation tasks:

- implement required-column validation for `chr`, `bp`, `a1`, `a2`, `z`, `n`;
- enforce required numeric-field validity: `z` and `n` must parse as finite
  numbers on all source rows; invalid, `NaN`, or infinite values fail load;
- implement optional field handling for `p`, `beta`, `se`, `eaf`, and `info`;
- join rows to the reference by exact `chr:bp:a1:a2`;
- split aligned vectors into `SumstatsShard`s matching the reference shards;
- represent absent variants and missing numeric values as `NaN`;
- derive `logpvec` from `p` only:
  - missing `p` -> `NaN`;
  - `p == 0` -> `Inf`;
  - `p < 0` or `p > 1` -> `NaN`;
- implement Python and MATLAB/Octave accessors for genome-wide concatenated
  vectors;
- implement `Sumstats.select_shards(shards)` returning a subset object;
- implement `create_sumstats(reference, zvec, nvec, optional pvec, optional
  beta_vec, optional se_vec, optional eaf_vec, optional info_vec)` in both
  languages with strict length/shape validation against `reference.num_snp`;
- implement cache save/load: `load_sumstats_cache(path, shards=None)` performs
  cache-internal validation only, supports optional `shards` subsetting, and
  retains per-shard checksums for explicit post-load compatibility checks.

Tests and acceptance criteria:

- Python and Octave aligned vectors match the fixture values and missingness;
- absent reference variants remain rows with `NaN`, not dropped rows;
- optional fields have consistent missing-field sentinels;
- malformed TSV files with missing required columns fail with clear errors;
- malformed TSV files with non-finite/non-numeric required `z` or `n` fail with
  clear errors;
- factory creation validates vector lengths/shapes and preserves optional-field
  sentinel semantics;
- post-load compatibility checks via `ReferencePanel.is_object_compatible`
  detect incompatible references.

## Phase 3: annotation painting

Goal: paint BED interval annotations onto a reference panel and expose binary
SNP-level feature-mask matrices.

Implementation tasks:

- implement BED loading using only columns 1-3;
- sort and merge overlapping or adjacent intervals per chromosome;
- convert BIM 1-based `bp` to 0-based positions for interval membership;
- implement binary-search interval painting in reference row order;
- derive annotation names deterministically from BED basenames (without
  extension); fail on duplicate names;
- implement Python `AnnotationShard` and `AnnotationPanel` classes with
  `annomat` and `annonames` accessors;
- implement `AnnotationPanel.select_shards(shards)` returning a subset panel;
- implement `AnnotationPanel.select_annotations(names)` preserving requested
  order and failing on unknown names;
- implement `AnnotationPanel.union_annotations(other, mode='by_name')`
  requiring `ReferencePanel.is_object_compatible(other)` and failing on name
  collisions;
- implement `create_annotations(reference, annomat, annonames)` and
  `create_annotation(reference, annovec, annoname)` with strict
  reference-alignment and binary-value validation;
- use sparse internal representation by default and avoid dense materialization
  unless explicitly requested by caller;
- implement MATLAB/Octave annotation painting with matching semantics;
- implement cache save/load: `load_annotations_cache(path, shards=None)` does
  not require a reference object; performs cache-internal validation only and
  supports optional `shards` subsetting; stores per-shard checksums for
  post-load compatibility checking via `ReferencePanel.is_object_compatible`.

Tests and acceptance criteria:

- boundary tests cover `[start, end)` membership exactly;
- adjacent BED intervals merge and produce the same mask as separate intervals;
- chromosomes absent from a BED file produce all-zero masks;
- Python and Octave `annomat` and `annonames` match for sharded and one-shard
  reference layouts;
- Python and Octave behavior matches for sparse-vs-dense representations
  (identical logical mask semantics);
- `select_annotations` preserves requested order and fails on unknown names;
- `union_annotations` enforces reference compatibility and fails on name
  collisions;
- factory creation validates shape, binary values, and unique names;
- post-load compatibility checks via `ReferencePanel.is_object_compatible`
  detect incompatible references for annotation caches;
- complement masks and LD-weighted matrices are left to user-space arrays.

## Phase 4a: LD artifact loading and validation

Goal: load runtime-native LD distribution artifacts, validate
manifest/metadata compatibility, and construct raw LD objects without exposing
panel-level numerical operations yet.

Implementation tasks:

- implement `ld_manifest.json` parsing for panel-root loads and per-file
  metadata parsing for single-shard loads;
- implement Python `.npz` shard loading for CSC32 arrays: `data`, `indices`,
  `indptr`, `shape`, `a1freq`, and UTF-8 JSON `metadata`;
- implement MATLAB/Octave `.mat` shard loading for native sparse `ld_r`,
  `a1freq`, and metadata struct;
- implement `load_ld(path, reference, default_chrX_sex=None)` from either a
  single shard file or a panel root; for panel roots, derive required shard
  files from `reference` shard labels and resolve them through the manifest;
- construct `LDShard` objects with raw fields: `chr`, `sex`, `num_snp`,
  `ld_r`, `a1freq`, and `reference_checksum`;
- construct `LDPanel` as ordered shard groups matching the supplied reference;
  store `default_chrX_sex` minimally but defer public selector semantics;
- validate schema/version, runtime format, manifest/per-file metadata
  agreement, `num_snp`, `nnz`, shape/count fields, sparse presence, explicit
  chrX sex labels, and exact reference checksum compatibility;
- implement `validate_ld_distribution(path, check_payload_structure=False)`;
  the default path validates file checksums and manifest/per-file agreement,
  and the payload-structure path performs expensive sparse checks such as
  bounds, diagonal, and symmetry;
- for MATLAB/Octave `.mat` shards, allow v5 sparse fixtures on load and in
  validation, but have `validate_ld_distribution` warn that v5 artifacts are
  for fixtures/local tests only and are not production distributions;
- keep Python `LDShard.ld_r` as SciPy CSC and MATLAB/Octave `LDShard.ld_r` as
  native sparse double;
- expose only shard-level/raw data needed by tests and artifact verification;
  Phase 4a tests should access `LDShard.ld_r` and `LDShard.a1freq` directly
  through panel shard iteration, not through partial genome-wide panel APIs;
- do not implement LD-specific public cache APIs (`save_ld_cache` or
  `load_ld_cache`); LD distribution artifacts are already the load-efficient
  representation.

Deferred to Phase 4b:

- public `LDPanel.a1freq(chrX_sex=None)` genome-wide accessor;
- full `default_chrX_sex` and chrX selector override semantics;
- `LDPanel.select_shards(shards)`;
- `LDPanel.multiply_r2`;
- `fast_prune`.

Tests and acceptance criteria:

- Python and Octave load paired `.npz`/`.mat` fixtures into matching raw sparse
  matrices and shard-level `a1freq` vectors within numeric tolerance;
- `validate_ld_distribution` is unit-tested and is usable by Phase 5a artifact
  writer tests;
- v5 sparse `.mat` fixtures load successfully, and
  `validate_ld_distribution` emits the non-production warning for them;
- when MATLAB is available, v7.3 sparse `.mat` artifacts load and validate
  without the v5 warning;
- chrX sex label parsing is validated in the load path: `"female"`, `"male"`,
  and `"combined"` are accepted; unknown labels are rejected; manifest and
  per-file `chr`/`sex` metadata must agree;
- invalid manifests, metadata mismatches, missing required files, malformed
  sparse payloads, and reference checksum mismatches fail clearly;
- single-chromosome reference panels (one shard) load correctly.

## Phase 5a: LD artifact writing and conversion

Goal: create Python `.npz` LD distribution artifacts from synthetic/internal
inputs and convert them into MATLAB/Octave `.mat` artifacts, verifying output
with Phase 4a loading and validation.

Implementation tasks:

- implement shared artifact-writing internals for symmetric CSC32 signed-`r`
  `.npz` shards with explicit unit diagonal, aligned `a1freq`, per-file
  metadata, and `ld_manifest.json`;
- unit-test artifact writing from synthetic PLINK-like LD/frequency tables
  without requiring PLINK2;
- implement `.npz` to `.mat` conversion with separate output policy by
  runtime: MATLAB writes production v7.3 sparse `.mat`; Octave may write v5
  sparse `.mat` for fixture-scale tests and local validation only;
- converter reads `.npz` shard files written by the artifact writer, including
  only the limited `.npy` member-array subset defined by the LD `.npz` schema,
  constructs native sparse matrices once, writes `.mat` shards plus a
  MATLAB/Octave manifest, and validates metadata consistency before writing;
- converter docs must state that Octave `.mat` output is not production
  distribution output;
- generated `.npz` and converted `.mat` artifacts must pass Phase 4a
  `validate_ld_distribution` and load successfully.

Tests and acceptance criteria:

- generated `.npz` artifacts validate and load through Phase 4a;
- unit-test `.npz` to `.mat` conversion under Octave when Octave is available;
- Octave-produced v5 fixture artifacts load and validate with the
  non-production warning;
- MATLAB-produced v7.3 artifacts validate without the v5 warning when MATLAB is
  available;
- `.npz` and `.mat` converted artifacts agree on `ld_r` values within
  tolerance, `a1freq`, and metadata fields `chr`, `sex`, `num_snp`, `nnz`, and
  `reference_checksum`.

## Phase 4b: LD panel operations

Goal: expose genome-wide LD panel accessors and numerical primitives on
artifacts that Phase 4a can already load and Phase 5a can already generate.

Implementation tasks:

- implement `LDPanel.a1freq(chrX_sex=None)` with full chrX default-sex and
  override semantics;
- implement public `LDPanel.default_chrX_sex` accessor;
- implement `LDPanel.select_shards(shards)` returning a subset panel;
- implement `LDPanel.multiply_r2(M, chrX_sex=None)` without building a dense
  genome-wide LD matrix;
- implement `fast_prune(logpvec, ld_panel, r2_threshold=0.2,
  chrX_sex=None)` with stable significance ordering and per-shard pruning;
- enforce documented dtype behavior for `a1freq` and `multiply_r2`.

Tests and acceptance criteria:

- Python and Octave `a1freq`, `multiply_r2`, and `fast_prune` outputs match on
  paired `.npz`/`.mat` fixtures and Phase 5a-generated artifacts within
  numeric tolerance;
- chrX default sex follows the loader default or caller override and is
  validated against loaded chrX shards;
- `chrX_sex` override affects only chrX and is ignored when chrX is absent;
- `select_shards` follows [contigs-and-shards.md](../spec/contigs-and-shards.md);
- `multiply_r2` accepts both vectors and matrices and returns the same shape.

## Phase 5b: PLINK LD build hardening

Goal: add production-facing PLINK2 orchestration on top of the Phase 5a
artifact writer.

Implementation tasks:

- implement `script/statgen_build_ld.py`;
- support sharded bfile input with `@` and non-sharded input;
- for non-sharded input, split by chromosome;
- run PLINK2 `--freq` and `--r` with `--keep-allele-order`;
- default to a 10,000 kb LD window and `r2 >= 0.05` storage threshold, with
  command-line overrides recorded in metadata;
- for chrX, build `male` and `female` outputs by default using FAM sex codes;
  `combined` output is opt-in and records the modeling rationale in metadata;
- record build command, PLINK version if available, sample count, window,
  threshold, reference checksum, and runtime/storage format in metadata.

Tests and acceptance criteria:

- integration tests that require PLINK2 are optional and skipped when PLINK2 is
  unavailable;
- generated PLINK-backed `.npz` and converted `.mat` LD distributions pass
  `validate_ld_distribution`, load successfully through Phase 4a/4b, and agree
  on `ld_r` values within tolerance, `a1freq`, and key metadata fields;
- chrX sex splitting rejects missing FAM sex values unless `--no-sex-split` is
  used.

## Phase 6: genotype metadata objects

Goal: implement the deferred, safe subset of genotype support: metadata loading
without dense genotype access.

Implementation tasks:

- implement `GenotypeShard` and `GenotypePanel` loaders for PLINK bfile
  prefixes;
- support single-prefix input and `@`-template sharded prefixes; shard
  discovery uses canonical-label substitution (no glob), per
  `spec/contigs-and-shards.md`;
- implement `load_genotype(path_or_paths, shards=None)` with optional `shards`
  subsetting;
- implement `GenotypePanel.select_shards(shards)` returning a subset panel;
- validate `.bed`, `.bim`, and `.fam` presence;
- load BIM and FAM metadata, retaining file paths and sample order;
- for sharded genotype panels, require identical FAM rows in the same order
  across all shards, reflecting standard per-chromosome PLINK bfile practice;
- expose read-only `num_snp`, `num_sample`, `bim`, and `fam` accessors;
- do not expose genotype matrix readers until slicing, missingness, and
  reference-alignment behavior are specified.

Tests and acceptance criteria:

- Python and Octave genotype metadata accessors match fixture files;
- sharded panels reject inconsistent FAM rows or sample order across shards;
- tests do not depend on dense genotype decoding.

## Phase 7: workflow-level validation

Goal: prove that the implemented object API can support intended downstream
statistical-genetics workflows without internal-field access.

Implementation tasks:

- implement small Python and Octave examples for:
  - univariate partitioned LDSC-style array preparation;
  - bivariate/cross-trait array preparation;
  - random-prune weight estimation;
- keep regression and model fitting outside `statgen`, using plain arrays from
  accessors;
- document common usage in `README.md`.

Tests and acceptance criteria:

- examples use only public loaders, accessors, `LDPanel.multiply_r2`, and
  `fast_prune`;
- Python and Octave fixture workflows produce matching arrays within tolerance;
- no test reaches into private shard internals to reimplement loading or LD
  logic.

Workflow test sketches:

```text
ref    = load_reference(ref_path)
ld     = load_ld(ld_path, ref)
annot  = load_annotations(bed_paths, ref)
ss     = load_sumstats(sumstats_path, ref)

a1freq = ld.a1freq
sig2_i = 2 * a1freq * (1 - a1freq)
annomat_ld = ld.multiply_r2(annot.annomat * sig2_i[:, None])

mask = isfinite(ss.zvec) & isfinite(ss.nvec) & (a1freq > maf_threshold)
randvec = uniform(0, 1, ref.num_snp)
randvec[~mask] = NaN
mask = mask & isfinite(fast_prune(randvec, ld, r2_threshold=0.1))
weights = 1.0 / ld.multiply_r2(mask.astype(float))[mask]
```

```text
ss1 = load_sumstats(sumstats_path_1, ref)
ss2 = load_sumstats(sumstats_path_2, ref)
mask = isfinite(ss1.zvec) & isfinite(ss2.zvec) & isfinite(ss1.nvec) & isfinite(ss2.nvec)
zprod = ss1.zvec[mask] * ss2.zvec[mask]
```

Random-prune weight estimation should be documented and tested as repeated
calls to `fast_prune` with excluded SNPs pre-set to `NaN`:

```text
weights = zeros(ref.num_snp)
for i in 1..n_iter:
    randvec = uniform(0, 1, ref.num_snp)
    randvec[~mask] = NaN
    weights += isfinite(fast_prune(randvec, ld, r2_threshold))
weights /= n_iter
```

## Phase 8: hardening, packaging, and documentation

Goal: make the package usable by project pipelines and stable enough for
external analysis scripts.

Implementation tasks:

- document installation and optional dependencies:
  Python numeric stack, Octave, PLINK2, and any cache libraries;
- finalize public API names and import paths in Python;
- finalize MATLAB/Octave path setup and class/function naming;
- add concise API reference pages or docstrings for each object family;
- add command-level documentation for `statgen_build_ld.py`;
- add schema-version checks for portable metadata and caches;
- define backwards-compatibility policy for schema `0.1`.

Tests and acceptance criteria:

- full pytest suite passes locally;
- Octave tests pass when Octave is installed and skip cleanly otherwise;
- documentation examples run against fixture data;
- project-specific pipelines can load independent object paths without a
  bundle-root manifest.

## Dependency ordering

The critical path is:

1. shared conventions and references;
2. objects that align to references without sparse LD: sumstats and annotations;
3. LD artifact loading and validation;
4. LD artifact writing and conversion;
5. LD panel numerical operations;
6. PLINK-backed LD build hardening;
7. genotype metadata;
8. workflow examples and hardening.

Genotype dense matrix access should remain out of scope until the genotype spec
is expanded. Collection manifests should also remain optional until object
loaders and caches are stable.

## Language implementation sequencing

Python and MATLAB/Octave should implement the same public object contracts, but
they should not be developed in lockstep while low-level details are still
moving. Python should lead each implementation phase, and MATLAB/Octave should
catch up at stable phase gates using the same portable fixtures where formats
are portable, paired LD runtime distributions where LD is involved, and
language-neutral expected outputs.

Recommended sequencing:

1. Implement Phase 0 and Phase 1 in Python first, including fixture generation,
   reference checksums, reference loading, and compatibility checks.
2. Record fixture expectations in language-neutral files or pytest helpers that
   compare plain arrays and metadata, not Python object internals.
3. Implement MATLAB/Octave Phase 1 against those same fixtures before advancing
   too far into higher-level object behavior.
4. Keep Python roughly one phase ahead for sumstats and annotations, then add
   MATLAB/Octave parity before treating each phase as complete.
5. Implement LD as `4a -> 5a -> 4b -> 5b`: first artifact loading and
   `validate_ld_distribution`, then artifact writing/conversion, then numerical
   panel operations, then PLINK-backed build hardening. Python may lead within
   each LD subphase, but MATLAB/Octave parity for a subphase should close
   before advancing to the next LD subphase, except for explicitly optional
   MATLAB-only v7.3 production-conversion coverage. This keeps the write path
   testable against the real load/validation contract without blocking on
   `multiply_r2` or `fast_prune`.
6. Use cross-language tests as phase gates: a phase is not complete until
   Python and MATLAB/Octave produce matching values from shared portable
   fixtures or paired LD `.npz`/`.mat` fixtures, except for explicitly
   deferred genotype matrix access.

This is Python-first, not Python-only. MATLAB/Octave parity remains part of the
acceptance criteria, but it should follow stable contracts instead of driving
early implementation choices before the storage and cache details are proven.

## Open implementation decisions

- Python cache format for reference, sumstats, and annotations.
- MATLAB/Octave cache field layout: it should optimize simple load/save and
  checksum validation, not mirror Python internals.
- Whether the Python package should be installable as part of the repository or
  via a nested `pyproject.toml` under `python`.
- Exact implementation detail for MATLAB-compatible large `.mat` writing when
  Octave cannot produce a target panel artifact.
- How much of existing `ld/*` code should be mined for tests or conversion
  details after the clean LD implementation is underway.
