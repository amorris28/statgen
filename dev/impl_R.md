# R implementation plan

This plan covers the R runtime for `statgen`. It is implementation guidance, not
normative contract. The source of truth remains `spec/SPEC.md` and the object
specs under `spec/`.

The older `dev/implementation_plan.md` is a legacy Python/MATLAB phase plan and
should not be treated as authoritative for R work.

## Goals

- Implement the language-agnostic object contracts in R for reference,
  annotations, genotype, summary statistics, and LD panels.
- Keep the R package CRAN-compatible: no required Python, MATLAB/Octave, PLINK,
  internet access, large bundled data, or generated code during `R CMD check`.
- Use R-native caches and R-native LD distribution files while preserving shared
  logical behavior and cross-runtime fixture parity.
- Keep spec documents minimal by moving scaffold details, test inventories, and
  sequencing notes into this plan.

## Non-goals

- R does not build LD panels from genotype data.
- R does not provide a general NumPy, MATLAB, or PLINK file conversion layer.
- R caches are not cross-runtime interchange formats.
- CRAN tests do not run Python, MATLAB/Octave, PLINK, shell scripts, or network
  workflows.

## Spec boundaries

- Keep the R LD storage and conversion contract in `spec/ld.md`.
- Keep `spec/R.md` limited to package layout, naming, dependencies, R type
  conventions, error handling, and the normative CRAN-facing test subset.
- Keep detailed R package scaffold examples, internal file names, expanded test
  lists, CRAN-readiness checklist items, and release-gating workflow ideas here.
- Keep implementation suggestions for the limited `.npz` handoff reader here
  unless they are required compatibility constraints.

## Phase 0: R package skeleton

Implementation tasks:

- Create `R-package/` as the CRAN package root with `DESCRIPTION`, `NAMESPACE`,
  `R/`, `man/`, `tests/testthat/`, and `inst/extdata/`.
- A reasonable starting layout is one domain file each for reference,
  annotations, genotype, sumstats, LD, plus internal I/O, validation, and hash
  helpers. This is illustrative; implementation may split helpers differently
  as long as non-public functions are not exported.
- Private helpers should not be exported and can use a leading dot when that is
  consistent with local R style. Internal files may group helpers by domain.
- Prefer S3 objects backed by ordinary lists unless implementation experience
  shows a clearer R-native representation. There are no public constructors;
  objects are created by exported loader/cache/converter functions after
  validation.
- Add minimal package load path and version/API smoke tests.
- Add required CRAN dependencies conservatively: `Matrix`, `jsonlite`,
  `digest`, and `bit64`; add `testthat` in `Suggests`.
- Use the locally available R version for development smoke work when practical
  (for example, R 4.3.3 is acceptable), but CI and release-gating checks should
  include the latest R release because current CRAN dependency behavior and
  `R CMD check` expectations follow current R, not only older local runtimes.
- Establish fixture copying or generation policy for tiny R package fixtures
  under `R-package/inst/extdata/`: portable-format duplicates of canonical
  repository fixtures are gitignored generated copies prepared by
  `make prepare-r-fixtures`, while R-native fixtures with no portable
  counterpart may be committed directly.
- Ensure repository R entry points that run package tests or `R CMD build`
  depend on prepared R fixtures, and add pytest coverage that prepared
  duplicate fixtures match their canonical sources.
- Add repository-level pytest helpers that can call `Rscript` when R is
  available and skip cleanly otherwise.

Acceptance criteria:

- `R CMD check` can run without Python, MATLAB/Octave, PLINK, internet access,
  or large data.
- `library(statgen)` succeeds after local package installation.
- Repository pytest skips R integration tests cleanly when `Rscript` is absent.

## Phase 1: Reference panels

Implementation tasks:

- Implement BIM loading for sharded and non-sharded inputs according to
  `spec/reference.md` and `spec/contigs-and-shards.md`.
- Implement allele syntax validation, row-order validation, shard checksums, and
  `a1_hash64`/`a2_hash64` using exact `bit64::integer64` storage.
- Implement `ReferencePanel` S3 objects, read-only accessors, shard offsets,
  `select_shards`, compatibility checks, and checksum validation.
- Implement R-native reference cache save/load as RDS unless the spec later
  defines another R-native format.
- Establish the shared S3 generic pattern documented in `spec/R.md`: public
  object access goes through accessors such as `bp(x)` and `shards(x)`;
  `is_object_compatible(reference, object)` dispatches on `ReferencePanel`; and
  future reference-aligned panel classes must register `shards.<Class>` for
  compatibility checks.

Acceptance criteria:

- R reference outputs match Python/MATLAB logical fixture outputs: shard labels,
  row order, offsets, hashes, and checksums.
- Cache round-trips preserve logical reference behavior and support shard
  subsetting.

## Phase 2: Summary statistics and annotations

Implementation tasks:

- Implement summary-statistic TSV loading, exact reference matching by
  `(shard, bp, a1_hash64, a2_hash64)`, missing-row representation, optional
  fields, `logpvec`, factory construction, and cache round-trips.
- Implement BED annotation loading, interval painting, sparse annotation
  matrices, annotation naming, annotation subsetting, unions, factory
  construction, and cache round-trips.
- Use R-native tabular readers with explicit schemas on default source-load
  paths.

Acceptance criteria:

- R, Python, and MATLAB/Octave agree on aligned vectors, masks, sparse matrices,
  and selected cache round-trips for shared fixtures.
- Malformed temporary inputs fail clearly without adding broad checked-in
  fixture growth.

## Phase 3: Genotype panels

Implementation tasks:

- Implement PLINK BIM/FAM metadata loading and reference alignment.
- Implement `fetch_genotypes_int8(panel, snp_indices, bed_path = NULL)` and
  `fetch_genotypes(panel, snp_indices, bed_path = NULL, haploid_mode = NULL)`
  as the R equivalents of the logical genotype fetch APIs in
  `spec/genotype.md`.
- Returned genotype matrices use samples as rows and requested SNPs as columns:
  `num_sample x length(snp_indices)`, matching the shared API contract and
  MATLAB/Python fixture expectations.
- Implement BED hardcall decoding vectorized across subjects for each SNP.
- Implement subject metadata, presence masks, chrX subject-subset semantics,
  source layout checks, and R-native genotype metadata caches.

Acceptance criteria:

- R agrees with Python/MATLAB on genotype metadata and selected genotype slices
  through public accessors.
- Cache load skips source-style revalidation but enforces schema and reference
  compatibility gates.

## Phase 4: LD RDS loading and operations

Implementation tasks:

- Implement RDS LD shard loading from CSC payloads and construction of
  `Matrix::dgCMatrix` in memory.
- Implement `load_ld`, `load_ld_reference`, `validate_ld_distribution`,
  `a1freq`, `multiply_r2`, `fast_prune`, shard subsetting, chrX sex selection,
  and optional `retain_ld_r = false` behavior.
- Implement `multiply_r2` with `Matrix` sparse matrix multiplication over
  `Matrix::dgCMatrix` as the default path. Compiled accelerators are optional.
- Document and implement the R S3 argument order as
  `multiply_r2(ld_panel, M, chrX_sex = NULL)`, with the dispatch object first.
- Emit compatibility warnings with `warning(..., call. = FALSE)`, including
  monomorphic SNP warnings when `num_monomorphic_snps > 0`.
- Keep default LD loading on the manifest-declared reference cache path; parse
  bundled BIM files only in explicit validation paths.

Acceptance criteria:

- R LD loading and operations match Python/MATLAB logical fixture outputs within
  documented tolerances.
- R validation rejects non-R runtime manifests and malformed RDS distributions
  clearly.

## Phase 5: Python NPZ to RDS LD conversion

Implementation tasks:

- Implement `convert_ld_npz_to_rds(npz_root, rds_root, shard)` as a limited
  handoff reader for documented LD `.npz` shard payloads only.
- Implement `create_ld_rds_manifest(npz_root, rds_root, shards)`.
- Convert one requested reference shard at a time so conversions can run as
  independent parallel jobs without manifest write races. For chrX, requested
  shard `X` converts every chrX sex-label shard present in the Python manifest.
- A likely implementation path is base R `unzip()` plus `readBin()` for the
  limited `.npy` payload set (`<f4`, `<i4`, `<i8`, and `|u1`), with `bit64`
  conversion/type handling for signed 64-bit integer payloads. This is guidance,
  not a required parsing strategy; any CRAN-compatible implementation is fine if
  it preserves integer precision and satisfies the LD handoff contract.
- Reject object arrays, pickled payloads, unsupported dtypes, unsupported
  byte-order encodings, malformed shapes, and archives with missing required
  members.
- Preserve the Python shard's logical signed-`r` CSC payload while using
  R-native vector types where appropriate. Do not write a precomputed `ld_r2`
  payload; `ld_r2` is derived by the R loader.
- The manifest finalizer reads the Python `ld_manifest.json` only to determine
  which shard files are expected for the requested reference shard labels. For
  each requested label, absence from the Python manifest is an error.
- For every expected R shard file, the finalizer should check that the `.rds`
  file exists, read its metadata, verify that metadata agrees with the expected
  `(chr, sex)` and corresponding `.rds` filename, compute `file_md5`, and write
  a manifest entry from R shard metadata plus manifest-only `file` and
  `file_md5`.
- Missing expected `.rds` files or bundled `reference_bim` files are errors
  before `ld_manifest.json` is written.
- The finalizer builds a `ReferencePanel` from bundled reference BIM files for
  requested labels, saves it with `ReferencePanel.save_cache(...)` under
  `rds_root`, and records the cache filename and MD5 in the R manifest.
- The converter may assume `.npz` shards were produced and validated by
  `statgen_build_ld.py`; it validates metadata consistency before writing but
  does not need to repeat expensive O(nnz) structure checks such as full
  symmetry or diagonal scans.
- Copy bundled reference `.bim` files named by each shard's `reference_bim`
  metadata from the `.npz` panel root into the `.rds` output directory unchanged
  and carry `reference_bim` values through to R shard metadata.

Acceptance criteria:

- Converted RDS shard payloads preserve signed-`r` CSC values and metadata from
  Python handoff artifacts.
- Finalized R manifests include lowercase MD5 values, copied bundled BIM files,
  and an R-native reference cache.

## Phase 6: CRAN and repository test hardening

Implementation tasks:

- Keep CRAN-facing `testthat` tests small and self-contained under
  `R-package/tests/testthat/`.
- Keep cross-runtime parity and optional external-tool workflows in repository
  pytest or CI jobs outside CRAN checks.
- Add `docs/TUTORIAL_4_R.md` as the R-language companion to the existing Python
  and MATLAB tutorials.
- Add at least one CRAN-safe package vignette under `R-package/vignettes/`
  covering the basic load -> cache -> LD -> operation workflow using tiny
  bundled fixtures.
- Add release-gating jobs only after local R package behavior is stable.

Suggested CRAN-facing coverage:

- package load and exported API availability;
- `load_reference` from sharded and non-sharded tiny BIM inputs;
- `load_annotations` from tiny BED intervals and reference projection;
- `load_genotype` from tiny sharded and non-sharded PLINK bfiles, including
  selected genotype slices through public accessors;
- `load_sumstats` from tiny `.tsv.gz` with exact reference matching and missing
  rows represented on the aligned SNP axis;
- cache round-trips for reference, annotations, genotype metadata, and sumstats;
- loading a tiny RDS LD distribution with autosomal and chrX shards;
- LD reference loading without LD matrix loading;
- manifest filename validation and runtime-format rejection for non-R manifests
  on the R load path;
- `validate_ld_distribution(..., check_payload_structure = FALSE)` on a valid
  tiny RDS distribution;
- one `check_payload_structure = TRUE` smoke test on a very small sparse LD
  matrix;
- `a1freq`, `multiply_r2`, and `fast_prune` numerical smoke/regression tests;
- `convert_ld_npz_to_rds` and `create_ld_rds_manifest` on a tiny committed
  Python `.npz` fixture.

Suggested repository-level coverage:

- Python/MATLAB/R parity on shared portable fixtures;
- Python `.npz` to RDS conversion parity at logical shard level;
- chrX `female`, `male`, and optional `combined` selection semantics;
- malformed handoff `.npz` and malformed source-input failure paths.

Acceptance criteria:

- `R CMD check --as-cran` passes on supported local/CI platforms without
  optional external tools.
- Repository integration tests compare R against Python and MATLAB/Octave when
  those runtimes are available and skip cleanly otherwise.
- `docs/TUTORIAL_4_R.md` and at least one R package vignette are present and run
  against tiny fixture data or guard optional external-tool sections.
