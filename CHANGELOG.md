# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.3] - 2026-05-08

### Added
- Python and MATLAB/Octave can load reference-aligned PLINK 1 genotype datasets
  with `load_genotype(bfile_prefix, reference)`, including non-sharded and
  `@`-sharded bfile layouts.
- `GenotypePanel.fetch_genotypes(...)` fetches selected genotype calls without
  loading the full genotype dataset.
- Genotype metadata caches can be saved and reloaded with
  `save_genotype_cache(...)` and `load_genotype_cache(...)` to speed up repeated
  genotype loading.
- chrX genotype loading supports `.ploidy` sidecars and FAM files with chrX-only
  sample subsets.
- `GenotypePanel.is_present` and `Sumstats.is_present` provide a boolean mask in
  reference coordinates for variants present after aligning input data to the
  reference.

### Changed
- Summary statistics now require `p` values and allow `z` and `n` to be omitted

## [0.2.2] - 2026-05-06

### Added
- Sumstats loaders now accept genomatch-style summary-statistic headers in
  Python and MATLAB/Octave.
- MATLAB/Octave flat cache layouts now expose user-inspectable top-level
  variables for reference, sumstats, and annotation caches.
- Reference objects now cache deterministic `a1_hash64` and `a2_hash64` vectors
  for fast variant matching across Python and MATLAB/Octave.
- MATLAB/Octave can save thin reference caches with
  `statgen.save_reference_cache(..., 'mode', 'thin')`; thin caches store only
  metadata, base-pair positions, and allele hashes, improving
  `load_reference_cache` performance.

### Fixed
- Improved MATLAB/Octave reference and sumstats loading performance.

## [0.2.1] - 2026-05-05

### Added
- LD reference builds now follow an explicit shard-then-finalize workflow. Build
  one Python LD shard at a time, then finalize the Python manifest:

  ```bash
  for shard in 1 2 X; do
    python script/statgen_build_ld.py \
      --bfile /path/to/ref/chr@ --out /path/to/ld_npz \
      --shard "$shard" --scratch /path/to/scratch/shard_"$shard"
  done

  python script/statgen_create_ld_manifest.py --ld /path/to/ld_npz
  ```

- MATLAB/Octave LD distributions use the same shard-then-finalize pattern from
  the Python `.npz` handoff:

  ```matlab
  shards = {'1', '2', 'X'};

  for i = 1:numel(shards)
      statgen.convert_ld_npz_to_mat(ld_npz, ld_mat, shards{i});
  end
  statgen.create_ld_mat_manifest(ld_npz, ld_mat, shards);
  ```

- LD distributions now bundle the matching reference BIM shards. `load_ld(...)`
  can load the paired `ReferencePanel` directly when a reference is not supplied,
  and loaded `LDPanel` objects carry that reference for compatibility checks.
- LD builders now detect monomorphic SNPs and report them in LD metadata. By
  default, builds fail when monomorphic SNPs are present; use
  `--allow-monomorphic-snps` to force output and write a monomorphic SNP report.

### Changed
- `statgen_build_ld.py` no longer exposes PLINK sample-QC controls such as
  `--mind`. Apply sample and variant QC upstream before building the LD
  reference supplied to `statgen`.

### Fixed
- Improved MATLAB/Octave annotation loading performance by replacing the slow
  row-by-row BED payload path with vectorized parsing and painting.

## [0.2.0] - 2026-05-04

### Added
- Summary-statistics loading, alignment, creation, shard subsetting, and cache
  support for Python and MATLAB/Octave.
- BED annotation painting into sparse reference-aligned annotation panels,
  including annotation/shard selection and panel union helpers.
- LD distribution loading and validation for Python `.npz` and MATLAB/Octave
  `.mat` artifacts, including manifest and reference-checksum checks.
- LD panel operations: allele-frequency accessors, shard selection,
  element-wise `r²` multiplication, and fast LD pruning.

## [0.1.0] - 2026-04-30

### Added
- Reference loading for Python and MATLAB/Octave.
- Python and MATLAB/Octave reference cache save/load.
