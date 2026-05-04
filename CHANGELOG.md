# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
