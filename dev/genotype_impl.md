# Genotype implementation plan

This plan implements the current contract in `spec/genotype.md`. The genotype
object is a reference-aligned accessor object: it caches metadata and source
row mappings, while hardcall payloads remain in PLINK `.bed` files and are read
on demand.

The existing Phase 6 text in `dev/implementation_plan.md` is older than the
current genotype spec. Treat this file and `spec/genotype.md` as the working
plan for genotype implementation.

## Goals

- Implement `load_genotype(bfile_prefix, reference)`.
- Implement `GenotypeShard` and `GenotypePanel` with immutable, read-only
  accessors.
- Align genotype source BIM rows to a required `ReferencePanel` using
  `(shard label, bp, a1_hash64, a2_hash64)`.
- Expose panel-level `num_snp` as the total reference-axis SNP count.
- Expose reference-axis `is_present`, `ploidy_male`, and `ploidy_female`.
- Expose panel-level FAM columns plus `is_male`, `is_female`, and
  `is_subject_present(shard)`.
- Implement on-demand `fetch_genotypes_int8(snp_indices, optional bed_path)` and
  `fetch_genotypes(snp_indices, optional bed_path)`.
- Implement genotype metadata cache save/load without saving genotype calls.
- Provide Python and MATLAB/Octave implementations with matching behavior.

## Non-goals

- No support for PGEN, BGEN, VCF, dosage, or non-PLINK bfile inputs.
- No allele normalization, allele swapping, strand repair, contig aliasing, or
  silent dropping of ambiguous variants.
- No dense genotype matrix cache.
- No public `bim` accessor on `GenotypePanel`.
- No sample filtering or subject remapping beyond the chrX subset expansion
  specified below.

## Shared design

### Input resolution

`bfile_prefix` is a PLINK bfile prefix without suffixes.

- Without `@`, load one non-sharded bfile triplet:
  `<bfile_prefix>.bed`, `<bfile_prefix>.bim`, `<bfile_prefix>.fam`, optional
  `<bfile_prefix>.ploidy`.
- With `@`, substitute each supplied reference shard label into the prefix and
  require every resolved shard triplet to exist.
- Missing requested source shards are errors. Users subset the reference first
  when loading a subset.
- Extra source files on disk are ignored for `@` templates because discovery is
  reference-label substitution, not globbing.

### Object structure

`GenotypePanel` holds ordered `GenotypeShard`s matching the supplied
`ReferencePanel`. It also stores `source_layout`, either `"non_sharded"` for a
single bfile source or `"sharded"` for an `@` source. `select_shards` and cache
subsetting preserve `source_layout`; do not infer source layout from equality of
stored `bed_path` values because single-shard subsets of sharded metadata still
have sharded physical layout.

`GenotypeShard` is the unit of `.bed` access. It is paired 1:1 with one
`ReferenceShard` and owns:

- `bed_path`
- `bed_file_size`
- `source_num_snp` (physical-source metadata: the total number of SNP rows in
  the `.bim`/`.bed` file serving this shard, not the number of rows matched to
  this reference shard)
- `source_num_sample`
- `source_row0`
- `subject_present`
- `source_subject_row0`
- `is_present`
- `ploidy_male`
- `ploidy_female`
- `reference_checksum`

For non-sharded bfile input, multiple `GenotypeShard`s may point to the same
source `.bed`, record the same full-file `source_num_snp`, and use different
`source_row0` vectors into that shared file. For `@`-sharded input, each shard
normally points to its resolved shard `.bed`.

### Variant alignment

For each source shard:

- parse BIM with explicit schema;
- validate contig labels and row order per `contigs-and-shards.md`;
- compute source allele hashes with the same algorithm as `ReferencePanel`;
- fail on duplicate source matching keys within a source shard;
- align source rows to the paired reference shard on
  `(bp, a1_hash64, a2_hash64)` after exact shard-label matching;
- set `is_present == true` for matched rows and `false` for unmatched
  reference rows;
- set `source_row0` to zero-based source BIM row indices and `-1` for
  unmatched rows.

For non-sharded input, parse the single BIM once, split source rows by `chr`
into canonical/reference shard groups, and keep `source_row0` as indices into
the original single source BIM/BED row order.

### Ploidy

- If `.ploidy` is absent, matched source rows have `(ploidy_male,
  ploidy_female) = (2, 2)`.
- If `.ploidy` is present, parse two integer columns with no header and require
  row count equal to `source_num_snp` (the row count of the BIM file it
  sidecars, not the number of rows matched to the reference shard).
- Allowed ploidy values are `0`, `1`, and `2`.
- Reference rows with `is_present == false` have `NaN` ploidy values.

### Sample axis

Parse FAM with explicit six-column schema. Expose:

- `fid`
- `iid`
- `father_id`
- `mother_id`
- `sex`
- `is_male`
- `is_female`

Validate:

- every FAM has unique `(fid, iid)` pairs;
- `sex` is exactly one of PLINK values `1`, `2`, or `0`;
- phenotype is parsed for row shape only and ignored for equality checks.

For non-sharded input, the single FAM defines the panel-level sample axis.

For sharded input:

- loaded autosomal shards must have identical FAM rows in the same order across
  `fid`, `iid`, `father_id`, `mother_id`, and `sex`;
- if any autosomal shard is loaded, its shared FAM axis defines the panel-level
  sample axis;
- chrX may have a FAM subset;
- chrX subjects match to the panel sample axis by unique `(fid, iid)`;
- chrX subject metadata must match the panel-level FAM exactly for matched
  subjects;
- chrX cannot contain subjects absent from the panel-level sample axis;
- if chrX is loaded without autosomes, chrX FAM defines the panel-level sample
  axis and `is_subject_present("X")` is all true.

Autosomal `subject_present` is all true. ChrX `subject_present` marks the
subset present in chrX FAM. `source_subject_row0` maps panel-level subject rows
to source FAM row indices and is `-1` for subjects absent from the shard.

### Performance rules

These rules apply to both runtimes and all phases. They restate the normative
requirements from `spec/performance-contract.md` for genotype-specific paths.

**No per-row loops over BIM, FAM, or PLOIDY rows.** Parsing, validation
(column types, allowed values, uniqueness), and source-to-reference alignment
must all use vectorized/columnar operations — pandas/NumPy in Python, native
numeric sort/merge in MATLAB/Octave. The performance-contract prohibition on
per-row loops applies in full to BIM/FAM/PLOIDY. Per-file loops (one iteration
per shard) are fine; per-row loops inside a file are not.

**The one approved SNP loop: `fetch_genotypes_int8`.** Each SNP occupies a
separate packed BED row requiring an individual disk seek; a loop over requested
SNPs is unavoidable. The inner decode must be vectorized across
`source_num_sample` subjects. Tag this loop in both runtimes:

```python
# PERF: loop over requested SNPs retained; each SNP is a separate BED row
#       requiring an individual seek. Inner decode is vectorized across subjects.
```
```matlab
% PERF: loop over requested SNPs retained; each SNP is a separate packed BED
%       row requiring an individual seek. Inner decode is vectorized across subjects.
```

**`select_shards` must not deep-copy shard payload arrays.** Build the subset
`GenotypePanel` from references to existing `GenotypeShard` objects. Only
metadata that must differ by subset (shard list, panel-wide index vectors) is
recomputed.

**MATLAB SNP-axis panel accessors must be lazy concatenations.**
`GenotypePanel` must not materialize full panel-wide SNP-axis vectors
(`is_present`, `ploidy_male`, `ploidy_female`, `source_row0`) in the
constructor. Each accessor call must build the panel-wide vector on demand by
concatenating shard payloads, without caching the result in the panel object.
Panel-level FAM metadata (`fid`, `iid`, `father_id`, `mother_id`, `sex`,
`is_male`, `is_female`) is sample-axis metadata and is stored once on
`GenotypePanel`, not concatenated from shards. Per-shard sample mappings
(`subject_present`, `source_subject_row0`) remain shard-owned; panel methods
return the selected shard's vector rather than concatenating across shards.

## Phase 1: Python metadata loading

**Status: implemented.** `python/statgen/genotype.py` now provides
`GenotypeShard`, `GenotypePanel`, `load_genotype`, `save_genotype_cache`, and
`load_genotype_cache` for metadata loading/cache paths. Shared Python helpers
live in `python/statgen/_bfile_utils.py` and
`python/statgen/_variant_match.py`.

### Implementation tasks

- Implement `python/statgen/genotype.py` with:
  - `GenotypeShard`
  - `GenotypePanel`
  - `load_genotype`
  - `save_genotype_cache`
  - `load_genotype_cache`
- Parse BIM, FAM, and PLOIDY using dataframe readers with explicit schemas.
  BIM and FAM use `pd.read_csv` with `sep=r'\s+'` because PLINK text sidecars
  may be tab- or whitespace-delimited. Always pass `header=None`, explicit
  `names`, and explicit `dtype` where needed:
  - BIM (whitespace-delimited, 6 columns):
    `names=['chr','snp','cm','bp','a1','a2']`,
    `dtype={'chr': str, 'snp': str, 'a1': str, 'a2': str}`;
    `chr` must be `str` to prevent `X` being coerced to `NaN`;
    retain `cm` for row-count validation only, then discard.
  - FAM (whitespace-delimited, 6 columns):
    `names=['fid','iid','father_id','mother_id','sex','pheno']`,
    `dtype={'fid': str, 'iid': str, 'father_id': str, 'mother_id': str}`;
    retain `pheno` for row-count validation only, then discard.
  - PLOIDY (`'\t'`-delimited, 2 columns):
    `names=['ploidy_male','ploidy_female']`, `dtype=int`.
  No per-row Python loops over BIM, FAM, or PLOIDY rows. All validation
  (dtype checks, allowed `sex` values, uniqueness of `(fid, iid)`) and
  alignment (building `source_row0`, `subject_present`, `source_subject_row0`)
  must use vectorized pandas/NumPy operations.
- Reuse shared validation helpers where available:
  - shard-label validation
  - requested-reference shard ordering
  - allele hash computation
  - shard offsets
  - cache metadata checks
- Resolve `bfile_prefix` into source file records:
  - one source record for non-sharded input;
  - one source record per reference shard for `@` input.
- Validate `.bed`, `.bim`, `.fam` existence and optional `.ploidy`.
- Validate `.bed` magic bytes and SNP-major mode.
- Validate `.bed` file size equals
  `3 + ceil(source_num_sample / 4) * source_num_snp`.
- Warn when chrX is loaded without a `.ploidy` sidecar.
- Build immutable `GenotypeShard` objects and `GenotypePanel` accessors.
- Implement `GenotypePanel.select_shards(shards)` preserving the panel-level
  sample axis, per-shard subject masks, and `source_layout`.

### Python tests

- Non-sharded bfile loads against multi-shard reference; multiple
  `GenotypeShard`s share one `.bed`.
- `@`-sharded bfile loads against the same reference.
- Missing requested source shard fails clearly.
- Reference subset plus matching `@` source subset loads successfully.
- Extra files on disk are ignored for `@` source discovery.
- Duplicate `(fid, iid)` fails.
- Invalid sex value fails.
- Autosomal FAM mismatch fails when exposed columns differ.
- Autosomal FAM phenotype mismatch does not fail.
- chrX subset FAM succeeds when `(fid, iid)` subset metadata matches.
- chrX subset FAM may be in a different order from autosomal FAM.
- chrX FAM containing subject absent from autosomes fails.
- chrX FAM metadata mismatch for a shared subject fails.
- chrX-only load defines panel sample axis from chrX.
- chrX load without `.ploidy` emits the documented warning and defaults matched
  rows to `(2, 2)`.
- `is_present`, `ploidy_male`, `ploidy_female`, `is_subject_present`,
  `is_male`, and `is_female` match expected fixture values.

## Phase 2: Python genotype fetch

**Status: implemented.** Python `GenotypePanel.fetch_genotypes_int8` and
`GenotypePanel.fetch_genotypes` now decode PLINK BED payloads on demand,
including mixed-shard requests, repeated indices, chrX sample-axis expansion,
empty SNP requests, and `bed_path` overrides.

### Implementation tasks

- Implement `GenotypePanel.fetch_genotypes_int8(snp_indices, bed_path=None)` as
  the core BED I/O path.
- Implement `GenotypePanel.fetch_genotypes(snp_indices, bed_path=None)` as a
  thin double-precision wrapper over `fetch_genotypes_int8`.
- Accept panel-global SNP indices in Python 0-based convention.
- Preserve requested SNP order, including repeated indices.
- If `bed_path` override is supplied, validate it before addressing any shard:
  - if the override contains `@`:
    - reject when `GenotypePanel.source_layout == "non_sharded"` (error: `@`
      override incompatible with non-sharded panel metadata);
    - otherwise resolve one path per addressed shard by substituting shard
      labels.
  - if the override is a flat path (no `@`):
    - accept when `GenotypePanel.source_layout == "non_sharded"`;
    - collect `source_num_snp` and `source_num_sample` across all addressed
      shards when `GenotypePanel.source_layout == "sharded"`; if any differ,
      reject (error: flat override requires equal `source_num_snp` and
      `source_num_sample` across addressed shards);
    - apply the flat path to all addressed shards.
- Split requested indices by `GenotypeShard`.
- For each shard:
  - validate all requested shard-local SNPs have `is_present == true`;
  - resolve effective `.bed` path from shard metadata or validated override;
  - validate `.bed` magic bytes, SNP-major mode, and that actual disk size
    equals stored `bed_file_size`;
  - read only requested source BED rows using `source_row0`;
  - decode vectorized across source subjects for each SNP;
  - expand to panel-level sample axis using `source_subject_row0`;
  - fill subjects absent from the source FAM with the integer missing code.
- `fetch_genotypes_int8` returns an `int8` matrix with shape
  `(num_sample, len(snp_indices))`, using `-1` for missing calls and `0`, `1`,
  `2` for observed `a1` counts.
- `fetch_genotypes` converts the int8 matrix to `float64`, mapping `-1` to
  `NaN`.

PLINK BED decoding under the statgen `a1` count contract:

```text
00 -> 2
01 -> missing
10 -> 1
11 -> 0
```

Implement the decoder in two layers:

- an internal helper that returns compact `int8` genotype codes, with ordinary
  calls encoded as `0`, `1`, and `2`, and `-1` for missing calls;
- a thin public/output wrapper that converts those integer codes to floating
  output with `NaN` for the missing sentinel.

`GenotypePanel.fetch_genotypes_int8` is the public compact hardcall accessor for
callers that want genotype hardcalls without immediate floating-point
materialization. `GenotypePanel.fetch_genotypes` remains the double-precision
convenience accessor.

Efficient floating conversion should avoid per-element scalar loops. In
MATLAB/Octave, use the mask-assignment pattern:

```matlab
geno = nan(size(geno_int8), 'single');
for code = int8([0, 1, 2])
    geno(geno_int8 == code) = single(code);
end
```

Use `double(geno)` at the public boundary when the API requires double output.
Python should use equivalent vectorized masking or table lookup.

`mostest/PlinkRead_binary2.m` is a useful source of BED mechanics, not a
definitive design. Retain these lessons:

- validate the first three BED bytes before decoding:
  - byte 1: `108`
  - byte 2: `27`
  - byte 3: `1` for SNP-major mode;
- compute `bytes_per_snp = ceil(source_num_sample / 4)`;
- for one zero-based source SNP row `source_row0`, seek to
  `3 + source_row0 * bytes_per_snp`;
- read exactly `bytes_per_snp` bytes and fail if fewer bytes are returned;
- decode packed bytes through a 256-by-4 lookup table, producing up to four
  subject calls per byte;
- trim the decoded vector to `source_num_sample`;
- keep decoding vectorized across source subjects for each SNP; a loop over
  requested SNPs is acceptable because each SNP is a separate packed BED row;
  tag the loop with the PERF comment from the Performance rules section above.

Generate the lookup table programmatically and cache it rather than maintaining
a 256-row literal table. In MATLAB/Octave this should be a `persistent` helper
so the table is generated at most once per session. In Python this can be a
module-level lazy constant. A third-party dependency may be considered only if
it preserves this exact contract and does not complicate MATLAB parity.

### Python tests

- Fetch one SNP, multiple SNPs within one shard, and SNPs across shards.
- Fetch preserves requested order.
- Repeated SNP indices produce repeated columns.
- Requesting an `is_present == false` SNP fails before reading payload rows.
- Missing hardcalls decode to `-1` in `fetch_genotypes_int8` and `NaN` in
  `fetch_genotypes`.
- All four PLINK two-bit states decode to expected values.
- chrX subset fetch expands rows to panel sample axis and fills absent subjects
  with missing calls.
- chrX subset with different FAM order still maps correctly.
- `bed_path` override works for a single non-sharded `.bed`.
- `bed_path` override with `@` works for sharded `.bed` payloads.
- `bed_path` `@` override against a non-sharded panel fails before reading any
  `.bed` (`source_layout == "non_sharded"`).
- Flat `bed_path` override against a sharded panel with unequal `source_num_snp`
  across addressed shards fails before reading any `.bed`.
- `@` override remains allowed for a single-shard subset loaded from sharded
  metadata because `source_layout` is preserved through subsetting.
- `bed_path` override size mismatch fails.
- `.bed` files with trailing bytes fail the exact-size check.
- Invalid magic bytes and non-SNP-major mode fail.

## Phase 3: Python genotype cache

**Status: implemented.** Python genotype metadata cache save/load, schema
validation, optional shard subsetting, lazy BED validation, alternate
`bed_path` fetch behavior, and bad-cache coverage are implemented.

### Implementation tasks

- Implement cache writer using NumPy `.npz` format, consistent with other
  Python caches in `statgen` (sumstats, reference).
- Store schema `"genotype_cache/0.1"`.
- Store panel-wide SNP-axis arrays:
  - `is_present`
  - `ploidy_male`
  - `ploidy_female`
  - `source_row0`
- Store per-shard sample-axis mappings:
  - `subject_present`
  - `source_subject_row0`
- Store panel-level FAM arrays:
  - `fid`
  - `iid`
  - `father_id`
  - `mother_id`
  - `sex`
  - `is_male`
  - `is_female`
- Store metadata:
  - `n_shards`
  - `shard_labels`
  - `shard_checksums`
  - `shard_start0`
  - `shard_stop0`
  - `source_layout`
  - `bed_paths`
  - `bed_file_sizes`
  - `source_num_snp`
  - `source_num_sample`
  - `num_sample`
- `load_genotype_cache(path, shards=None)` performs cache-internal validation
  and optional shard subsetting without reparsing `.bim`, `.fam`, `.ploidy`,
  or `.bed`. Cache load validates `source_layout` is one of `"non_sharded"` or
  `"sharded"` and validates
  `3 + ceil(source_num_sample / 4) * source_num_snp == bed_file_size` for
  each shard as a metadata integrity check.
- BED existence, magic bytes, and actual-disk-size-versus-stored-`bed_file_size`
  checks remain lazy in genotype fetch accessors.

### Python tests

- Cache round-trip preserves all accessors and fetch behavior.
- Cache load does not require original `.bed` to exist until fetch.
- Cache load with `shards=["X"]` preserves the original panel-level sample
  axis and chrX subject mask.
- Bad schema, inconsistent vector lengths, bad shard offsets, and unknown
  requested cache shard fail clearly.
- Cache loaded with alternate `bed_path` fetches correctly.

## Phase 4: MATLAB/Octave metadata loading

**Status: implemented.** MATLAB/Octave metadata loading is implemented with
shared BFILE utilities, `load_genotype`, `GenotypeShard`, `GenotypePanel`,
cache APIs, and the flat loader wrapper.

### Implementation tasks

- Add MATLAB/Octave files:
  - `matlab/+statgen/load_genotype.m`
  - `matlab/+statgen/GenotypeShard.m`
  - `matlab/+statgen/GenotypePanel.m`
  - `matlab/+statgen/save_genotype_cache.m`
  - `matlab/+statgen/load_genotype_cache.m`
  - flat wrapper `matlab/statgen_load_genotype.m`
- Parse BIM, FAM, and PLOIDY with `textscan` using explicit format strings:
  - BIM: reuse the `read_bim_tabular_` pattern from `load_reference.m`
    (`'%s%s%f%f%s%s'` format accepting tabs or runs of spaces as delimiters);
    `chr` must be `%s` to preserve `X`.
  - FAM (whitespace-delimited): `'%s%s%s%s%d%d'` with default whitespace
    delimiter; the first four fields (FID, IID, father, mother) are strings,
    the last two (sex, pheno) are integers. Discard `pheno` after reading.
  - PLOIDY (tab-delimited, 2 integer columns): `'%d%d'` with
    `'Delimiter', '\t'`, `'Whitespace', ''`, `'MultipleDelimsAsOne', false`.
  No per-row MATLAB loops over BIM, FAM, or PLOIDY rows. All validation and
  alignment must use vectorized native operations (numeric sort/merge matching
  the approach in reference and sumstats loaders).
- Store FID/IID and other FAM string fields as cell arrays of character
  vectors in public accessors.
- Use column vectors for per-SNP and per-sample accessors.
- Match Python validation behavior and error semantics as closely as practical.
- Use `statgen.internal.match_shard_numeric` for reference-to-genotype variant
  alignment so genotype matching shares the same numeric sort/merge helper as
  MATLAB/Octave sumstats. Compute source BIM allele hashes with
  `statgen.internal.allele_hash64`, then match each reference shard on
  `(bp, a1_hash64, a2_hash64)` within exact shard labels.
- Implement `select_shards` preserving panel-level sample axis; do not
  deep-copy shard payload arrays.
- Panel-wide accessors (`is_present`, `ploidy_male`, `ploidy_female`,
  `source_row0`, `fid`, `iid`, etc.) must be lazy concatenations per the
  Performance rules above; do not materialize them in the `GenotypePanel`
  constructor.

### MATLAB/Octave tests

- Mirror Python metadata loading fixtures.
- Verify accessors match Python expected values.
- Verify FAM sex masks and chrX subject masks.
- Verify chrX load without `.ploidy` emits the documented warning.
- Verify failure paths for missing shards, duplicate subjects, invalid sex,
  FAM mismatch, and `.ploidy` row-count mismatch.

## Phase 5: MATLAB/Octave genotype fetch

**Status: implemented.** MATLAB/Octave genotype fetch is implemented with
on-demand BED decoding, mixed-shard requests, repeated indices, chrX
sample-axis expansion, empty SNP requests, and `bed_path` overrides.

### Implementation tasks

- Implement `GenotypePanel.fetch_genotypes_int8(snp_indices, optional bed_path)`
  and `GenotypePanel.fetch_genotypes(snp_indices, optional bed_path)` in
  MATLAB/Octave convention.
- Accept MATLAB/Octave 1-based panel-global SNP indices.
- Split requests by `GenotypeShard`, decode source `.bed` rows, and reassemble
  requested column order.
- Use a generated and cached lookup-table helper equivalent to the BED
  mechanics in `mostest/PlinkRead_binary2.m`.
- Prefer an internal integer-code decoder plus a thin floating conversion
  wrapper, matching the Python implementation structure.
- Decode each SNP vectorized across source subjects.
- `fetch_genotypes_int8` returns an `int8` matrix with shape
  `num_sample x length(snp_indices)`.
- `fetch_genotypes` returns a double matrix with shape
  `num_sample x length(snp_indices)`.
- Expand chrX subset rows to panel sample axis and fill absent subjects with
  missing calls.

### MATLAB/Octave tests

- Cross-runtime fixture fetches match Python values within exact/NaN-aware
  comparison.
- Mixed-shard fetch preserves requested order.
- Repeated SNP indices produce repeated columns.
- chrX subset expansion matches Python.
- `bed_path` override with and without `@` works.
- Invalid `.bed` payload checks fail clearly.

## Phase 6: MATLAB/Octave cache

**Status: implemented.** MATLAB/Octave genotype caches are implemented as
native `.mat` metadata caches with format options, cache-only load, optional
shard subsetting, preserved panel sample axis, preserved `source_layout`, and
fetch behavior after cache load.

### Implementation tasks

- Write `.mat` genotype caches with top-level variables described in
  `spec/genotype.md`.
- Support cache `format` options consistent with other cache writers.
- Load cache without reparsing source BIM/FAM/PLOIDY/BED.
- Support optional cache `shards` subsetting.
- Preserve panel-level sample axis and `source_layout` for all shard subsets.

### MATLAB/Octave tests

- Cache round-trip preserves metadata and fetch behavior.
- Cache `shards` subsetting behaves like Python.
- Cross-runtime logical values match fixture source loads. Cache files do not
  need to be cross-runtime compatible.

## Fixture updates

**Status: partially implemented.** Deterministic genotype fixtures now include
sharded and non-sharded PLINK bfiles with `.ploidy` coverage. Additional fetch
fixtures for all PLINK two-bit states and malformed payload cases remain tied
to the fetch phases.

Extend `tests/fixtures/generate.py` with deterministic PLINK bfile fixtures:

- one non-sharded bfile spanning at least one autosome and chrX;
- one `@`-sharded bfile with matching autosomal FAM files;
- one chrX shard with a FAM subset that is not contiguous and not in autosomal
  order;
- PLINK BED rows covering hardcall values `0`, `1`, `2`, and missing under the
  statgen decoding convention;
- a reference containing variants absent from the genotype source to exercise
  `is_present == false`;
- optional `.ploidy` for chrX and at least one shard with no `.ploidy`.

Malformed/failure fixtures should generally be generated on the fly in tests:

- missing source shard;
- bad BED magic bytes;
- non-SNP-major BED mode;
- BED size mismatch;
- duplicate `(fid, iid)`;
- invalid FAM sex;
- chrX subject absent from autosomal axis;
- chrX metadata mismatch;
- malformed or wrong-length `.ploidy`.

## Definition of done

- Python and MATLAB/Octave implement all public genotype APIs in
  `spec/genotype.md`.
- Python and MATLAB/Octave agree on fixture metadata, masks, ploidy vectors,
  and selected genotype slices.
- Genotype cache round-trips preserve object behavior without storing genotype
  calls.
- `fetch_genotypes_int8` and `fetch_genotypes` support mixed-shard requests,
  repeated indices, chrX FAM subsets, and `bed_path` overrides.
- Narrow genotype tests pass in Python; Octave tests skip cleanly when Octave
  is unavailable.
