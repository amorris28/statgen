# Contigs, Shards, and Panels

## Canonical contig set

`statgen` supports contigs `1`–`22` and `X`, in that order. This is the
canonical contig order used for all shard sequences, row-order keys, and
shard-label validation.

- `Y` and `MT` are recognized non-supported contigs. They are ignored when
  outside a loader's supported output/projection rather than normalized into a
  supported shard. Source loaders drop `Y`/`MT` rows before row-order and
  duplicate-key validation for supported output shards.
- Labels such as `chr1`/`chrX` indicate non-preprocessed input and are rejected.
- Contig labels are never normalized: `chr1`, `1`, and `NC_000001.11` are
  distinct labels; mismatches are not resolved silently.

`X` is a first-class chromosome. Implementations must not hard-code
chromosomes 1–22 in APIs, manifests, or test fixtures.

## Row-order contract

Source-like SNP tables must already be sorted by:

```text
(chr_rank, bp_numeric, a1_lexicographic, a2_lexicographic)
```

where `chr_rank` follows canonical contig order. Duplicate `(chr, bp, a1, a2)`
tuples are not allowed.

This contract applies to both non-sharded files (enforced across the full file)
and sharded files (shard sequence must follow canonical contig order; rows
within each shard must follow row order).

For performance, reference loaders validate only the ordering fields relied on
by in-memory matching: `chr_rank` and `bp_numeric`. They do not inspect
`a1_lexicographic`/`a2_lexicographic` ordering because sumstats/reference
matching uses `(shard label, bp, a1_hash64, a2_hash64)` rather than source row
order among alleles at the same base-pair coordinate. The full tuple sort order
remains a normative input contract. Loaders must not reorder or deduplicate
rows.

## Shard discovery

Sharded paths use `@` as the shard-label placeholder.

- **Reference loaders**: substitute `@` with canonical contig labels in order
  (no glob-based discovery) and load existing matches.
- **Genotype loaders**: substitute `@` with reference shard labels in order
  when a reference is supplied. Every requested reference shard must have a
  matching bfile shard; users who want a shard subset should first subset the
  reference and load genotype against that smaller reference.
- **LD loaders**: always resolve shard files through `ld_manifest.json`.
  `load_ld` and `load_ld_reference` load the manifest-declared reference cache.
  The optional `shards` parameter further subsets which shards are loaded.
  `load_ld` errors if a requested shard is absent from the reference cache or
  LD panel; `load_ld_reference` errors if a requested shard is absent from the
  reference cache.

Non-sharded single-file inputs are split by the `chr` column into
per-chromosome shards in canonical order. For source loaders with an explicit
`reference` argument, the split is a reference-driven projection: validated
source rows outside the supplied reference shard set are ignored for alignment.
Invalid or ambiguous contig labels remain errors.

## Shard subsetting

All panel objects expose `select_shards(shards)`.
Cache loaders that operate without a required reference
(`load_reference_cache`, `load_annotations_cache`, `load_sumstats_cache`,
`load_genotype_cache`) accept an optional `shards` parameter.
Source loaders with an explicit `reference` argument
(`load_annotations`, `load_sumstats`, `load_genotype`) use the supplied
reference shard structure; subsetting is done via `select_shards` on the
reference before passing it in. Such loaders do not accept an additional
`shards` parameter.
`load_ld` follows the cache-loader pattern because the LD distribution includes
its own manifest-declared reference cache. The `shards` argument is honored as
a subset of both that reference cache and the LD panel. `load_ld_reference`
honors the same `shards` argument for the manifest-declared reference cache.
The rules are uniform:

- Omitting `shards` loads or returns all available shards present in the input
  or cache.
- Provided `shards` must be a non-empty list of unique contig labels in
  canonical subsequence order.
- Requesting a shard not present in the input or cache is an error.
