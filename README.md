# statgen

`statgen` provides Python and MATLAB/Octave tools for working with
reference-aligned statistical genetics data: reference variants, genotypes, LD,
annotations, and GWAS summary statistics.

All panel objects use reference coordinates and are sharded internally by
chromosome. User-facing APIs work at the panel level and preserve input contig
labels and allele fields exactly. `statgen` expects post-harmonization inputs;
genome build, contig naming, and allele orientation should be resolved upstream
by [genomatch](https://github.com/precimed/genomatch) or an equivalent pipeline.

The format and behavior contracts are specified in [spec/SPEC.md](spec/SPEC.md).

## Tutorials

- [Prepare LD and object caches](docs/TUTORIAL_1_PREPARE_DATA.md): build an LD
  distribution and save Python and MATLAB/Octave caches from source files.
- [Python analysis from caches](docs/TUTORIAL_2_PYTHON.md): load caches, run LD
  operations, prune variants, and fetch genotypes.
- [MATLAB/Octave analysis from caches](docs/TUTORIAL_3_MATLAB.md): run the same
  cached analysis workflow in MATLAB or Octave.

## Object Overview

| Object | Meaning | Typical input | Cache format |
| --- | --- | --- | --- |
| `ReferencePanel` | Ordered reference variant table | PLINK `.bim`, either one file or an `@`-sharded path | Python `.npz`; MATLAB/Octave `.mat` |
| `LDPanel` | Sparse LD distribution tied to a matching reference | Built LD distribution directory | No panel cache; build or convert the distribution |
| `AnnotationPanel` | Binary annotation matrix painted onto the reference | BED files | Python `.npz`; MATLAB/Octave `.mat` |
| `GenotypePanel` | Reference-aligned PLINK genotype metadata with on-demand hardcall access | PLINK 1 bfile prefix, either one bfile or an `@`-sharded prefix | Metadata cache: Python `.npz`; MATLAB/Octave `.mat` |
| `Sumstats` | One aligned summary-statistics trait or source | `.tsv.gz` with `chr`, `bp`, `a1`, `a2`, and `p` columns | Python `.npz`; MATLAB/Octave `.mat` |

## Key APIs

### ReferencePanel

- `load_reference(path)` loads a `.bim` file or an `@`-sharded `.bim` pattern.
- `panel.chr`, `panel.snp`, `panel.bp`, `panel.a1`, and `panel.a2` expose
  reference-coordinate variant fields.
- `panel.select_shards(shards)` subsets by shard label.
- `panel.is_object_compatible(other)` checks reference compatibility.
- `panel.save_cache(path, mode="full")` and `load_reference_cache(...)` save and
  reload reference caches. Thin caches use `mode="thin"`.

### LDPanel

- `load_ld(path, reference=None, shards=None, default_chrX_sex=None)` loads a
  sparse LD distribution.
- `panel.a1freq(chrX_sex=None)` returns reference-aligned allele frequencies.
- `panel.multiply_r2(M, chrX_sex=None)` multiplies by LD `r²`.
- `panel.select_shards(shards)` subsets by shard label.
- `fast_prune(logpvec, ld_panel, r2_threshold=0.2, chrX_sex=None)` performs LD
  pruning using aligned scores.

### Building LD Distributions

- `python script/statgen_build_ld.py ... --shard SHARD` builds one LD shard.
- `python script/statgen_create_ld_manifest.py --ld PATH` finalizes a Python LD
  distribution.
- `statgen.convert_ld_npz_to_mat(input_root, output_root, shard)` converts a
  Python LD shard for MATLAB/Octave.
- `statgen.create_ld_mat_manifest(input_root, output_root, shards)` finalizes a
  MATLAB/Octave LD distribution after conversion.

### AnnotationPanel

- `load_annotations(bed_paths, reference)` paints BED intervals onto a
  `ReferencePanel`.
- `create_annotations(reference, annomat, annonames)` creates an annotation
  panel from an already aligned matrix.
- `panel.annomat` exposes the reference-aligned annotation matrix.
- `panel.select_shards(shards)` and `panel.select_annotations(names)` subset an
  annotation panel.
- `panel.union_annotations(other, mode="by_name")` combines annotation panels.
- `panel.save_cache(path)` and `load_annotations_cache(...)` save and reload
  painted annotation caches.

### GenotypePanel

- `load_genotype(bfile_prefix, reference)` loads PLINK 1 genotype metadata and
  aligns it to a `ReferencePanel`. `bfile_prefix` may be a single prefix or an
  `@`-sharded prefix.
- `panel.fetch_genotypes(snp_indices)` returns selected genotype dosages without
  loading the full BED file.
- `panel.fetch_genotypes_int8(snp_indices)` returns selected hardcalls as compact
  integer calls.
- `panel.is_present` is a boolean mask in reference coordinates for variants
  present after alignment.
- `panel.fid`, `panel.iid`, `panel.sex`, `panel.is_male`, and
  `panel.is_female` expose sample metadata.
- `panel.save_cache(path)` and `load_genotype_cache(...)` save and reload
  genotype metadata caches for faster repeated loading.

### Sumstats

- `load_sumstats(path, reference)` loads a `.tsv.gz` summary-statistics file and
  aligns rows by `chr:bp:a1:a2`.
- Input files require `p`; `z` and `n` are optional.
- `create_sumstats(reference, pvec, zvec=None, nvec=None, ...)` creates a
  `Sumstats` object from already aligned vectors.
- `sumstats.logpvec`, `sumstats.zvec`, and `sumstats.nvec` expose aligned
  vectors. Optional vectors are `None` in Python or empty in MATLAB/Octave when
  absent.
- `sumstats.is_present` is a boolean mask in reference coordinates for variants
  present after alignment.
- `sumstats.select_shards(shards)` subsets by shard label.
- `sumstats.save_cache(path)` and `load_sumstats_cache(...)` save and reload
  aligned summary-statistics caches.

## Setup

### Python

Install from a cloned repository:

```sh
pip install -e python/
```

Then import APIs from the relevant submodule, for example:

```python
from statgen.reference import load_reference
from statgen.sumstats import load_sumstats
```

### MATLAB/Octave

From a cloned repository, add the MATLAB package folder to the path:

```matlab
addpath('/path/to/statgen/matlab')
```

Alternatively, download the MATLAB/Octave bundle from a GitHub release and add
that bundle's `matlab` folder to the path.
