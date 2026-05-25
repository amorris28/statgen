# Cache fixture generation

This directory stores committed runtime-native cache fixtures used to validate
cache-loader compatibility. These files are not portable interchange formats;
they are language-specific accelerators generated from the canonical source
fixtures under `tests/fixtures/`.

Cache fixture generators are intentionally separate from
`tests/fixtures/generate.py` and are not run by `make fixtures`.

Run all current-schema generators from the repository root:

```sh
python tests/fixtures/cache/generate_reference_python_reference_cache_0_1.py
python tests/fixtures/cache/generate_sumstats_python_sumstats_cache_0_1.py
python tests/fixtures/cache/generate_annotations_python_annotations_cache_0_2.py
python tests/fixtures/cache/generate_genotype_python_genotype_cache_0_1.py

Rscript tests/fixtures/cache/generate_reference_r_reference_cache_0_1.R
Rscript tests/fixtures/cache/generate_sumstats_r_sumstats_cache_0_1.R
Rscript tests/fixtures/cache/generate_annotations_r_annotations_cache_0_2.R
Rscript tests/fixtures/cache/generate_genotype_r_genotype_cache_0_1.R

octave --no-gui --quiet tests/fixtures/cache/generate_reference_matlab_reference_cache_0_1.m
octave --no-gui --quiet tests/fixtures/cache/generate_sumstats_matlab_sumstats_cache_0_1.m
octave --no-gui --quiet tests/fixtures/cache/generate_annotations_matlab_annotations_cache_0_2.m
octave --no-gui --quiet tests/fixtures/cache/generate_genotype_matlab_genotype_cache_0_1.m
```

The `matlab` runtime fixtures are written by Octave in normal CI because
Octave is the standard MATLAB compatibility proxy. The generators write MATLAB
`.mat` caches with `format = "v7"`. When MATLAB is available, run the cache
fixture tests with `STATGEN_MATLAB=1` to validate MATLAB can read the same
Octave-written fixtures.

Legacy schema generators may require old release artifacts or an old git tag.
Do not regenerate legacy fixtures casually; those bytes are compatibility
anchors for old public cache schemas.
