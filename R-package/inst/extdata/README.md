Tiny CRAN-safe fixtures for R package tests are included here in built package
tarballs.

Portable-format fixtures copied into this directory must be small and derived
from the repository-level canonical fixtures under `tests/fixtures/` when
possible. They are generated working-tree copies, not independent committed
sources. Broader cross-runtime parity tests belong in repository-level pytest
tests.

R-native fixtures with no portable counterpart may be committed directly here.

Current copied fixtures:

- `reference_chr1.bim` is a byte-for-byte copy of
  `tests/fixtures/reference/sharded/1.bim`.

Prepare copied R fixtures with:

```sh
make prepare-r-fixtures
```

Repository pytest checks that copied fixtures remain synchronized with their
canonical source files.
