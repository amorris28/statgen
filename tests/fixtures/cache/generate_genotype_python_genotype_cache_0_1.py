#!/usr/bin/env python3
"""Generate tests/fixtures/cache/genotype_python_genotype_cache_0_1.npz.

Run from the repository root:
    python tests/fixtures/cache/generate_genotype_python_genotype_cache_0_1.py
"""

from cache_fixture_helpers import generate_python_cache_fixture

generate_python_cache_fixture("genotype")
