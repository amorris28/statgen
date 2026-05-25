#!/usr/bin/env python3
"""Generate tests/fixtures/cache/annotations_python_annotations_cache_0_2.npz.

Run from the repository root:
    python tests/fixtures/cache/generate_annotations_python_annotations_cache_0_2.py
"""

from cache_fixture_helpers import generate_python_cache_fixture

generate_python_cache_fixture("annotations")
