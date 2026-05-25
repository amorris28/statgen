#!/usr/bin/env python3
"""Generate tests/fixtures/cache/annotations_python_annotations_cache_0_1.npz.

Run from a v0.3.2 checkout at the repository root:
    python tests/fixtures/cache/generate_annotations_python_annotations_cache_0_1.py

This uses the official v0.3.2 Python implementation, whose annotation writer
produces annotations_cache/0.1.
"""

from pathlib import Path
import sys

CACHE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CACHE_DIR.parent.parent.parent
PYTHON_ROOT = REPO_ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from statgen.annotations import load_annotations, save_annotations_cache
from statgen.reference import load_reference


def main() -> None:
    reference = load_reference("tests/fixtures/reference/sharded/@.bim")
    annotations = load_annotations(
        [
            "tests/fixtures/annotations/anno1.bed",
            "tests/fixtures/annotations/anno2.bed",
        ],
        reference,
    )
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    save_annotations_cache(
        annotations,
        CACHE_DIR / "annotations_python_annotations_cache_0_1.npz",
    )


if __name__ == "__main__":
    main()
