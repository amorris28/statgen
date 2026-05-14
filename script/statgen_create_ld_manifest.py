#!/usr/bin/env python3
"""Create a statgen Python LD manifest from existing shard files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from statgen._ld_writer import create_ld_npz_manifest


def main(argv=None) -> int:
    args = _parse_args(argv)
    manifest_path = Path(args.ld) / "ld_manifest.json"
    if manifest_path.exists():
        print(
            f"Refusing to overwrite existing LD manifest: {manifest_path}",
            file=sys.stderr,
        )
        return 1
    create_ld_npz_manifest(args.ld, validate=not args.no_validate)
    return 0


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ld", required=True, help="LD distribution directory containing .npz shard files")
    parser.add_argument("--no-validate", action="store_true", help="Write the manifest without validating the full distribution")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
