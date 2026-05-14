"""R phase 0 acceptance tests: package skeleton and optional R smoke tests."""

import json
import os
import subprocess

import pytest

from tests.conftest import (
    FIXTURES_DIR,
    R_PACKAGE_DIR,
    REPO_ROOT,
    run_rscript,
    skipif_no_rscript,
)


def _r_package_version() -> str:
    for line in (R_PACKAGE_DIR / "DESCRIPTION").read_text().splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    raise AssertionError("R package DESCRIPTION is missing Version")


def test_r_package_skeleton_present():
    required = [
        "DESCRIPTION",
        "NAMESPACE",
        "R/version.R",
        "R/verbosity.R",
        "man/statgen-package.Rd",
        "man/version.Rd",
        "man/verbosity.Rd",
        "tests/testthat.R",
        "tests/testthat/test-ld.R",
        "tests/testthat/test-smoke.R",
        "inst/extdata/README.md",
    ]
    for rel in required:
        assert (R_PACKAGE_DIR / rel).is_file(), f"missing R package file: {rel}"


def test_r_extdata_reference_fixture_matches_canonical_source():
    copies = {
        "reference_chr1.bim": "reference/sharded/1.bim",
        "reference_chrX.bim": "reference/sharded/X.bim",
        "traits.tsv.gz": "sumstats/traits.tsv.gz",
        "anno1.bed": "annotations/anno1.bed",
        "anno2.bed": "annotations/anno2.bed",
        "genotype_1.bed": "genotype/sharded/1.bed",
        "genotype_1.bim": "genotype/sharded/1.bim",
        "genotype_1.fam": "genotype/sharded/1.fam",
        "genotype_X.bed": "genotype/sharded/X.bed",
        "genotype_X.bim": "genotype/sharded/X.bim",
        "genotype_X.fam": "genotype/sharded/X.fam",
        "genotype_X.ploidy": "genotype/sharded/X.ploidy",
        "ld/python/ld_manifest.json": "ld/python/ld_manifest.json",
        "ld/python/reference_cache.npz": "ld/python/reference_cache.npz",
        "ld/python/reference_chr1.bim": "ld/python/reference_chr1.bim",
        "ld/python/reference_chrX.bim": "ld/python/reference_chrX.bim",
        "ld/python/ld_chr1.npz": "ld/python/ld_chr1.npz",
        "ld/python/ld_chrX_female.npz": "ld/python/ld_chrX_female.npz",
        "ld/python/ld_chrX_male.npz": "ld/python/ld_chrX_male.npz",
        "ld/python/ld_chrX_combined.npz": "ld/python/ld_chrX_combined.npz",
    }
    for r_name, canonical_rel in copies.items():
        r_fixture = R_PACKAGE_DIR / f"inst/extdata/{r_name}"
        canonical = FIXTURES_DIR / canonical_rel
        assert r_fixture.is_file(), "run `make prepare-r-fixtures` to create R extdata copies"
        assert r_fixture.read_bytes() == canonical.read_bytes()


def test_r_and_python_versions_match():
    import statgen

    assert _r_package_version() == statgen.__version__


def test_r_package_license_matches_repository_license():
    root_license = (REPO_ROOT / "LICENSE").read_text()
    r_license = (R_PACKAGE_DIR / "LICENSE").read_text().splitlines()

    assert "MIT License" in root_license
    assert "Copyright (c) 2026 Oleksandr Frei and contributors" in root_license
    assert "YEAR: 2026" in r_license
    assert "COPYRIGHT HOLDER: Oleksandr Frei and contributors" in r_license


@pytest.mark.r
@skipif_no_rscript
def test_rscript_phase0_source_smoke():
    script = (
        f"source({str(R_PACKAGE_DIR / 'R/verbosity.R')!r}); "
        "stopifnot(identical(get_verbosity(), 'info')); "
        "set_verbosity('quiet'); "
        "stopifnot(identical(get_verbosity(), 'quiet')); "
        "ok <- FALSE; tryCatch(set_verbosity('verbose'), error = function(e) ok <<- TRUE); "
        "stopifnot(ok)"
    )
    result = run_rscript(script)
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_package_install_and_library_smoke(tmp_path):
    deps = run_rscript(
        "pkgs <- c('Matrix', 'jsonlite', 'digest', 'bit64', 'testthat'); "
        "missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; "
        "cat(paste(missing, collapse = ','))"
    )
    assert deps.returncode == 0, deps.stderr
    missing = deps.stdout.strip()
    if missing:
        pytest.skip(f"R package dependencies unavailable: {missing}")

    lib = tmp_path / "r-lib"
    lib.mkdir()
    install = subprocess.run(
        ["R", "CMD", "INSTALL", "-l", str(lib), str(R_PACKAGE_DIR)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert install.returncode == 0, install.stderr

    env = os.environ.copy()
    existing_libs = env.get("R_LIBS")
    env["R_LIBS"] = str(lib) if not existing_libs else str(lib) + os.pathsep + existing_libs
    expected_version = _r_package_version()
    result = run_rscript(
        "library(statgen); "
        f"stopifnot(identical(version(), {expected_version!r})); "
        "set_verbosity('quiet'); "
        "stopifnot(identical(get_verbosity(), 'quiet'))",
        env=env,
    )
    assert result.returncode == 0, result.stderr

    testthat_dir = R_PACKAGE_DIR / "tests/testthat"
    result = run_rscript(
        "library(testthat); "
        "library(statgen); "
        f"testthat::test_dir({json.dumps(str(testthat_dir))}, reporter = 'summary')",
        env=env,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
