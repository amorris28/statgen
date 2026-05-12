R_PACKAGE_VERSION := $(shell awk '/^Version:/ {print $$2}' R-package/DESCRIPTION)
R_PACKAGE_TARBALL := /tmp/statgen_$(R_PACKAGE_VERSION).tar.gz

.PHONY: install fixtures prepare-r-fixtures test test-python test-octave test-matlab test-r rcmd-check

install:
	pip install -e python/

fixtures:
	python tests/fixtures/generate.py
	$(MAKE) prepare-r-fixtures

prepare-r-fixtures:
	mkdir -p R-package/inst/extdata
	cp tests/fixtures/reference/sharded/1.bim R-package/inst/extdata/reference_chr1.bim

test: prepare-r-fixtures
	pytest tests/

test-python: prepare-r-fixtures
	pytest tests/ -m "not octave"

test-octave:
	pytest tests/ -m octave

test-matlab:
	STATGEN_MATLAB=1 pytest tests/ -m octave

test-r: prepare-r-fixtures
	pytest tests/ -m r

rcmd-check: prepare-r-fixtures
	cd /tmp && R CMD build $(CURDIR)/R-package
	R CMD check --as-cran --output=/tmp $(R_PACKAGE_TARBALL)
