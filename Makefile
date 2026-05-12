R_PACKAGE_VERSION := $(shell awk '/^Version:/ {print $$2}' R-package/DESCRIPTION)
R_PACKAGE_TARBALL := /tmp/statgen_$(R_PACKAGE_VERSION).tar.gz

.PHONY: install fixtures prepare-r-fixtures test test-python test-octave test-matlab test-r r-manual rcmd-check

install:
	pip install -e python/

fixtures:
	python tests/fixtures/generate.py
	$(MAKE) prepare-r-fixtures

prepare-r-fixtures:
	mkdir -p R-package/inst/extdata
	cp tests/fixtures/reference/sharded/1.bim R-package/inst/extdata/reference_chr1.bim
	cp tests/fixtures/reference/sharded/X.bim R-package/inst/extdata/reference_chrX.bim
	cp tests/fixtures/sumstats/traits.tsv.gz R-package/inst/extdata/traits.tsv.gz
	cp tests/fixtures/annotations/anno1.bed R-package/inst/extdata/anno1.bed
	cp tests/fixtures/annotations/anno2.bed R-package/inst/extdata/anno2.bed

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

r-manual:
	$(RM) /tmp/statgen-manual.pdf
	R CMD Rd2pdf R-package --output=/tmp/statgen-manual.pdf

rcmd-check: prepare-r-fixtures
	cd /tmp && R CMD build $(CURDIR)/R-package
	R CMD check --as-cran --output=/tmp $(R_PACKAGE_TARBALL)
