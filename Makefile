R := R
PACKAGE := epigrowthfit
VERSION := $(shell sed -n "/^Version: /s/Version: // p" DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
MANUAL := $(PACKAGE)-manual.pdf
CHECKDIR := $(PACKAGE).Rcheck

all: clean install-github

replace-cout-rcout:
	sed -i.bak "s/std::cout/Rcout/" src/*.{cpp,h}
	sed -i.bak "s/\([^R]\)printf/\1Rprintf/" src/*.{cpp,h}
	rm -f src/*.bak

replace-rcout-cout:
	sed -i.bak "s/Rcout/std::cout/" src/*.{cpp,h}
	sed -i.bak "s/Rprintf/printf/" src/*.{cpp,h}
	rm -f src/*.bak

copy-headers:
	cp src/*.h inst/testsrc

update-enums: utils/update_enums.R src/enums.h
	cd $(dir $<) && $(R) --quiet -f $(notdir $<)

update-docs: R/*.R
	$(R) --quiet -e "devtools::document(\".\")"

build-cran: replace-cout-rcout build replace-rcout-cout
	$(MAKE) copy-headers

build: copy-headers update-enums update-docs $(TARBALL)

$(TARBALL):
	$(R) CMD build --no-manual .

install-deps:
	$(R) --quiet -e "devtools::install_deps(\".\")"

install-cran: build-cran install

install-github: build install

install:
	$(R) CMD INSTALL --preclean $(TARBALL)

pdf: $(MANUAL)

$(MANUAL): update-docs
	$(R) CMD Rd2pdf -o $@ --force --no-preview .

check-cran: clean build-cran
	export NOT_CRAN=false && $(R) CMD check --as-cran $(TARBALL)

check-test: clean build
	export NOT_CRAN=true && $(R) CMD check $(TARBALL)

check: clean build
	export NOT_CRAN=true && $(R) CMD check --no-tests $(TARBALL)

clean:
	find . \( -name "#*" -o -name "*~" \) -exec rm {} +
	rm -rf $(TARBALL) $(MANUAL) $(CHECKDIR)
	rm -f src/*.{o,so,tmp}
	rm -f tests/testthat/*.{o,so,tmp}
	rm -rf inst/exdata
