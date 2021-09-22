R := R
PACKAGE := epigrowthfit
VERSION := $(shell sed -n "/^Version: /s/Version: // p" DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
MANUAL := $(PACKAGE)-manual.pdf
CHECKDIR := $(PACKAGE).Rcheck

all: clean install
.PHONY = clean test install-deps

build: $(TARBALL)

install: $(TARBALL)
	export NOT_CRAN=true
	$(R) CMD INSTALL --preclean $<

install-deps:
	$(R) --quiet -e "devtools::install_deps(\".\")"

$(TARBALL): enums docs src/*.cpp src/*.h DESCRIPTION NAMESPACE
	$(R) CMD build --no-manual .

manual: $(MANUAL)

$(MANUAL): enums docs
	$(R) CMD Rd2pdf -o $@ --force --no-preview .

enums: utils/update_enums.R src/enums.h
	cd $(dir $<) && $(R) --quiet -f $(notdir $<)

docs: R/*.R
	$(R) --quiet -e "devtools::document(\".\")"

check: $(TARBALL)
	$(R) CMD check --no-tests $<

check-as-cran: $(TARBALL)
	$(R) CMD check --as-cran --no-tests $<

test:
	$(R) --quiet -e "devtools::test(\".\")"

clean:
	rm -fr $(TARBALL) $(MANUAL) $(CHECKDIR)
	find . \( -name "#*" -o -name "*~" \) -exec rm {} +
	rm -f src/*.{o,so,tmp} tests/testthat/*.{o,so,tmp}
