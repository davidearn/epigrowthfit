R := R
PACKAGE := epigrowthfit
VERSION := $(shell sed -n "/^Version: /s/Version: // p" DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
MANUAL := $(PACKAGE)-manual.pdf
CHECKDIR := $(PACKAGE).Rcheck

all: clean install
.PHONY = clean install-deps

build: $(TARBALL)

install: $(TARBALL)
	export NOT_CRAN=true
	$(R) CMD INSTALL --preclean $<

install-deps:
	$(R) --quiet -e "devtools::install_deps(\".\")"

$(TARBALL): update-enums update-docs test-setup src/*.cpp DESCRIPTION NAMESPACE
	$(R) CMD build --no-manual .

pdf: $(MANUAL)

$(MANUAL): update-enums update-docs
	$(R) CMD Rd2pdf -o $@ --force --no-preview .

update-enums: utils/update_enums.R src/enums.h
	cd $(dir $<) && $(R) --quiet -f $(notdir $<)

update-docs: R/*.R
	$(R) --quiet -e "devtools::document(\".\")"

test: test-setup
	$(R) --quiet -e "devtools::test(\".\")"

test-setup: src/*.h
	cp $^ inst/testsrc

check: $(TARBALL)
	$(R) CMD check --no-tests $<

check-as-cran: $(TARBALL) test-setup
	$(R) CMD check --as-cran $<

clean:
	rm -fr $(TARBALL) $(MANUAL) $(CHECKDIR)
	find . \( -name "#*" -o -name "*~" \) -exec rm {} +
	rm -f src/*.{o,so,tmp} tests/testthat/*.{o,so,tmp} inst/testsrc/*.h
