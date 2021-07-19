MAKE := make
R := R
PACKAGE := epigrowthfit
VERSION := $(shell sed -n '/^Version: /s/Version: // p' DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
MANUAL := $(PACKAGE)-manual.pdf
CHECKDIR := $(PACKAGE).Rcheck

all: clean install
.PHONY = deps 

install: install-deps build
	export NOT_CRAN=true
	$(R) CMD INSTALL $(TARBALL)

install-deps:
	$(R) --quiet -e 'devtools::install_deps(".")'

build: $(TARBALL)

$(TARBALL): enums docs src/$(PACKAGE).cpp src/*.h DESCRIPTION NAMESPACE
	$(R) CMD build --no-manual .

enums: utils/update_enums.R src/enums.h
	cd $(dir $<) && $(R) --quiet -f $(notdir $<)

docs: R/*.R
	$(R) --quiet -e 'devtools::document(".")'

manual: $(MANUAL)

$(MANUAL): enum docs
	$(R) CMD Rd2pdf -o $@ --force --no-preview .

check: build
	$(R) CMD check --no-tests $(TARBALL)

test:
	$(R) --quiet -e 'devtools::test(".")'

clean:
	rm -fr $(TARBALL) $(PACKAGE)-manual.pdf $(CHECKDIR)
	find . \( -name "#*" -o -name "*~" -o -name ".Rhistory" \) \
		-exec rm {} +
	rm -f src/*.{o,so}
