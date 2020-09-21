SHELL = bash
MAKE = make
R = R
PACKAGE = epigrowthfitTMB
VERSION := $(shell sed -n '/^Version: /s/Version: // p' DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz

all: clean install

fancy: clean install check

basic:
	$(R) --no-echo -e 'devtools::document(".")'
	$(R) CMD build --no-manual .
	$(R) CMD INSTALL $(TARBALL)

install: build
	export NOT_CRAN=true
	$(R) CMD INSTALL ../$(TARBALL)

build: $(TARBALL)

$(TARBALL): R/*.R rd DESCRIPTION NAMESPACE
	$(R) CMD build --no-manual .
	mv $@ ..

rd: man/*.Rd

man/*.Rd: R/*.R inst/REFERENCES.bib
	$(R) --no-echo -e 'devtools::document(".", roclets = "rd")'
	touch man/*.Rd

DESCRIPTION: R/*.R
	$(R) --no-echo -e 'devtools::document(".", roclets = "collate")'
	touch $@

NAMESPACE: R/*.R
	$(R) --no-echo -e 'devtools::document(".", roclets = "namespace")'
	touch $@

check:
	$(R) --no-echo -e 'devtools::check(".")'

test:
	$(R) --no-echo -e 'devtools::test(".")'

clean:
	rm -f $(TARBALL) ../$(TARBALL)
	find . \( -name "#*" -o -name "*~" -o -name ".Rhistory" \) \
		-exec rm {} +

new: clean
