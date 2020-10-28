SHELL = bash
MAKE = make
R = R
PACKAGE = epigrowthfit
VERSION := $(shell sed -n '/^Version: /s/Version: // p' DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz

all: clean install

fancy: clean install check

install: deps build
	export NOT_CRAN=true
	$(R) CMD INSTALL ../$(TARBALL)

deps:
	$(R) --no-echo -e 'devtools::install_deps(".", dependencies = TRUE)'

build: $(TARBALL)

$(TARBALL): src/$(PACKAGE).so R/*.R rd DESCRIPTION NAMESPACE
	$(R) CMD build --no-manual .
	mv $@ ..

src/$(PACKAGE).so:
	$(MAKE) -C src

## FIXME: I'm not sure this is right way  to do this? Can't we just have a single target 'doc-update' (or 'rd')?

rd: man/*.Rd

man/*.Rd: R/*.R inst/REFERENCES.bib
	$(R) --no-echo -e 'devtools::document(".", roclets = "rd")'
	touch man/*.Rd  ## FIXME: why?
	sed -i -e 's/USCORE/_/' man/compute_R0.Rd  ## FIXME do this for all .Rd files?


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
	$(MAKE) -C vignettes clean
	$(MAKE) -C src clean
	$(MAKE) -C scrap clean
