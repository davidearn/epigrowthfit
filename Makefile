SHELL = bash
MAKE = make
R = R
PACKAGE = epigrowthfit
VERSION := $(shell sed -n '/^Version: /s/Version: // p' DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
MANUAL := $(PACKAGE)-manual.pdf

all: clean install

install: deps build
	export NOT_CRAN=true
	$(R) CMD INSTALL ../$(TARBALL)

deps:
	$(R) --no-echo -e 'devtools::install_deps(".", dependencies = TRUE)'

build: $(TARBALL)

$(TARBALL): rd DESCRIPTION src/$(PACKAGE).so
	$(R) CMD build .
	mv $@ ..

rd: R/*.R inst/REFERENCES.bib
	$(R) --no-echo -e 'devtools::document(".")'
	sed -i 's/USCORE/_/ g' man/*.Rd
	sed -i 's/\\text{/\\textrm{/ g' man/*.Rd

src/$(PACKAGE).so:
	$(MAKE) -C src

manual: $(MANUAL)

$(MANUAL): man/*.Rd
	$(R) CMD Rd2pdf --force --output=$@ --no-preview .
	mv $@ ..

check:
	$(R) --no-echo -e 'devtools::check(".")'

test:
	$(R) --no-echo -e 'devtools::test(".")'

clean:
	rm -f ../$(TARBALL) ../$(PACKAGE)-manual.pdf
	find . \( -name "#*" -o -name "*~" -o -name ".Rhistory" \) \
		-exec rm {} +
	$(MAKE) -C vignettes clean
	$(MAKE) -C src clean
	$(MAKE) -C scrap clean
