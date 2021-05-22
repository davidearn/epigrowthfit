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
	$(R) -e 'devtools::install_deps(".", dependencies = TRUE)'

build: $(TARBALL)

$(TARBALL): flags rd src/$(PACKAGE).so DESCRIPTION
	$(R) CMD build .
	mv $@ ..

flags: utils/update_flags.R
	$(R) -e 'setwd("utils"); source("update_flags.R")'

rd: R/*.R inst/REFERENCES.bib
	$(R) -e 'devtools::document(".")'
	sed -i.bak 's/USCORE/_/ g' man/*.Rd 
	sed -i.bak 's/\\text{/\\textrm{/ g' man/*.Rd
	rm man/*.Rd.bak

src/$(PACKAGE).so: src/$(PACKAGE).cpp src/*.h
	$(MAKE) -C src

manual: $(MANUAL)

$(MANUAL): man/*.Rd
	$(R) CMD Rd2pdf --force --output=$@ --no-preview .
	mv $@ ..

check:
	$(R) -e 'devtools::check(".")'

test:
	$(R) -e 'devtools::test(".")'

clean:
	rm -f ../$(TARBALL) ../$(PACKAGE)-manual.pdf
	find . \( -name "#*" -o -name "*~" -o -name ".Rhistory" \) \
		-exec rm {} +
	#$(MAKE) -C vignettes clean
	$(MAKE) -C src clean
	$(MAKE) -C scrap clean
