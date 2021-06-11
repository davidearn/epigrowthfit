MAKE := make
R := R
PACKAGE := epigrowthfit
VERSION := $(shell sed -n '/^Version: /s/Version: // p' DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
MANUAL := $(PACKAGE)-manual.pdf

all: clean install

install: deps build
	export NOT_CRAN=true
	$(R) CMD INSTALL ../$(TARBALL)

deps:
	$(R) --quiet -e 'devtools::install_deps(".")'

build: $(TARBALL)

$(TARBALL): enum rd src/$(PACKAGE).so DESCRIPTION
	$(R) CMD build .
	mv $@ ..

enum: utils/update_enum.R
	$(MAKE) -C utils enum

rd: R/*.R inst/REFERENCES.bib
	$(R) --quiet -e 'devtools::document(".")'
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
	$(R) --quiet -e 'devtools::check(".")'

clean:
	rm -f ../$(TARBALL) ../$(PACKAGE)-manual.pdf
	find . \( -name "#*" -o -name "*~" -o -name ".Rhistory" \) \
		-exec rm {} +
	$(MAKE) -C src clean
	$(MAKE) -C scrap clean
