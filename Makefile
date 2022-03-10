R := $(shell R RHOME)/bin/R
MAKE := make
PACKAGE := $(shell sed -n "/^Package: /s/Package: // p" DESCRIPTION)
VERSION := $(shell sed -n "/^Version: /s/Version: // p" DESCRIPTION)
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
CHECKDIR := $(PACKAGE).Rcheck

all: install

build: $(TARBALL)
$(TARBALL):
	$(R) CMD build .

install: build
	$(R) CMD INSTALL $(TARBALL)

check: build
	$(R) CMD check $(TARBALL)

clean:
	rm -rf $(CHECKDIR) $(TARBALL) *~
	$(MAKE) -C vignettes clean
