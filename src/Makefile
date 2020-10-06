PACKAGE = epigrowthfit

all: $(PACKAGE).so

$(PACKAGE).so: $(PACKAGE).cpp
	R --no-echo -e "TMB::compile('epigrowthfit.cpp', '-O0 -g', safebounds = FALSE, safeunload = FALSE)"  # CXX_FLAGS = "-O0 -g"
