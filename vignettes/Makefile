R = $(R_HOME)/bin/R

all: plotmodels.pdf epigrowthfit-intro.pdf

plotmodels.tex: plotmodels.R
	$(R) CMD BATCH --vanilla $<

%.pdf: %.Rnw *.bib *preamble.tex
	$(R) --no-echo -e "knitr::knit2pdf('$<')"

%.pdf: %.tex
	pdflatex $<
