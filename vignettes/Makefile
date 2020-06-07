Rscript = $(R_HOME)/bin/Rscript

all: plotmodels.pdf epigrowthfit-intro.pdf

plotmodels.tex: plotmodels.R
	R CMD BATCH --vanilla $<

%.pdf: %.Rnw
	$(Rscript) -e "knitr::knit2pdf('$<')"

%.pdf: %.tex
	pdflatex $<
