all: plotmodels.pdf epigrowthfit-intro.pdf clean

plotmodels.tex: plotmodels.R
	R --no-echo --vanilla -f $<
plotmodels.pdf: plotmodels.tex
	pdflatex $<

epigrowthfit-intro.tex: epigrowthfit-intro.Rnw epigrowthfit-intro.bib vignettepreamble.tex
	R --no-echo --vanilla -e "knitr::knit(\"$<\")"

epigrowthfit-intro.pdf: epigrowthfit-intro.tex plotmodels.pdf
	pdflatex -interaction=batchmode $<
	bibtex epigrowthfit-intro
	pdflatex -interaction=batchmode $<
	pdflatex $<

clean:
	rm -f *.{aux,bbl,blg,log,out,toc}
	rm -f plotmodels.{tex,pdf}
	rm -f epigrowthfit-intro.{R,tex}
	rm -rf cache figure
