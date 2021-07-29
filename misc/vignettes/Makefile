SHELL = bash
R = R

vignettes = \
	model
pdf = $(addsuffix .pdf,$(vignettes))
tex = $(addsuffix .tex,$(vignettes))
bib = epigrowthfit

build: $(pdf)
	#$(R) --no-echo -e 'tools::compactPDF("$<", gs_quality = "ebook")'
	rm -f $(addsuffix .tex,$(vignettes))
	rm -f *.{aux,log,out,toc} *~

%.pdf: %.tex logo.png
	pdflatex -interaction=batchmode $<
	#bibtex $(bib)
	pdflatex -interaction=batchmode $<
	pdflatex $<

%.tex: %.Rnw preamble.tex knit_config.R knit_theme.css
	$(R) --no-echo -e 'knitr::knit("$<")'

all: $(pdf)

new: clean
	rm -f $(pdf)

clean:
	rm -f $(tex)
	rm -f *.{aux,log,out,toc} *~
