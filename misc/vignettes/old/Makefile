SHELL = bash
R = R
NAME = epigrowthfit-vignette

build: $(NAME).pdf
	#$(R) --no-echo -e 'tools::compactPDF("$<", gs_quality = "ebook")'
	rm -rf cache figure
	rm -f *.{bbl,blg,log,out,synctex.gz,tex,toc,aux}

$(NAME).pdf: $(NAME).tex $(NAME).bib logo.png
	pdflatex -interaction=batchmode $<
	bibtex $(NAME)
	pdflatex -interaction=batchmode $<
	pdflatex $<

$(NAME).tex: $(NAME).Rnw knit_theme.css
	$(R) --no-echo -e 'knitr::knit("$<")'

all: $(NAME).pdf

new: clean
	rm -rf cache figure $(NAME).{aux,pdf}

clean:
	rm -f *.{bbl,blg,log,out,synctex.gz,tex,toc}
