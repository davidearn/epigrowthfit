all: gi.pdf models.pdf density_cov.pdf

%.pdf: %.tex
	pdflatex -interaction=batchmode $<
	pdflatex $<

clean:
	rm -f *.out *.log *.aux *~
