all: gi.pdf models.pdf ll.pdf

%.pdf: %.tex
	pdflatex -interaction=batchmode $<
	pdflatex $<

clean:
	rm -f *.out *.log *.aux *~
