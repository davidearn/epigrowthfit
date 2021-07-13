all: mre.Rout

%.Rout: %.R %.cpp %.RData
	R CMD BATCH $<

clean:
	rm -f *.{o,so,Rout}
