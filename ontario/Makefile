MACID ?= jaganmn
DATE ?= 2021-04-26
FILENAME = IPHIS_REPORT_2years

all: *.pdf

*.pdf: sandbox.R endpoints
	R -f $<	

endpoints: endpoints.R census.R
	R -f $<

download: tidy.R
	scp -p $(MACID)@ms.mcmaster.ca:/home/earn/$(DATE)/$(FILENAME).CSV.gz .
	gzip -dc $(FILENAME).CSV.gz > $(FILENAME).csv
	R -f $<

new: clean
	rm -f *.{csv,pdf,RData}

clean:
	rm -f *.gz *~
