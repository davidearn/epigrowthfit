all: download publish

download:
	wget -O tmp.csv https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv
	R -e 'canadacovid_raw <- read.csv("tmp.csv"); save(canadacovid_raw, file = "canadacovid_raw.RData")'
	rm tmp.csv

publish: publish.R *.RData
	R -f $<
