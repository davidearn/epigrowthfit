all: covid plague

covid: publish_covid.R
	wget -O canadacovid_raw.csv https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv
	R -f $<
	rm canadacovid_raw.csv

plague: publish_plague.R *.RData
	R -f $<

