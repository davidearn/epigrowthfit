MAKE = make
CSV = csv

all: download tidy weather world endpoints

download: setup covid coords population mobility npi devel equity

tidy:
	$(MAKE) -C R tidy

metoffice:
	$(MAKE) -C R metoffice

world:
	$(MAKE) -C R world

endpoints:
	$(MAKE) -C R endpoints

clean:
	rm -r $(CSV) *.RData

setup:
	mkdir -p $(CSV)

# https://github.com/CSSEGISandData/COVID-19
covid:
	wget -O $(CSV)/covid.csv https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv

# https://simplemaps.com/data/world-cities
coords:
	wget -O $(CSV)/tmp.zip https://simplemaps.com/static/data/world-cities/basic/simplemaps_worldcities_basicv1.73.zip
	unzip -p $(CSV)/tmp.zip worldcities.csv > $(CSV)/coords.csv
	rm $(CSV)/tmp.zip

# https://population.un.org/wpp/Download/Standard/$(CSV)/
population:
	wget -O $(CSV)/population.csv https://population.un.org/wpp/Download/Files/1_Indicators%20\(Standard\)/CSV_FILES/WPP2019_TotalPopulationBySex.csv

# https://www.google.com/covid19/mobility/
mobility:
	wget -O $(CSV)/mobility.csv https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv

# https://covidtracker.bsg.ox.ac.uk/
npi:
	wget -O $(CSV)/npi.csv https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest.csv

# https://datacatalog.worldbank.org/dataset/world-development-indicators
devel:
	wget -O $(CSV)/tmp.zip https://databank.worldbank.org/data/download/WDI_csv.zip
	unzip -p $(CSV)/tmp.zip WDIData.csv > $(CSV)/devel.csv
	rm $(CSV)/tmp.zip

# https://datacatalog.worldbank.org/dataset/poverty-and-equity-database
equity:
	wget -O $(CSV)/tmp.zip https://databank.worldbank.org/data/download/PovStats_csv.zip
	unzip -p $(CSV)/tmp.zip PovStatsData.csv > $(CSV)/equity.csv
	rm $(CSV)/tmp.zip
