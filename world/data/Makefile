all: download tidy world endpoints

download: covid coords population mobility npi

tidy:
	R --no-echo -f tidy.R

world:
	R --no-echo -f make_world.R

endpoints:
	R --no-echo -f make_endpoints.R
	R --no-echo -f augment_endpoints.R

clean:
	rm *.csv *.RData

# https://github.com/CSSEGISandData/COVID-19
covid:
	wget -O covid.csv https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv

# https://simplemaps.com/data/world-cities
coords:
	wget -O tmp.zip https://simplemaps.com/static/data/world-cities/basic/simplemaps_worldcities_basicv1.73.zip
	unzip -p tmp.zip worldcities.csv > coords.csv
	rm tmp.zip

# https://population.un.org/wpp/Download/Standard/CSV/
population:
	wget -O population.csv https://population.un.org/wpp/Download/Files/1_Indicators%20\(Standard\)/CSV_FILES/WPP2019_TotalPopulationBySex.csv

# https://www.google.com/covid19/mobility/
mobility:
	wget -O mobility.csv https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv

# https://covidtracker.bsg.ox.ac.uk/
npi:
	wget -O npi.csv https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest.csv
