SHELL := /bin/bash
MAKE := make
R := R

dir_R = R
dir_csv = csv
dir_weather = weather

datasets = covid coords population mobility npi devel equity
csv = $(addprefix $(dir_csv)/,$(addsuffix .csv,$(datasets)))
rda = $(addsuffix .RData,$(datasets))

all: world endpoints

.PHONY: weather clean

world: $(dir_R)/make_world.R covid.RData
	$(R) --quiet -f $<

endpoints: $(addprefix $(dir_R)/,$(addsuffix .R,make_endpoints augment_endpoints utils)) \
  $(rda) weather.RData
	$(R) --quiet -f $(word 1, $^)
	$(R) --quiet -f $(word 2, $^)

$(rda): $(dir_R)/tidy.R $(csv)
	$(R) --quiet -f $<

weather.RData: $(dir_R)/make_weather.R coords.RData weather
	$(R) --quiet -f $<

weather:
	mkdir -p $(dir_weather)
	$(R) --quiet -f $(dir_R)/get_weather.R >& get_weather.Rout

new: clean
	rm -rf $(dir_csv) $(dir_weather)

clean:
	rm -f *.{RData,Rout}

download:
	mkdir -p $(dir_csv)
	# https://github.com/CSSEGISandData/COVID-19
	wget -O $(dir_csv)/covid.csv https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv
	# https://simplemaps.com/data/world-cities
	wget -O $(dir_csv)/tmp.zip https://simplemaps.com/static/data/world-cities/basic/simplemaps_worldcities_basicv1.73.zip
	unzip -p $(dir_csv)/tmp.zip worldcities.csv > $(dir_csv)/coords.csv
	rm $(dir_csv)/tmp.zip
	# https://population.un.org/wpp/Download/Standard/CSV/
	wget -O $(dir_csv)/population.csv https://population.un.org/wpp/Download/Files/1_Indicators%20\(Standard\)/CSV_FILES/WPP2019_TotalPopulationBySex.csv
	# https://www.google.com/covid19/mobility/
	wget -O $(dir_csv)/mobility.csv https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv
	# https://covidtracker.bsg.ox.ac.uk/
	wget -O $(dir_csv)/npi.csv https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest.csv
	# https://datacatalog.worldbank.org/dataset/world-development-indicators
	wget -O $(dir_csv)/tmp.zip https://databank.worldbank.org/data/download/WDI_csv.zip
	unzip -p $(dir_csv)/tmp.zip WDIData.csv > $(dir_csv)/devel.csv
	rm $(dir_csv)/tmp.zip
	# https://datacatalog.worldbank.org/dataset/poverty-and-equity-database
	wget -O $(dir_csv)/tmp.zip https://databank.worldbank.org/data/download/PovStats_csv.zip
	unzip -p $(dir_csv)/tmp.zip PovStatsData.csv > $(dir_csv)/equity.csv
	rm $(dir_csv)/tmp.zip
