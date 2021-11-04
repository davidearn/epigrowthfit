SHELL := /bin/bash
MAKE := make
R := R

DIR_CSV := csv
DIR_RDS := rds

DS_COVID := covid
DS_COORDS := coords
DS_POPULATION := population
DS_MOBILITY := mobility
DS_NPI := npi
DS_VACCINE := vaccine
DS_DEVEL := devel
DS_EQUITY := equity
DS_WEATHER := weather
DS_WORLD := world
DS_WINDOWS := windows
DS_WINDOWS_DEPS := $(filter-out %_WORLD %_WINDOWS,$(filter DS_%,$(.VARIABLES)))

all: world windows

weather: $(DIR_RDS)/$(DS_WEATHER).rds

world: $(DIR_RDS)/$(DS_WORLD).rds

windows: $(DIR_RDS)/$(DS_WINDOWS).rds

$(DIR_RDS):
	mkdir -p $@

$(DIR_RDS)/%.rds: tidy_%.R $(DIR_RDS) $(DIR_CSV)/%.csv
	$(R) --quiet -f $< >& $(<:.R=.Rout)
$(DIR_RDS)/%.rds: export PATH_CSV = $(lastword $^)
$(DIR_RDS)/%.rds: export PATH_RDS = $@

$(DIR_RDS)/$(DS_WEATHER).rds: make_$(DS_WEATHER).R utils_$(DS_WEATHER).R $(DIR_RDS)/$(DS_COORDS).rds
	mkdir -p $(DIR_WEATHER)
	$(R) --quiet -f $< >& $(<:.R=.Rout)
$(DIR_RDS)/$(DS_WEATHER).rds: export DIR_WEATHER = weather
$(DIR_RDS)/$(DS_WEATHER).rds: export PATH_UTILS = $(word 2,$^)
$(DIR_RDS)/$(DS_WEATHER).rds: export PATH_COORDS = $(word 3,$^)
$(DIR_RDS)/$(DS_WEATHER).rds: export PATH_RDS = $@

$(DIR_RDS)/$(DS_WORLD).rds: make_$(DS_WORLD).R $(DIR_RDS)/$(DS_COVID).rds
	$(R) --quiet -f $< >& $(<:.R=.Rout)
$(DIR_RDS)/$(DS_WORLD).rds: export PATH_COVID = $(word 2,$^)
$(DIR_RDS)/$(DS_WORLD).rds: export PATH_RDS = $@

$(DIR_RDS)/$(DS_WINDOWS).rds: make_$(DS_WINDOWS).R augment_$(DS_WINDOWS).R utils_$(DS_WINDOWS).R $(foreach name,$(DS_WINDOWS_DEPS),$(DIR_RDS)/$($(name)).rds) 
	$(R) --quiet -f $< >& $(<:.R=.Rout)
	$(R) --quiet -f $(word 2,$^) >& $((word 2,$^),.R=.Rout)
$(DIR_RDS)/$(DS_WINDOWS).rds: export PATH_UTILS = $(word 3,$^)
$(DIR_RDS)/$(DS_WINDOWS).rds: export PATH_RDS = $@
$(foreach name,$(DS_WINDOWS_DEPS),$(eval $(DIR_RDS)/$(DS_WINDOWS).rds: export $(patsubst DS_%,PATH_%,$(name)) = $(DIR_RDS)/$($(name)).rds))

$(DIR_CSV):
	mkdir -p $@

$(DIR_CSV)/$(DS_COVID).csv: $(DIR_CSV)
	# https://github.com/CSSEGISandData/COVID-19
	wget -O $@ https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv

$(DIR_CSV)/$(DS_COORDS).csv: $(DIR_CSV)
	# https://simplemaps.com/data/world-cities
	wget -O tmp.zip https://simplemaps.com/static/data/world-cities/basic/simplemaps_worldcities_basicv1.73.zip
	unzip -p tmp.zip worldcities.csv > $@
	rm tmp.zip

$(DIR_CSV)/$(DS_POPULATION).csv: $(DIR_CSV)
	# https://population.un.org/wpp/Download/Standard/CSV/
	wget -O $@ https://population.un.org/wpp/Download/Files/1_Indicators%20\(Standard\)/CSV_FILES/WPP2019_TotalPopulationBySex.csv

$(DIR_CSV)/$(DS_MOBILITY).csv: $(DIR_CSV)
	# https://www.google.com/covid19/mobility/
	wget -O $@ https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv

$(DIR_CSV)/$(DS_NPI).csv: $(DIR_CSV)
	# https://covidtracker.bsg.ox.ac.uk/
	wget -O $@ https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest.csv

$(DIR_CSV)/$(DS_VACCINE).csv: $(DIR_CSV)
	# https://github.com/owid/covid-19-data/raw/master/public/data/vaccinations/vaccinations.csv
	wget -O $@ https://github.com/owid/covid-19-data/raw/master/public/data/vaccinations/vaccinations.csv

$(DIR_CSV)/$(DS_DEVEL).csv: $(DIR_CSV)
	# https://datacatalog.worldbank.org/dataset/world-development-indicators
	wget -O tmp.zip https://databank.worldbank.org/data/download/WDI_csv.zip
	unzip -p tmp.zip WDIData.csv > $@
	rm tmp.zip

$(DIR_CSV)/$(DS_EQUITY).csv: $(DIR_CSV)
	# https://datacatalog.worldbank.org/dataset/poverty-and-equity-database
	wget -O tmp.zip https://databank.worldbank.org/data/download/PovStats_csv.zip
	unzip -p tmp.zip PovStatsData.csv > $@
	rm tmp.zip

clean:
	rm -rf $(DIR_RDS) *.Rout *~

new: clean
	rm -rf $(DIR_CSV) $(DIR_WEATHER)
