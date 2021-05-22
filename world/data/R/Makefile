R = R

all: tidy weather world endpoints

tidy:
	$(R) --quiet -f tidy.R

metoffice:
	$(R) --quiet -f populate_weather.R
	$(R) --quiet -f merge_weather.R

world:
	$(R) --quiet -f make_world.R

endpoints:
	$(R) --quiet -f make_endpoints.R
	$(R) --quiet -f augment_endpoints.R

