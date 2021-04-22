R = R

all: tidy weather world endpoints

tidy:
	$(R) -f tidy.R

metoffice:
	$(R) -f populate_weather.R
	$(R) -f merge_weather.R

world:
	$(R) -f make_world.R

endpoints:
	$(R) -f make_endpoints.R
	$(R) -f augment_endpoints.R

