R = R

all: tidy world endpoints

tidy:
	$(R) -f tidy.R

world:
	$(R) -f make_world.R

endpoints:
	$(R) -f make_endpoints.R
	$(R) -f augment_endpoints.R

