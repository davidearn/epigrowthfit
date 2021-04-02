cities <- read.csv("worldcities.csv")
cities <- cities[c("city", "iso3", "population", "lat", "lng")]
cities <- cities[complete.cases(cities), , drop = FALSE]
save(cities, file = "cities.RData")
