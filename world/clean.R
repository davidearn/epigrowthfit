## Retrieve full DELVE COVID-19 data set
loc <- paste0(
  "https://raw.githubusercontent.com/",
  "rs-delve/covid19_datasets/master/dataset/",
  "combined_dataset_latest.csv"
)
world_full <- read.csv(url(loc))
world <- world_full[c("DATE", "cases_new", "ISO", "stats_population")]
names(world) <- c("date", "cases", "country_iso3", "population")

## Clean
world$date <- as.Date(world$date)
world$cases[world$cases < 0] <- NA
world$country_iso3 <- factor(world$country_iso3)
f <- function(d) {
  if (sum(d$cases, na.rm = TRUE) > 1000) {
    i <- seq.int(which.min(is.na(d$cases)), nrow(d))
  } else {
    i <- integer(0L)
  }
  d[i, , drop = FALSE]
}
world <- droplevels(do.call(rbind, lapply(split(world, world$country_iso3), f)))
row.names(world) <- NULL

## Retrieve country name, region and continent from ISO codes
library("countrycode")
country_iso3 <- levels(world$country_iso3)
country_name <- factor(countrycode(
  sourcevar = country_iso3,
  origin = "iso3c",
  destination = "country.name"
))
region <- factor(countrycode(
  sourcevar = country_iso3,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = c(TWN = "Eastern Asia")
))
continent <- factor(countrycode(
  sourcevar = country_iso3,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = c(TWN = "Asia")
))

## Retrieve population sizes
population <- tapply(world$population, world$country_iso3, `[`, 1L)

## Compute population-weighted average latitude and longitude of cities
load("cities/cities.RData")
f <- function(d) {
  c(latitude  = weighted.mean(d$lat, d$population),
    longitude = weighted.mean(d$lng, d$population))
}
ll <- lapply(split(cities, cities$iso3), f)
m <- match("Jerusalem", cities$city, 0L)
ll$PSE <- c(latitude = cities$lat[m], longitude = cities$lng[m])
ll <- do.call(rbind, ll[country_iso3])

## Save time series
world <- world[c("date", "cases", "country_iso3")]
save(world, file = "world.RData")

## Save time invariant variables in a separate data frame
worldstats <- data.frame(
  country_iso3,
  country_name,
  region,
  continent,
  population,
  ll,
  row.names = NULL,
  stringsAsFactors = TRUE
)
save(worldstats, file = "worldstats.RData")
