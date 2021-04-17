library("countrycode")
source("utils.R")
load("endpoints.RData")

## ISO 3166-1 alpha-3 code
country_iso_alpha3 <- levels(endpoints$country_iso_alpha3)
n <- length(country_iso_alpha3)
r <- match(endpoints$country_iso_alpha3, country_iso_alpha3, 0L)

## ISO 3166-1 alpha-2 code
country_iso_alpha2 <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "iso2c"
)
endpoints$country_iso_alpha2 <- gl(n, 1L, labels = country_iso_alpha2)[r]

## ISO 3166-1 numeric code
country_iso_numeric <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "iso3n"
)
endpoints$country_iso_numeric <- gl(n, 1L, labels = country_iso_numeric)[r]

## Name
country <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "country.name",
  custom_match = c(MMR = "Myanmar",
                   PSE = "Palestine",
                   TTO = "Trinidad and Tobago")
)
endpoints$country <- gl(n, 1L, labels = country)[r]

## Subregion
subregion <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = c(TWN = "Eastern Asia")
)
endpoints$subregion <- factor(subregion)[r]

## Region
region <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = c(TWN = "Asia")
)
endpoints$region <- factor(region)[r]

## Population-weighted latitude and longitude
load("coords.RData")
coords$country_iso_alpha3 <- factor(coords$country_iso_alpha3, levels = country_iso_alpha3)
f <- function(d) {
  c(weighted.mean(d$latitude,  d$population),
    weighted.mean(d$longitude, d$population))
}
ll_rows <- by(coords, coords$country_iso_alpha3, f, simplify = FALSE)
ll <- matrix(unlist(ll_rows), ncol = 2L, byrow = TRUE)
endpoints[c("latitude", "longitude")] <- ll[r, , drop = FALSE]

## Latitude and longitude bands
latitude_band <- cut(ll[, 1L],
  breaks = seq.int(-90, 90, by = 20),
  include.lowest = TRUE,
  ordered_result = TRUE
)
longitude_band <- cut(ll[, 2L],
  breaks = seq.int(-180, 180, by = 40),
  include.lowest = TRUE,
  ordered_result = TRUE
)
endpoints$latitude_band <- latitude_band[r]
endpoints$longitude_band <- longitude_band[r]

## Population size
load("population.RData")
m <- match(country_iso_numeric, population$country_iso_numeric, 0L)
endpoints$population <- population$population[m][r]

## Mobility
load("mobility.RData")
s <- grep("^mobility_", names(mobility), value = TRUE)
mobility$country_iso_alpha2 <- factor(mobility$country_iso_alpha2, levels = country_iso_alpha2)
M <- get_summary(
  d = mobility[s],
  d_Date = mobility$Date,
  d_index = mobility$country_iso_alpha2,
  start = endpoints$start,
  start_index = endpoints$country_iso_alpha2,
  func = geom_mean,
  na.rm = TRUE,
  zero.rm = FALSE,
  lag = 14L,
  k = 14L,
  method = "linear",
  x0 = 1,
  x1 = NULL,
  period = 7L,
  geom = TRUE
)
endpoints[s] <- M[, s]

## NPI
load("npi.RData")
s <- grep("^npi_", names(npi), value = TRUE)
s_ordered <- grep("^npi_(flag|indic_(?!(E3|E4|H4|H5))).*$", names(npi), value = TRUE, perl = TRUE)
npi$country_iso_alpha3 <- factor(npi$country_iso_alpha3, levels = country_iso_alpha3)
npi[s_ordered] <- lapply(npi[s_ordered], function(x) as.integer(as.character(x)))
M <- get_summary(
  d = npi[s],
  d_Date = npi$Date,
  d_index = npi$country_iso_alpha3,
  start = endpoints$start,
  start_index = endpoints$country_iso_alpha3,
  func = function(x) min(stat_mode(x, na.rm = TRUE)),
  lag = 14L,
  k = 14L,
  method = "locf",
  x0 = 0L
)
endpoints[s] <-  M[, s]
endpoints[s_ordered] <- lapply(endpoints[s_ordered], ordered)

save(endpoints, file = "endpoints.RData")
