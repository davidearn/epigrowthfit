(path_utils <- Sys.getenv("PATH_UTILS"))
utils <- new.env()
source(path_utils, local = utils)

(path_rds <- Sys.getenv("PATH_RDS"))
windows <- readRDS(path_rds)

datasets <- c("covid", "coords", "population", "mobility", "npi", "vaccine", "devel", "equity")
for (s in paste0("path_", datasets)) {
  print(assign(s, Sys.getenv(toupper(s))))
  assign(sub("path_", "", s), readRDS(get(s)))
}

## ISO 3166-1 alpha-3 code
country_iso_alpha3 <- levels(windows[["country_iso_alpha3"]])
n <- length(country_iso_alpha3)
i <- unclass(windows[["country_iso_alpha3"]])

## Country code, ISO 3166-1 alpha-3
country_iso_alpha2 <- countrycode::countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "iso2c"
)
windows[["country_iso_alpha2"]]  <- gl(n, 1L, labels = country_iso_alpha2)[i]

## Country code, ISO 3166-1 numeric
country_iso_numeric <- countrycode::countrycode(
  sourcevar = al3,
  origin = "iso3c",
  destination = "iso3n"
)
windows[["country_iso_numeric"]] <- gl(n, 1L, labels = country_iso_numeric)[i]

## Country name
country <- countrycode::countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "country.name",
  custom_match = c(
    BIH = "Bosnia and Herzegovina",
    CIV = "Côte d'Ivoire",
    MMR = "Myanmar",
    PSE = "Palestine",
    TTO = "Trinidad and Tobago"
  )
)
windows[["country"]] <- gl(n, 1L, labels = country)[i]

## Subregion
subregion <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = c(TWN = "Eastern Asia")
)
windows[["subregion"]] <- factor(subregion)[i]

## Region
region <- countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = c(TWN = "Asia")
)
windows[["region"]] <- factor(region)[i]

## City population-weighted latitude and longitude
coords[["country_iso_alpha3"]] <- factor(coords[["country_iso_alpha3"]], levels = country_iso_alpha3)
f <- function(d) {
  c(weighted.mean(d[["longitude"]], d[["population"]]),
    weighted.mean(d[["latitude"]],  d[["population"]]))
}
longitude_latitude_rows <- c(by(coords, coords[["country_iso_alpha3"]], f, simplify = FALSE))
longitude_latitude <- matrix(unlist(ll_rows), ncol = 2L, byrow = TRUE)
windows[c("longitude", "latitude")] <- ll[i, , drop = FALSE]

## Longitude and latitude bands
longitude_band_20deg <- cut(
  x = ll[, 1L],
  breaks = seq.int(0, 360, by = 20),
  right = FALSE,
  include.lowest = FALSE,
  ordered_result = TRUE
)
latitude_band_20deg <- cut(
  x = ll[, 2L],
  breaks = seq.int(-90, 90, by = 20),
  right = FALSE,
  include.lowest = TRUE,
  ordered_result = TRUE
)
windows[["longitude_band_20deg"]] <- longitude_band_20deg[i]
windows[["latitude_band_20deg"]] <- latitude_band_20deg[i]

## Population size
population <- population[population[["year"]] == "2020", , drop = FALSE]
m <- match(country_iso_numeric, population[["country_iso_numeric"]], 0L)
windows[["population"]] <- population[["population"]][m][i]

## Mobility
s <- grep("^mobility_", names(mobility), value = TRUE)
mobility[["country_iso_alpha2"]] <- factor(mobility[["country_iso_alpha2"]], levels = country_iso_alpha2)
M <- get_summary(
  d = mobility[s],
  d_Date = mobility$Date,
  d_index = mobility$country_iso_alpha2,
  start = windows$start,
  start_index = windows$country_iso_alpha2,
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
windows[s] <- M[, s]
rm(mobility)

## NPI
npi <- readRDS("rds/npi.rds")
s <- grep("^npi_", names(npi), value = TRUE)
s_ordered <- grep("^npi_(flag|indic_(?!(E3|E4|H4|H5))).*$", names(npi), value = TRUE, perl = TRUE)
s_numeric <- setdiff(s, s_ordered)
npi$country_iso_alpha3 <- factor(npi$country_iso_alpha3, levels = country_iso_alpha3)
npi[s_ordered] <- lapply(npi[s_ordered], function(x) as.integer(as.character(x)))
M_ordered <- get_summary(
  d = npi[s_ordered],
  d_Date = npi$Date,
  d_index = npi$country_iso_alpha3,
  start = windows$start,
  start_index = windows$country_iso_alpha3,
  func = function(x) min(stat_mode(x, na.rm = TRUE)),
  lag = 14L,
  k = 14L,
  method = "locf",
  x0 = 0L
)
windows[s_ordered] <- lapply(as.data.frame(M_ordered[, s_ordered]), ordered)
M <- get_summary(
  d = npi[s_numeric],
  d_Date = npi$Date,
  d_index = npi$country_iso_alpha3,
  start = windows$start,
  start_index = windows$country_iso_alpha3,
  func = mean,
  na.rm = TRUE,
  lag = 14L,
  k = 14L,
  method = "linear",
  x0 = NULL,
  x1 = NULL,
  period = 1L,
  geom = FALSE
)
windows[s_numeric] <- M[, s_numeric]
rm(npi)

## Vaccination
vaccine <- readRDS("rds/vaccine.rds")
s <- grep("^vaccinated_", names(vaccine), value = TRUE)
vaccine$country_iso_alpha3 <- factor(vaccine$country_iso_alpha3, levels = country_iso_alpha3)
M <- get_summary(
  d = vaccine[s],
  d_Date = vaccine$Date,
  d_index = vaccine$country_iso_alpha3,
  start = windows$start,
  start_index = windows$country_iso_alpha3,
  func = mean,
  na.rm = TRUE,
  lag = 14L,
  k = 14L,
  method = "linear",
  x0 = 0,
  x1 = NULL,
  period = 1L,
  geom = FALSE
)
windows[sub("_per_100$", "", s)] <- M[, s] / 100
rm(vaccine)

## Weather
weather <- readRDS("rds/weather.rds")
s <- grep("^weather_", names(weather), value = TRUE)
weather$country_iso_alpha3 <- factor(weather$country_iso_alpha3, levels = country_iso_alpha3)
M <- get_summary(
  d = weather[s],
  d_Date = weather$Date,
  d_index = weather$country_iso_alpha3,
  start = windows$start,
  start_index = windows$country_iso_alpha3,
  func = mean,
  na.rm = TRUE,
  lag = 14L,
  k = 14L,
  method = "linear",
  x0 = NULL,
  x1 = NULL,
  period = 1L,
  geom = FALSE
)
windows[s] <- M[, s]
rm(weather)

## Economic indicators
devel <- readRDS("rds/devel.rds")
v <- c(
  econ_gdp_per_capita = "GDP per capita (constant 2010 US$)",
  econ_gini           = "Gini index (World Bank estimate)"
)
devel$country_iso_alpha3 <- factor(devel$country_iso_alpha3, levels = country_iso_alpha3)
devel$indic <- factor(devel$indic, levels = v, labels = names(v))
devel$year <- factor(devel$year, levels = 2001:2020)
i <- complete.cases(devel[c("country_iso_alpha3", "indic", "year")])
devel <- devel[i, , drop = FALSE]
lo <- function(x) {
  if (all(argna <- is.na(x))) NA else x[max(which(!argna))]
}
lodevel <- aggregate(devel["value"], by = devel[c("country_iso_alpha3", "indic")], lo, drop = FALSE)
windows[levels(lodevel$indic)] <- c(tapply(lodevel$value, lodevel$indic, `[`, r, simplify = FALSE))
windows$econ_gini <- windows$econ_gini / 100
rm(devel)

## Number of days since:
## * 50 persons reported infected
## * 1 in 500,000 persons reported infected
covid19 <- readRDS("rds/covid19.rds")
covid19$country_iso_alpha3 <- factor(covid19$country_iso_alpha3, levels = country_iso_alpha3)
f <- function(d, min = 1L, population = NULL) {
  if (!is.null(population)) {
    min <- min * population
  }
  if (is.na(min) | !any(geq <- d$cases_total >= min)) {
    return(.Date(NA_real_))
  }
  d$Date[which.max(geq)]
}
Date_50 <- .Date(c(by(covid19, covid19$country_iso_alpha3, f, min = 50L)))
Date_2in1m <- .Date(mapply(f,
  d = split(covid19, covid19$country_iso_alpha3),
  min = 2e-06,
  population = windows$population[match(country_iso_alpha3, windows$country_iso_alpha3)]
))
windows$days_since_50 <- as.numeric(windows$start - Date_50[r])
windows$days_since_2in1m <- as.numeric(windows$start - Date_2in1m[r])
rm(covid19)

saveRDS(windows, file = path_rds)
