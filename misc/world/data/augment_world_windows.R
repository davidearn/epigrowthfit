(path_utils <- Sys.getenv("PATH_UTILS"))
utils <- new.env()
source(path_utils, local = utils)

(path_rds <- Sys.getenv("PATH_RDS"))
windows <- readRDS(path_rds)

datasets <- c("covid", "coords", "population", "mobility", "npi", "vaccine", "devel", "equity", "weather")
for (s in paste0("path_", datasets)) {
  print(assign(s, Sys.getenv(toupper(s))))
  assign(sub("path_", "", s), readRDS(get(s)))
}

## ISO 3166-1 alpha-3 code
country_iso_alpha3 <- levels(windows[["country_iso_alpha3"]])
n <- length(country_iso_alpha3)
k <- unclass(windows[["country_iso_alpha3"]])

## Country code, ISO 3166-1 alpha-3
country_iso_alpha2 <- countrycode::countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "iso2c"
)
windows[["country_iso_alpha2"]]  <- gl(n, 1L, labels = country_iso_alpha2)[k]

## Country code, ISO 3166-1 numeric
country_iso_numeric <- countrycode::countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "iso3n"
)
windows[["country_iso_numeric"]] <- gl(n, 1L, labels = country_iso_numeric)[k]

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
windows[["country"]] <- gl(n, 1L, labels = country)[k]

## Subregion
subregion <- countrycode::countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = c(TWN = "Eastern Asia")
)
windows[["subregion"]] <- factor(subregion)[k]

## Region
region <- countrycode::countrycode(
  sourcevar = country_iso_alpha3,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = c(TWN = "Asia")
)
windows[["region"]] <- factor(region)[k]

## City population-weighted latitude and longitude
coords[["country_iso_alpha3"]] <- factor(coords[["country_iso_alpha3"]], levels = country_iso_alpha3)
f <- function(d) {
  c(weighted.mean(d[["longitude"]], d[["population"]]),
    weighted.mean(d[["latitude"]],  d[["population"]]))
}
longitude_latitude_rows <- c(by(coords, coords[["country_iso_alpha3"]], f, simplify = FALSE))
longitude_latitude <- matrix(unlist(longitude_latitude_rows), ncol = 2L, byrow = TRUE)
windows[c("longitude", "latitude")] <- longitude_latitude[k, , drop = FALSE]

## Longitude and latitude bands
longitude_band_20deg <- cut(
  x = longitude_latitude[, 1L],
  breaks = seq.int(0, 360, by = 20),
  right = FALSE,
  include.lowest = FALSE,
  ordered_result = TRUE
)
latitude_band_20deg <- cut(
  x = longitude_latitude[, 2L],
  breaks = seq.int(-90, 90, by = 20),
  right = FALSE,
  include.lowest = TRUE,
  ordered_result = TRUE
)
windows[["longitude_band_20deg"]] <- longitude_band_20deg[k]
windows[["latitude_band_20deg"]] <- latitude_band_20deg[k]

## Population size
population <- population[population[["year"]] == "2020", , drop = FALSE]
m <- match(country_iso_numeric, population[["country_iso_numeric"]], 0L)
windows[["population"]] <- population[["population"]][m][k]

## Mobility
s <- grep("^mobility_", names(mobility), value = TRUE)
windows[s] <- utils$get_summary(
  data = mobility[s],
  date = mobility[["date"]],
  start = windows[["start"]] - 21,
  end = windows[["start"]] - 7,
  f1 = mobility[["country_iso_alpha2"]],
  f2 = windows[["country_iso_alpha2"]],
  presummarize = function(x) utils$linear(x, x0 = 1, x1 = NULL, period = 7L, geom = TRUE),
  summarize = function(x) utils$Mean(x, na.rm = TRUE, geom = TRUE, zero.rm = FALSE)
)[["X"]]

## NPI
s <- grep("^npi_", names(npi), value = TRUE)
s_ordered <- grep("^npi_(flag|indic_(?!(E3|E4|H4|H5))).*$", names(npi), value = TRUE, perl = TRUE)
s_numeric <- setdiff(s, s_ordered)
windows[s_ordered] <- utils$get_summary(
  data = npi[s_ordered],
  date = npi[["date"]],
  start = windows[["start"]] - 21,
  end = windows[["start"]] - 7,
  f1 = npi[["country_iso_alpha3"]],
  f2 = windows[["country_iso_alpha3"]],
  presummarize = function(x) utils$locf(x, x0 = 0),
  summarize = function(x) min(utils$Mode(x, na.rm = TRUE))
)[["X"]]
windows[s_numeric] <- utils$get_summary(
  data = npi[s_numeric],
  date = npi[["date"]],
  start = windows[["start"]] - 21,
  end = windows[["start"]] - 7,
  f1 = npi[["country_iso_alpha3"]],
  f2 = windows[["country_iso_alpha3"]],
  presummarize = function(x) utils$linear(x, x0 = NULL, x1 = NULL, period = 1L, geom = FALSE),
  summarize = function(x) utils$Mean(x, na.rm = TRUE, geom = FALSE)
)[["X"]]

## Vaccination
s <- grep("^vaccinated_", names(vaccine), value = TRUE)
windows[s] <- utils$get_summary(
  data = vaccine[s],
  date = vaccine[["date"]],
  start = windows[["start"]] - 21,
  end = windows[["start"]] - 7,
  f1 = vaccine[["country_iso_alpha3"]],
  f2 = windows[["country_iso_alpha3"]],
  presummarize = function(x) utils$linear(x, x0 = 0, x1 = NULL, period = 1L, geom = FALSE),
  summarize = function(x) utils$Mean(x, na.rm = TRUE, geom = FALSE)
)[["X"]]

## Economic indicators
map <- c(
  econ_gdp_per_capita = "GDP per capita (constant 2015 US$)",
  econ_gini           = "Gini index (World Bank estimate)"
)
devel[["country_iso_alpha3"]] <- factor(devel[["country_iso_alpha3"]], levels = country_iso_alpha3)
devel[["indic"]] <- factor(devel[["indic"]], levels = map, labels = names(map))
devel[["year"]] <- factor(devel[["year"]], levels = 2001:2020)
i <- complete.cases(devel[c("country_iso_alpha3", "indic", "year")])
devel <- devel[i, , drop = FALSE]
lodevel <- aggregate(devel["value"], by = devel[c("country_iso_alpha3", "indic")], utils$lo, value = TRUE, drop = FALSE)
windows[levels(lodevel[["indic"]])] <- c(tapply(lodevel[["value"]], lodevel[["indic"]], `[`, k, simplify = FALSE))
windows[["econ_gini"]] <- 0.01 * windows[["econ_gini"]]

## Weather
s <- grep("^weather_", names(weather), value = TRUE)
windows[s] <- utils$get_summary(
  data = weather[s],
  date = weather[["date"]],
  start = windows[["start"]] - 21,
  end = windows[["start"]] - 7,
  f1 = weather[["country_iso_alpha3"]],
  f2 = windows[["country_iso_alpha3"]],
  presummarize = function(x) utils$linear(x, x0 = NULL, x1 = NULL, period = 1L, geom = FALSE),
  summarize = function(x) utils$Mean(x, na.rm = TRUE, geom = FALSE)
)[["X"]]

## Number of days since:
## * 50 persons reported infected
## * 2 in 1 million persons reported infected
covid[["country_iso_alpha3"]] <- factor(covid[["country_iso_alpha3"]], levels = country_iso_alpha3)
f <- function(d, min = 0L, scale = NULL) {
  if (!is.null(scale)) {
    min <- min * scale
  }
  x <- d[["cases_total"]]
  n <- length(x)
  i <- 1L
  while (i <= n && x[i] < min) {
    i <- i + 1L
  }
  if (i == n + 1L) .Date(NA_real_) else d[i, "date"]
}
scale <- windows[match(country_iso_alpha3, windows[["country_iso_alpha3"]]), "population"]
date_50 <- .Date(c(by(covid, covid[["country_iso_alpha3"]], f, min = 50L)))
date_2in1m <- .Date(mapply(f, d = split(covid, covid[["country_iso_alpha3"]]), min = 2e-06, scale = scale))
windows[["days_since_50"]] <- as.numeric(windows[["start"]] - date_50[k])
windows[["days_since_2in1m"]] <- as.numeric(windows[["start"]] - date_2in1m[k])

saveRDS(windows, file = path_rds)
