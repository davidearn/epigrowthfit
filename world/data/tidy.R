library("countrycode")

covid <- read.csv("covid.csv")
covid$country_iso_alpha3 <- countrycode(
  sourcevar = covid$Country.Region,
  origin = "country.name",
  destination = "iso3c"
)
s <- grep("^X[0-9]+\\.[0-9]+\\.[0-9]+$", names(covid), value = TRUE)
f <- function(d) {
  if (nrow(d) < 2L) {
    return(d)
  }
  blank <- d$Province.State == ""
  if (sum(blank) == 1L) {
    d[blank, , drop = FALSE]
  }
  d$Province.State <- ""
  r <- d[1L, , drop = FALSE]
  r[s] <- lapply(d[s], sum, na.rm = FALSE)
  r
}
covid <- do.call(rbind, lapply(split(covid, covid$country_iso_alpha3), f))
covid <- data.frame(
  country_iso_alpha3 = rep(factor(covid$country_iso_alpha3), times = length(s)),
  Date = rep(as.Date(s, format = "X%m.%d.%y"), each = nrow(covid)),
  total_cases = unlist(covid[s], recursive = FALSE, use.names = FALSE),
  cases = NA_integer_
)
o <- do.call(order, covid[c("country_iso_alpha3", "Date")])
covid <- covid[o, , drop = FALSE]
split(covid$cases, covid$country_iso_alpha3) <-
  tapply(covid$total_cases, covid$country_iso_alpha3, function(x) c(NA, diff(x)), simplify = FALSE)
covid$total_cases <- NULL
row.names(covid) <- NULL
save(covid, file = "covid.RData")

coords <- read.csv("coords.csv")
s <- c("iso3", "city", "population", "lat", "lng")
i <- complete.cases(coords[s])
coords <- coords[i, s, drop = FALSE]
names(coords) <- c("country_iso_alpha3", "city", "population", "latitude", "longitude")
m <- match("Jerusalem", coords$city, 0L)
coords <- rbind(coords, coords[m, , drop = FALSE])
coords$country_iso_alpha3[nrow(coords)] <- "PSE"
coords$country_iso_alpha3 <- factor(coords$country_iso_alpha3)
o <- do.call(order, coords[c("country_iso_alpha3", "population")])
coords <- coords[o, , drop = FALSE]
row.names(coords) <- NULL
save(coords, file = "coords.RData")

population <- read.csv("population.csv")
s <- c("LocID", "PopTotal")
i <- population$Time == 2020L & population$Variant == "Medium"
population <- population[i, s, drop = FALSE]
names(population) <- c("country_iso_numeric", "population")
population$population <- 1000 * population$population
o <- order(population$country_iso_numeric)
population <- population[o, , drop = FALSE]
row.names(population) <- NULL
save(population, file = "population.RData")

mobility <- read.csv("mobility.csv")
varnames <- grep("_percent_change_from_baseline$", names(mobility), value = TRUE)
s <- c("country_region_code", "date", varnames)
i <-
  !is.na(mobility$country_region_code) &
  mobility$sub_region_1 == "" &
  mobility$sub_region_2 == "" &
  mobility$metro_area   == ""
mobility <- mobility[i, s, drop = FALSE]
mobility[varnames] <- 1 + mobility[varnames] / 100
varnames <- sub("_percent_change_from_baseline$", "", varnames)
names(mobility) <- c("country_iso_alpha2", "Date", sprintf("mobility_%s", varnames))
mobility$country_iso_alpha2 <- factor(mobility$country_iso_alpha2)
mobility$Date <- as.Date(mobility$Date)
o <- do.call(order, mobility[c("country_iso_alpha2", "Date")])
mobility <- mobility[o, , drop = FALSE]
row.names(mobility) <- NULL
save(mobility, file = "mobility.RData")

npi <- read.csv("npi.csv")
s_indic <- grep("^(C|E|H)[0-9]+_(?!Flag)", names(npi), value = TRUE, perl = TRUE)
s_flag  <- grep("^(C|E|H)[0-9]+_Flag$", names(npi), value = TRUE)
s_index <- grep("^(?!.*Legacy).*Index$", names(npi), value = TRUE, perl = TRUE)
s <- c("CountryCode", "Date", s_indic, s_flag, s_index)
i <- npi$RegionCode == ""
npi <- npi[i, s, drop = FALSE]
s_indic <- sub("^((C|E|H)[0-9]+)_.*", "npi_indic_\\1", s_indic)
s_flag <- sub("^((C|E|H)[0-9]+)_.*", "npi_flag_\\1", s_flag)
s_index <- sub("(.*)_Index", "npi_index_\\L\\1", gsub("([A-Z])(?<!^.)", "_\\1", s_index, perl = TRUE), perl = TRUE)
names(npi) <- c("country_iso_alpha3", "Date", s_indic, s_flag, s_index)
npi$country_iso_alpha3 <- factor(npi$country_iso_alpha3)
npi$Date <- as.Date(as.character(npi$Date), format = "%Y%m%d")
s_ordered <- c(grep("(?<!(E3|E4|H4|H5))$", s_indic, value = TRUE, perl = TRUE), s_flag)
npi[s_ordered] <- lapply(npi[s_ordered], ordered)
o <- do.call(order, npi[c("country_iso_alpha3", "Date")])
npi <- npi[o, , drop = FALSE]
row.names(npi) <- NULL
save(npi, file = "npi.RData")
