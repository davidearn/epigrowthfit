library("countrycode")

covid <- read.csv("csv/covid.csv")
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
    return(d[blank, , drop = FALSE])
  }
  d$Province.State <- ""
  r <- d[1L, , drop = FALSE]
  r[s] <- lapply(d[s], sum, na.rm = FALSE)
  r
}
covid <- do.call(rbind, by(covid, covid$country_iso_alpha3, f, simplify = FALSE))
covid <- data.frame(
  country_iso_alpha3 = rep(factor(covid$country_iso_alpha3), times = length(s)),
  Date = rep(as.Date(s, format = "X%m.%d.%y"), each = nrow(covid)),
  cases_total = unlist(covid[s], recursive = FALSE, use.names = FALSE),
  cases_new = NA_integer_
)
o <- do.call(order, unname(covid[c("country_iso_alpha3", "Date")]))
covid <- covid[o, , drop = FALSE]
split(covid$cases_new, covid$country_iso_alpha3) <-
  tapply(covid$cases_total, covid$country_iso_alpha3, function(x) c(NA, diff(x)), simplify = FALSE)
row.names(covid) <- NULL
save(covid, file = "covid.RData")
rm(covid)

coords <- read.csv("csv/coords.csv")
s <- c(
  country_iso_alpha3 = "iso3",
  city               = "city",
  population         = "population",
  latitude           = "lat",
  longitude          = "lng"
)
i <- complete.cases(coords[s])
coords <- coords[i, s, drop = FALSE]
names(coords) <- names(s)
coords$longitude <- coords$longitude %% 360
m <- match("Jerusalem", coords$city, 0L)
coords <- rbind(coords, coords[m, , drop = FALSE])
coords$country_iso_alpha3[nrow(coords)] <- "PSE"
coords$country_iso_alpha3 <- factor(coords$country_iso_alpha3)
o <- do.call(order, unname(coords[c("country_iso_alpha3", "city", "population")]))
coords <- coords[o, , drop = FALSE]
row.names(coords) <- NULL
save(coords, file = "coords.RData")
rm(coords)

population <- read.csv("csv/population.csv")
s <- c(
  country_iso_numeric = "LocID",
  year                = "Time",
  population          = "PopTotal"
)
i <- population$Variant == "Medium"
population <- population[i, s, drop = FALSE]
names(population) <- names(s)
population$year <- ordered(population$year)
population$population <- 1000 * population$population
o <- do.call(order, unname(population[c("country_iso_numeric", "year")]))
population <- population[o, , drop = FALSE]
row.names(population) <- NULL
save(population, file = "population.RData")
rm(population)

mobility <- read.csv("csv/mobility.csv")
varnames <- grep("_percent_change_from_baseline$", names(mobility), value = TRUE)
names(varnames) <- sub("^(.*)_percent_change_from_baseline$", "mobility_\\1", varnames)
s <- c(
  country_iso_alpha2 = "country_region_code",
  Date               = "date",
  varnames
)
i <-
  !is.na(mobility$country_region_code) &
  mobility$sub_region_1 == "" &
  mobility$sub_region_2 == "" &
  mobility$metro_area   == ""
mobility <- mobility[i, s, drop = FALSE]
mobility[varnames] <- 1 + mobility[varnames] / 100
names(mobility) <- names(s)
mobility$country_iso_alpha2 <- factor(mobility$country_iso_alpha2)
mobility$Date <- as.Date(mobility$Date)
o <- do.call(order, unname(mobility[c("country_iso_alpha2", "Date")]))
mobility <- mobility[o, , drop = FALSE]
row.names(mobility) <- NULL
save(mobility, file = "mobility.RData")
rm(mobility)

npi <- read.csv("csv/npi.csv")
s_indic <- grep("^(C|E|H)[0-9]+_(?!Flag)", names(npi), value = TRUE, perl = TRUE)
names(s_indic) <- sub("^((C|E|H)[0-9]+)_.*", "npi_indic_\\1", s_indic)
s_flag  <- grep("^(C|E|H)[0-9]+_Flag$", names(npi), value = TRUE)
names(s_flag) <- sub("^((C|E|H)[0-9]+)_.*", "npi_flag_\\1", s_flag)
s_index <- grep("^(?!.*Legacy).*Index$", names(npi), value = TRUE, perl = TRUE)
names(s_index) <- sub("(.*)_Index", "npi_index_\\L\\1", gsub("([A-Z])(?<!^.)", "_\\1", s_index, perl = TRUE), perl = TRUE)
s_ordered <- c(grep("(?<!(E3|E4|H4|H5))$", names(s_indic), value = TRUE, perl = TRUE), names(s_flag))
s <- c(
  country_iso_alpha3 = "CountryCode",
  Date               = "Date",
  s_indic,
  s_flag,
  s_index
)
i <- npi$RegionCode == ""
npi <- npi[i, s, drop = FALSE]
names(npi) <- names(s)
npi$country_iso_alpha3 <- factor(npi$country_iso_alpha3)
npi$Date <- as.Date(as.character(npi$Date), format = "%Y%m%d")
npi[s_ordered] <- lapply(npi[s_ordered], ordered)
o <- do.call(order, unname(npi[c("country_iso_alpha3", "Date")]))
npi <- npi[o, , drop = FALSE]
row.names(npi) <- NULL
save(npi, file = "npi.RData")
rm(npi)

devel <- read.csv("csv/devel.csv")
s <- grep("^X[0-9]{4}$", names(devel), value = TRUE)
devel <- data.frame(
  country_iso_alpha3 = rep(factor(devel$Country.Code), times = length(s)),
  indic = rep(factor(devel$Indicator.Name), times = length(s)),
  indic_code = rep(factor(devel$Indicator.Code), times = length(s)),
  year = rep(ordered(substr(s, 2L, 5L)), each = nrow(devel)),
  value = unlist(devel[s], recursive = FALSE, use.names = FALSE)
)
m <- match(levels(devel$indic), devel$indic)
devel$indic_code <- factor(devel$indic_code, levels = devel$indic_code[m])
o <- do.call(order, unname(devel[c("country_iso_alpha3", "indic", "year")]))
devel <- devel[o, , drop = FALSE]
row.names(devel) <- NULL
save(devel, file = "devel.RData")
rm(devel)

equity <- read.csv("csv/equity.csv")
s <- grep("^X[0-9]{4}$", names(equity), value = TRUE)
equity <- data.frame(
  country_iso_alpha3 = rep(factor(equity$Country.Code), times = length(s)),
  indic = rep(factor(equity$Indicator.Name), times = length(s)),
  indic_code = rep(factor(equity$Indicator.Code), times = length(s)),
  year = rep(factor(ordered(s, 2L, 5L)), each = nrow(equity)),
  value = unlist(equity[s], recursive = FALSE, use.names = FALSE)
)
m <- match(levels(equity$indic), equity$indic)
equity$indic_code <- factor(equity$indic_code, levels = equity$indic_code[m])
o <- do.call(order, unname(equity[c("country_iso_alpha3", "indic", "year")]))
equity <- equity[o, , drop = FALSE]
row.names(equity) <- NULL
save(equity, file = "equity.RData")
rm(equity)
