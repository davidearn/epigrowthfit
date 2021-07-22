library("countrycode")

covid19 <- read.csv("csv/covid19.csv")
covid19$country_iso_alpha3 <- countrycode(
  sourcevar = covid19$Country.Region,
  origin = "country.name",
  destination = "iso3c"
)
s <- grep("^X[0-9]+\\.[0-9]+\\.[0-9]+$", names(covid19), value = TRUE)
f <- function(d) {
  if (nrow(d) < 2L) {
    return(d)
  }
  blank <- !nzchar(d$Province.State)
  if (sum(blank) == 1L) {
    return(d[blank, , drop = FALSE])
  }
  d$Province.State <- ""
  r <- d[1L, , drop = FALSE]
  r[s] <- lapply(d[s], sum, na.rm = FALSE)
  r
}
covid19 <- do.call(rbind, by(covid19, covid19$country_iso_alpha3, f, simplify = FALSE))
covid19 <- data.frame(
  country_iso_alpha3 = rep(factor(covid19$country_iso_alpha3), times = length(s)),
  Date = rep(as.Date(s, format = "X%m.%d.%y"), each = nrow(covid19)),
  cases_total = unlist(covid19[s], FALSE, FALSE),
  cases_new = NA_integer_
)
o <- do.call(order, unname(covid19[c("country_iso_alpha3", "Date")]))
covid19 <- covid19[o, , drop = FALSE]
split(covid19$cases_new, covid19$country_iso_alpha3) <-
  c(tapply(covid19$cases_total, covid19$country_iso_alpha3, function(x) c(NA, diff(x)), simplify = FALSE))
row.names(covid19) <- NULL
saveRDS(covid19, file = "rds/covid19.rds")
rm(covid19)

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
saveRDS(coords, file = "rds/coords.rds")
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
saveRDS(population, file = "rds/population.rds")
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
saveRDS(mobility, file = "rds/mobility.rds")
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
npi[names(s_index)] <- npi[names(s_index)] / 100
o <- do.call(order, unname(npi[c("country_iso_alpha3", "Date")]))
npi <- npi[o, , drop = FALSE]
row.names(npi) <- NULL
saveRDS(npi, file = "rds/npi.rds")
rm(npi)

vaccine <- read.csv("csv/vaccine.csv")
s <- c(
  country_iso_alpha3       = "iso_code",
  Date                     = "date",
  vaccinated_geq1d_per_100 = "people_vaccinated_per_hundred",
  vaccinated_fully_per_100 = "people_fully_vaccinated_per_hundred"
)
i <- grep("^OWID_", vaccine$iso_code, invert = TRUE)
vaccine <- vaccine[i, s, drop = FALSE]
names(vaccine) <- names(s)
vaccine$country_iso_alpha3 <- factor(vaccine$country_iso_alpha3)
vaccine$Date <- as.Date(vaccine$Date)
o <- do.call(order, unname(vaccine[c("country_iso_alpha3", "Date")]))
vaccine <- vaccine[o, , drop = FALSE]
row.names(vaccine) <- NULL
saveRDS(vaccine, file = "rds/vaccine.rds")
rm(vaccine)

devel <- read.csv("csv/devel.csv")
s <- grep("^X[0-9]{4}$", names(devel), value = TRUE)
devel <- data.frame(
  country_iso_alpha3 = rep(factor(devel$Country.Code), times = length(s)),
  indic = rep(factor(devel$Indicator.Name), times = length(s)),
  indic_code = rep(factor(devel$Indicator.Code), times = length(s)),
  year = rep(ordered(substr(s, 2L, 5L)), each = nrow(devel)),
  value = unlist(devel[s], FALSE, FALSE)
)
m <- match(levels(devel$indic), devel$indic)
devel$indic_code <- factor(devel$indic_code, levels = devel$indic_code[m])
o <- do.call(order, unname(devel[c("country_iso_alpha3", "indic", "year")]))
devel <- devel[o, , drop = FALSE]
row.names(devel) <- NULL
saveRDS(devel, file = "rds/devel.rds")
rm(devel)

equity <- read.csv("csv/equity.csv")
s <- grep("^X[0-9]{4}$", names(equity), value = TRUE)
equity <- data.frame(
  country_iso_alpha3 = rep(factor(equity$Country.Code), times = length(s)),
  indic = rep(factor(equity$Indicator.Name), times = length(s)),
  indic_code = rep(factor(equity$Indicator.Code), times = length(s)),
  year = rep(ordered(substr(s, 2L, 5L)), each = nrow(equity)),
  value = unlist(equity[s], FALSE, FALSE)
)
m <- match(levels(equity$indic), equity$indic)
equity$indic_code <- factor(equity$indic_code, levels = equity$indic_code[m])
o <- do.call(order, unname(equity[c("country_iso_alpha3", "indic", "year")]))
equity <- equity[o, , drop = FALSE]
row.names(equity) <- NULL
saveRDS(equity, file = "rds/equity.rds")
rm(equity)
