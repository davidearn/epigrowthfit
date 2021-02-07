## Retrieve full DELVE COVID-19 data set
loc <- paste0(
  "https://raw.githubusercontent.com/",
  "rs-delve/covid19_datasets/master/dataset/",
  "combined_dataset_latest.csv"
)
world_full <- read.csv(url(loc))

## Extract incidence
world <- world_full[c("DATE", "cases_new", "ISO", "country_name")]

## Make nice
names(world) <- c("date", "cases", "country_code", "country")
world$date <- as.Date(world$date)
world$cases <- as.integer(world$cases)
world$country_code <- factor(world$country_code)
world$country <- factor(world$country)

## Prefer familiar country names
lev <- levels(world$country)
old <- c(
  "Aruba",
  "Bolivia, Plurinational State of",
  "Brunei Darussalam",
  "Congo, The Democratic Republic of the",
  "Czechia",
  "Faroe Islands",
  "Greenland",
  "Guam",
  "Iran, Islamic Republic of",
  "Korea, Republic of",
  "Lao People's Democratic Republic",
  "Moldova, Republic of",
  "Palestine, State of",
  "Puerto Rico",
  "Republic of Kosovo",
  "Russian Federation",
  "Syrian Arab Republic",
  "Taiwan, Province of China",
  "Tanzania, United Republic of",
  "Venezuela, Bolivarian Republic of",
  "Viet Nam",
  "Virgin Islands, U.S."
)
new <- c(
  "Aruba (Netherlands)",
  "Bolivia",
  "Brunei",
  "D. R. Congo",
  "Czech Republic",
  "Faroe Islands (Denmark)",
  "Greenland (Denmark)",
  "Guam (United States)",
  "Iran",
  "South Korea",
  "Laos",
  "Moldova",
  "Palestine",
  "Puerto Rico (United States)",
  "Kosovo",
  "Russia",
  "Syria",
  "Taiwan",
  "Tanzania",
  "Venezuela",
  "Vietnam",
  "Virgin Islands (United States)"
)
lev[match(old, lev)] <- new
levels(world$country) <- lev

## Group countries by continent
library("countrycode")
world$continent <- factor(countrycode(
  sourcevar = world$country_code,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = `names<-`(c("Europe", "Asia"), c("RKS", "TWN"))
))
world$region <- factor(countrycode(
  sourcevar = world$country_code,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = `names<-`(c("Southern Europe", "Eastern Asia"), c("RKS", "TWN"))
))

## Replace negative observations with NA
world$cases[world$cases < 0L] <- NA

## Omit leading NA from time series
world_split <- lapply(split(world, world$country_code), function(d) {
  i <- which(!is.na(d$cases))
  if (length(i) > 0L) {
    d[seq.int(i[1L], nrow(d)), , drop = FALSE]
  } else {
    NULL
  }
})
world <- droplevels(do.call(rbind, world_split))
row.names(world) <- NULL

## Save
save(world, file = "world.RData")
