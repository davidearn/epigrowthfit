## Retrieve full DELVE COVID-19 data set
loc <- paste0(
  "https://raw.githubusercontent.com/",
  "rs-delve/covid19_datasets/master/dataset/",
  "combined_dataset_latest.csv"
)
world_full <- read.csv(url(loc))
world <- world_full[c("DATE", "cases_new", "ISO", "country_name")]
names(world) <- c("date", "cases", "country_code", "country")

## Clean
world$date <- as.Date(world$date)
world$cases[world$cases < 0] <- NA
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
lev[match(old, lev, 0L)] <- new
levels(world$country) <- lev

## Group countries by region and continent
library("countrycode")
world$region <- factor(countrycode(
  sourcevar = world$country_code,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = `names<-`(c("Southern Europe", "Eastern Asia"), c("RKS", "TWN"))
))
world$continent <- factor(countrycode(
  sourcevar = world$country_code,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = `names<-`(c("Europe", "Asia"), c("RKS", "TWN"))
))

## Clean again
f <- function(x) {
  if (length(i <- which(!is.na(x))) == 0L) {
    return(i)
  }
  seq.int(i[1L], length(x))
}
keep <- tapply(world$cases, world$country, f, simplify = FALSE)
world_split <- split(world, world$country, drop = TRUE)
world_split <- Map(`[`, world_split, keep, list(seq_along(world)), drop = FALSE)
world <- droplevels(do.call(rbind, world_split))
row.names(world) <- NULL

## Save
save(world, file = "world.RData")
