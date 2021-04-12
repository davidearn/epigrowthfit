## Retrieve full DELVE COVID-19 data set
loc <- paste0(
  "https://raw.githubusercontent.com/",
  "rs-delve/covid19_datasets/master/dataset/",
  "combined_dataset_latest.csv"
)
world_full <- read.csv(url(loc))
world <- world_full[c("ISO", "DATE", "cases_new")]
names(world) <- c("country_iso3", "date", "cases")

## Use friendly classes
world$country_iso3 <- factor(world$country_iso3)
world$date <- as.Date(world$date)

## Treat negative numbers as missing
world$cases[world$cases < 0] <- NA

## Discard time series with fewer than 1000 cases
tot <- tapply(world$cases, world$country_iso3, sum, na.rm = TRUE)
tot <- tot[tot >= 1000]
world$country_iso3 <- factor(world$country_iso3, levels = names(tot))
world <- world[!is.na(world$country_iso3), , drop = FALSE]

## Discard leading NA
f <- function(x) {
  n <- length(x)
  replace(rep_len(FALSE, n), seq.int(which.min(is.na(x)), n), TRUE)
}
keep <- unlist(tapply(world$cases, world$country_iso3, f, simplify = FALSE))
world <- world[keep, , drop = FALSE]
world_split <- split(world, world$country_iso3)

## Aggregate weekly
g <- function(d) {
  n <- nrow(d)
  n <- n - n %% 7L
  d7 <- d[seq.int(7L, n, by = 7L), , drop = FALSE]
  d7$cases <- tapply(d$cases[seq_len(n)], gl(n / 7L, 7L), sum)
  r0 <- d[1L, , drop = FALSE]
  r0$date <- r0$date - 1
  r0$cases <- NA_real_
  rbind(r0, d7)
}
world7_split <- lapply(world_split, g)

## Delete spurious zeros
p0 <- tapply(world$cases, world$country_iso3,
             function(x) sum(x == 0, na.rm = TRUE) / sum(!is.na(x)))
h_ <- function(x, b, tol) {
  zero <- !is.na(x) & x == 0
  if (any(zero)) {
    if (b < 3 || b %% 2 != 1) {
      stop("`b` must be an odd number greater than 1.")
    }
    p <- rep_len(NA_real_, (b - 1) / 2)
    X <- embed(c(p, x, p), b)
    zero[zero] <- apply(X[zero, , drop = FALSE] > tol, 1L, any, na.rm = TRUE)
  }
  !zero
}
h <- function(d, b, tol) {
  ok <- h_(d$cases, b, tol)
  d[ok, , drop = FALSE]
}
world_split   <- lapply(world_split,  h, b = 15, tol = 15)
world7_split  <- lapply(world7_split, h, b = 3,  tol = 90)
world  <- do.call(rbind, world_split)
world7 <- do.call(rbind, world7_split)
row.names(world)  <- NULL
row.names(world7) <- NULL

## Retrieve country name, region, continent, and population
## from ISO codes
library("countrycode")
country_iso3 <- levels(world$country_iso3)
country_name <- countrycode(
  sourcevar = country_iso3,
  origin = "iso3c",
  destination = "country.name",
  custom_match = c(MMR = "Myanmar",
                   PSE = "Palestine",
                   TTO = "Trinidad and Tobago")
)
region <- countrycode(
  sourcevar = country_iso3,
  origin = "iso3c",
  destination = "un.regionsub.name",
  custom_match = c(TWN = "Eastern Asia")
)
continent <- countrycode(
  sourcevar = country_iso3,
  origin = "iso3c",
  destination = "un.region.name",
  custom_match = c(TWN = "Asia")
)
population <- tapply(world_full$stats_population, world_full$ISO, `[`, 1L)
population <- population[country_iso3]

## Compute population-weighted average latitude and longitude of cities
load("cities/cities.RData")
do_wmll <- function(d) {
  c(latitude  = weighted.mean(d$lat, d$population),
    longitude = weighted.mean(d$lng, d$population))
}
ll_list <- by(cities, cities$iso3, do_wmll, simplify = FALSE)
m <- match("Jerusalem", cities$city, 0L)
ll_list$PSE <- c(latitude = cities$lat[m], longitude = cities$lng[m])
ll_mat <- matrix(unlist(ll_list), ncol = 2L, byrow = TRUE,
             dimnames = list(names(ll_list), names(ll_list[[1L]])))

## Keep all time invariant variables in a separate data frame
worldstats <- data.frame(
  country_iso3,
  country_name,
  region,
  continent,
  population,
  ll_mat[country_iso3, , drop = FALSE],
  row.names = NULL,
  stringsAsFactors = TRUE
)

## Save everything
tot <- sort(tot)
p0 <- sort(p0, decreasing = TRUE)
save(world, world7, worldstats, tot, p0, file = "world.RData")

