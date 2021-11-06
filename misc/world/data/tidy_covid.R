(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  stopifnot(requireNamespace("countrycode"))
  d0[["country_iso_alpha3"]] <- factor(countrycode::countrycode(
    sourcevar = d0[["Country.Region"]],
    origin = "country.name",
    destination = "iso3c"
  ))
  i <- grep("^X\\d+\\.\\d+\\.\\d+$", names(d0))
  d1 <- aggregate(d0[i], d0["country_iso_alpha3"], sum, na.rm = TRUE)

  d2 <- reshape(
    data = d1,
    direction = "long",
    varying = names(d1)[-1L],
    v.names = "cases_total",
    timevar = "date",
    idvar = names(d1)[1L],
    times = as.Date(names(d1)[-1L], format = "X%m.%d.%y")
  )
  attr(d2, "reshapeLong") <- NULL

  o <- do.call(order, unname(d2[c("country_iso_alpha3", "date")]))
  d3 <- d2[o, , drop = FALSE]
  row.names(d3) <- NULL
  d3
}

covid_raw <- read.csv(path_csv)
covid <- tidy(covid_raw)
saveRDS(covid, file = path_rds)
str(covid)
