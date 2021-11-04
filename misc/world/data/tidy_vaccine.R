(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  map <- c(
    country_iso_alpha3       = "iso_code",
    date                     = "date",
    vaccinated_geq1d_per_100 = "people_vaccinated_per_hundred",
    vaccinated_fully_per_100 = "people_fully_vaccinated_per_hundred"
  )
  i <- grep("^OWID_", d0[["iso_code"]], invert = TRUE)
  d1 <- d0[i, map, drop = FALSE]
  names(d1) <- names(map)

  d1[["country_iso_alpha3"]] <- factor(d1[["country_iso_alpha3"]])
  d1[["date"]] <- as.Date(d1[["date"]])

  o <- do.call(order, unname(d1[c("country_iso_alpha3", "date")]))
  d2 <- d1[o, , drop = FALSE]
  row.names(d2) <- NULL
  d2
}

vaccine_raw <- read.csv(path_csv)
vaccine <- tidy(vaccine_raw)
saveRDS(vaccine, file = path_rds)
str(vaccine)
