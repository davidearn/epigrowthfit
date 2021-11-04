(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  map <- c(
    country_iso_numeric = "LocID",
    year                = "Time",
    population          = "PopTotal"
  )
  i <- d0[["Variant"]] == "Medium"
  d1 <- d0[i, map, drop = FALSE]
  names(d1) <- names(map)

  d1[["year"]] <- ordered(d1[["year"]])
  d1[["population"]] <- trunc(1000 * d1[["population"]])

  o <- do.call(order, unname(d1[c("country_iso_numeric", "year")]))
  d2 <- d1[o, , drop = FALSE]
  row.names(d2) <- NULL
  d2
}

population_raw <- read.csv(path_csv)
population <- tidy(population_raw)
saveRDS(population, file = path_rds)
str(population)
