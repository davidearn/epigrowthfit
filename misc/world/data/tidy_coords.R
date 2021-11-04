(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  map <- c(
    country_iso_alpha3 = "iso3",
    city               = "city",
    population         = "population",
    latitude           = "lat",
    longitude          = "lng"
  )
  i <- complete.cases(d0[map])
  d1 <- d0[i, map, drop = FALSE]
  names(d1) <- names(map)

  m <- match("Jerusalem", d1[["city"]], 0L)
  r <- d1[m, , drop = FALSE]
  r[["country_iso_alpha3"]] <- "PSE"
  d2 <- rbind(d1, r)

  d2[["country_iso_alpha3"]] <- factor(d2[["country_iso_alpha3"]])
  d2[["population"]] <- trunc(d2[["population"]])
  d2[["longitude"]] <- d2[["longitude"]] %% 360

  o <- do.call(order, unname(d2[c("country_iso_alpha3", "city", "population")]))
  d3 <- d2[o, , drop = FALSE]
  row.names(d3) <- NULL
  d3
}

coords_raw <- read.csv(path_csv)
coords <- tidy(coords_raw)
saveRDS(coords, file = path_rds)
str(coords)
