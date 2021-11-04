(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  map_ratio <- grep("_percent_change_from_baseline$", names(d0), value = TRUE)
  names(map_ratio) <- sub("^(.*)_percent_change_from_baseline$", "mobility_\\1", map_ratio)
  map <- c(
    country_iso_alpha2 = "country_region_code",
    date               = "date",
    map_ratio
  )
  i <- !Reduce(`|`, lapply(d0[c("sub_region_1", "sub_region_2", "metro_area")], nzchar), is.na(d0[["country_region_code"]]))
  d1 <- d0[i, map, drop = FALSE]
  names(d1) <- names(map)

  d1[["country_iso_alpha2"]] <- factor(d1[["country_iso_alpha2"]])
  d1[["date"]] <- as.Date(d1[["date"]])
  s <- names(map_ratio)
  d1[s] <- 1 + 0.01 * d1[s]

  o <- do.call(order, unname(d1[c("country_iso_alpha2", "date")]))
  d2 <- d1[o, , drop = FALSE]
  row.names(d2) <- NULL
  d2
}

mobility_raw <- read.csv(path_csv)
mobility <- tidy(mobility_raw)
saveRDS(mobility, file = path_rds)
str(mobility)
