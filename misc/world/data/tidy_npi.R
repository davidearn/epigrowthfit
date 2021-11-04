(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  map_indic <- grep("^(C|E|H)[0-9]+_(?!Flag)", names(d0), value = TRUE, perl = TRUE)
  names(map_indic) <- sub("^((C|E|H)[0-9]+)_.*", "npi_indic_\\1", map_indic)
  map_flag  <- grep("^(C|E|H)[0-9]+_Flag$", names(d0), value = TRUE)
  names(map_flag) <- sub("^((C|E|H)[0-9]+)_.*", "npi_flag_\\1", map_flag)
  map_index <- grep("^(?!.*Legacy).*Index$", names(d0), value = TRUE, perl = TRUE)
  names(map_index) <- sub("(.*)_Index", "npi_index_\\L\\1", gsub("([A-Z])(?<!^.)", "_\\1", map_index, perl = TRUE), perl = TRUE)
  map <- c(
    country_iso_alpha3 = "CountryCode",
    date               = "Date",
    map_indic,
    map_flag,
    map_index
  )
  i <- !nzchar(d0[["RegionCode"]])
  d1 <- d0[i, map, drop = FALSE]
  names(d1) <- names(map)

  d1[["country_iso_alpha3"]] <- factor(d1[["country_iso_alpha3"]])
  d1[["date"]] <- as.Date(sprintf("%d", d1[["date"]]), format = "%Y%m%d")
  s <- c(grep("(?<!(E3|E4|H4|H5))$", names(map_indic), value = TRUE, perl = TRUE), names(map_flag))
  d1[s] <- lapply(d1[s], ordered)
  s <- names(map_index)
  d1[s] <- 0.01 * d1[s]

  o <- do.call(order, unname(d1[c("country_iso_alpha3", "date")]))
  d2 <- d1[o, , drop = FALSE]
  row.names(d2) <- NULL
  d2
}

npi_raw <- read.csv(path_csv)
npi <- tidy(npi_raw)
saveRDS(npi, file = path_rds)
str(npi)
