(path_csv <- Sys.getenv("PATH_CSV"))
(path_rds <- Sys.getenv("PATH_RDS"))

tidy <- function(d0) {
  map <- c(
    country_iso_alpha3 = "Country.Code",
    indic              = "Indicator.Name",
    indic_code         = "Indicator.Code"
  )
  d0[names(map)] <- lapply(d0[map], factor)
  m <- match(seq_len(nlevels(d0[["indic"]])), unclass(d0[["indic"]]), 0L)
  d0[["indic_code"]] <- factor(d0[["indic_code"]], levels = d0[["indic_code"]][m])

  x <- grep("^X[0-9]{4}", names(d0), value = TRUE)
  d1 <- reshape(
    data = d0,
    direction = "long",
    varying = x,
    v.names = "value",
    timevar = "year",
    idvar = names(map),
    drop = setdiff(names(d0), c(names(map), x)),
    times = ordered(sub("^X", "", x)),
  )
  attr(d1, "reshapeLong") <- NULL

  o <- do.call(order, unname(d1[c("country_iso_alpha3", "indic", "year")]))
  d2 <- d1[o, , drop = FALSE]
  row.names(d2) <- NULL
  d2
}

equity_raw <- read.csv(path_csv)
equity <- tidy(equity_raw)
saveRDS(equity, file = path_rds)
str(equity)
