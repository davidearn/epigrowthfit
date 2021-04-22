paths <- list.files("../weather", pattern = "[0-9]{8}\\.RData$", recursive = TRUE)
r <- range(as.Date(sub("^.*([0-9]{8})\\.RData$", "\\1", paths), format = "%Y%m%d"))

varnames <- c(
  # me = "metoffice"
  temperature         = "t1o5m",
  specific_humidity   = "sh",
  shortwave_radiation = "sw",
  precipitation       = "precip",
  wind_speed          = "windspeed"
)
dirnames <- vapply(varnames, sprintf, "", fmt = "../weather/%s_mean")

f <- function(dirname, Dates_full = seq(r[1L], r[2L], by = 1)) {
  paths <- list.files(dirname, pattern = "[0-9]{8}\\.RData", full.names = TRUE)
  if (length(paths) == 0L) {
    return(NULL)
  }
  Dates <- as.Date(sub("^.*([0-9]{8})\\.RData$", "\\1", paths), format = "%Y%m%d")
  lx <- lapply(paths, function(path) {load(path); x})
  x <- unlist(lx)
  d <- data.frame(
    country_iso_alpha3 = factor(names(x)),
    Date = rep.int(Dates, lengths(lx)),
    value = x
  )
  m <- match(Dates_full, Dates, 0L) == 0L
  if (any(m)) {
    M <- sum(m)
    cia3 <- levels(d$country_iso_alpha3)
    N <- length(cia3)
    d_fill <- data.frame(
      country_iso_alpha3 = gl(N, M, labels = cia3),
      Date = rep.int(Dates_full[m], N),
      value = NA_real_
    )
    d <- rbind(d, d_fill)
  }
  o <- do.call(order, d[1:2])
  d <- d[o, , drop = FALSE]
  row.names(d) <- NULL
  d
}

ld <- lapply(dirnames, f)
stopifnot(length(unique(lapply(ld, `[`, 1:2))) == 1L)
weather <- data.frame(
  ld[[1L]][1:2],
  lapply(ld, `[[`, 3L)
)
names(weather)[-(1:2)] <- sprintf("weather_%s", names(weather)[-(1:2)])
save(weather, file = "../weather.RData")
