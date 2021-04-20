library("parallel")
library("ncdf4")
library("httr")
RNGkind("L'Ecuyer-CMRG")

## Data frame supplying latitude, longitude, population
## of major cities
load("../coords.RData")
coords <- coords[order(coords$country_iso_alpha3), , drop = FALSE] # just in case

## Parent directory for Met Office downloads
parent <- "https://metdatasa.blob.core.windows.net/covid19-response/metoffice_global_daily"

## Path to output
path_to_weather <- "../weather.RData"

## Variables on which country means are desired
varnames <- c(
# mikael                 metoffice    ncdf4
  "temperature",         "t1o5m",     "air_temperature",
  "specific_humidity",   "sh",        "specific_humidity",
  "shortwave_radiation", "sw",        "m01s01i202",
  "precipitation",       "precip",    "precipitation_flux",
  "wind_speed",          "windspeed", "wind_speed"
)
varnames <- matrix(varnames, ncol = 3L, byrow = TRUE,
                   dimnames = list(NULL, c("mikael", "metoffice", "ncdf4")))
varnames <- as.data.frame(varnames, stringsAsFactors = FALSE)

## Dates on which country means are desired
yesterday <- Sys.Date() - 1
Dates <- seq(as.Date("2020-01-01"), yesterday, by = 1)

## Initializing/augmenting output
if (file.exists(path_to_weather)) {
  load(path_to_weather) # weather X Y dX dY
  lastday <- max(weather$Date)
  append <- data.frame(
    country_iso_alpha3 = gl(nlevels(coords$country_iso_alpha3), yesterday - lastday, labels = levels(coords$country_iso_alpha3)),
    Date = rep.int(Dates[Dates > lastday], N)
  )
  s0 <- names(append)
  s1 <- names(weather)
  s2 <- c(s0, varnames$mikael)
  append[setdiff(s1, s0)] <- NA_real_
  weather <- rbind(weather, append[s1])
  weather[setdiff(s2, s1)] <- NA_real_
  weather <- weather[do.call(order, weather[s0]), s2, drop = FALSE]
  row.names(weather) <- NULL
  rm(lastday, append, s0, s1, s2)
} else {
  ## Data frame to fill with country means
  weather <- data.frame(
    country_iso_alpha3 = gl(nlevels(coords$country_iso_alpha3), length(Dates), labels = levels(coords$country_iso_alpha3)),
    Date = rep.int(Dates, N)
  )
  weather[varnames$mikael] <- NA_real_

  ## Grid points common to all ncdf4 objects
  url <- sprintf("%s/t1o5m_mean/global_daily_t1o5m_mean_20200101.nc", parent)
  path <- "tmp.nc"
  download.file(url, path)
  nc <- nc_open(path)
  X <- c(ncvar_get(nc, "longitude"))
  Y <- c(ncvar_get(nc, "latitude"))
  dX <- unique(diff(X))
  dY <- unique(diff(Y))
  X <- seq.int(  0, 360, by = dX)
  Y <- seq.int(-90,  90, by = dY)
  nc_close(nc)
  file.remove(path)
  rm(url, path, nc)
}

## Dates on which country means haven't already been computed
Dates_missing <- lapply(weather[varnames$mikael], function(x) unique(weather$Date[is.na(x)]))

#' Monte Carlo integrator in (longitude, latitude) domain
#' @param n sample size
#' @param x0,x1,y0,y1 limits of box to be sampled uniformly
#' @param z function z(x, y) over domain
#' @param X,Y,Z replacement for `z` in special case: `X`, `Y` are
#'   vectors of grid points in longitude, latitude dimensions;
#'   `Z` is a matrix such that `Z[i, j]` is the value of `z` in box
#'   bounded by `X[i+0:1]`, `Y[i+0:1]`
#' @param f function returning sample statistic
#' @param ... optional arguments to `f`
#' @details
#' longitude is assumed to be `[0, 360)`-valued,
#' latitude is assumed to be `[-90, 90]`-valued
#' @value
#' `f(z(x, y), ...)` where `cbind(x, y)` are the sampled points
mc <- function(n, x0, x1, y0, y1, z = NULL, X, Y, Z, f = mean, ...) {
  x <- runif(n, x0, x1) %% 360
  y <- runif(n, y0, y1)
  abs_y_g90 <- abs(y) > 90
  x[abs_y_g90] <- (x[abs_y_g90] + 180) %% 360
  y[abs_y_g90] <- sign(y[abs_y_g90]) * 180 - y[abs_y_90]
  if (is.null(z)) {
    i <- .bincode(x, breaks = X, right = FALSE, include.lowest = FALSE)
    j <- .bincode(y, breaks = Y, right = FALSE, include.lowest = TRUE)
    f(Z[cbind(i, j)], ...)
  } else {
    f(z(x, y), ...)
  }
}


### Do stuff ###

set.seed(960850250L)
n <- 1000L
scale <- 1.5

for (i in seq_len(nrow(varnames))) {
  v_mikael <- varnames$mikael[i]
  v_metoffice <- varnames$metoffice[i]
  v_ncdf4 <- varnames$ncdf4[i]

  Dm <- Dates_missing[[v_mikael]]
  Dm_Ymd <- format(Dm, "%Y%m%d")
  urls <- sprintf("%s/%s_mean/global_daily_%s_mean_%s.nc", parent, v_metoffice, v_metoffice, Dm_Ymd)
  paths <- sprintf("../weather/%s_mean/%s.nc", v_metoffice, basename(urls))

  ## Subset netCDF files that actually exist at the URL ...
  ## to save time, only test files from the past 7 days
  ## (on April 20, 2021, files up to April 18, 2021 were available)
  e <- rep_len(TRUE, length(urls))
  check <- Dm > yesterday - 7
  e[check] <- vapply(urls[check], function(x) isTRUE(HEAD(x)$status_code == 200L), FALSE)
  if (!any(e)) {
    next
  }
  Dm <- Dm[e]
  urls <- urls[e]
  paths <- paths[e]

  ## Download
  download.file(urls, paths, method = "libcurl")

  for (j in seq_along(paths)) {
    nc <- nc_open(paths[j])
    Z <- ncvar_get(nc, v_ncdf4)
    coords$v <- mcmapply(mc,
      n  = n,
      x0 = coords$longitude - scale * dX,
      x1 = coords$longitude + scale * dX,
      y0 = coords$latitude  - scale * dY,
      y1 = coords$latitude  + scale * dY,
      X  = list(X),
      Y  = list(Y),
      Z  = list(Z),
      mc.preschedule = TRUE,
      mc.cores = 2L
    )
    f <- function(d) weighted.mean(d$v, d$population)
    weather[weather$Date == Dm[j], v_mikael] <-
      c(by(coords, coords$country_iso_alpha3, f))
  }

  ## Clean up
  file.remove(paths)

  ## Save
  save(weather, X, Y, dX, dY, file = path_to_weather)
}
