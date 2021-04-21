library("parallel")
library("ncdf4")
library("httr")
RNGkind("L'Ecuyer-CMRG")

## Data frame supplying latitude, longitude, population
## of major cities
load("../coords.RData")
coords <- coords[order(coords$country_iso_alpha3), , drop = FALSE] # just in case
cia3 <- levels(coords$country_iso_alpha3)
N <- length(cia3)

## Parent directory for Met Office downloads
parent <- "https://metdatasa.blob.core.windows.net/covid19-response/metoffice_global_daily"

## Dates on which country means are desired
yesterday <- Sys.Date() - 1
Dates <- seq(as.Date("2020-01-01"), yesterday, by = 1)
Dates_Ymd <- format(Dates, "%Y%m%d")

## Variables on which country means are desired
varnames <- c(
  # metoffice = "ncdf4"
  t1o5m     = "air_temperature",
  sh        = "specific_humidity",
  sw        = "m01s01i202",
  precip    = "precipitation_flux",
  windspeed = "wind_speed"
)

## Directory structure
dirnames <- sprintf("../weather/%s_mean", names(varnames))
for (d in dirnames[!dir.exists(dirnames)]) {
  dir.create(d, recursive = TRUE)
}

## Grid points common to all ncdf4 objects
path_to_grid <- "../weather/grid.RData"
if (file.exists(path_to_grid)) {
  load(path_to_grid) # X Y dX dY
} else {
  url <- sprintf("%s/t1o5m_mean/global_daily_t1o5m_mean_20200101.nc", parent)
  path <- "../weather/tmp.nc"
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
  y[abs_y_g90] <- sign(y[abs_y_g90]) * 180 - y[abs_y_g90]
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
threads <- 2L
scale <- 1.5

for (i in seq_along(varnames)) {
  v_metoffice <- names(varnames)[i]
  v_ncdf4 <- varnames[[i]]

  urls <- sprintf("%s/%s_mean/global_daily_%s_mean_%s.nc", parent, v_metoffice, v_metoffice, Dates_Ymd)
  paths_to_nc <- sprintf("../weather/%s_mean/%s", v_metoffice, basename(urls))
  paths_to_RData <- sub("nc$", "RData", paths_to_nc)

  ## Subset netCDF files that haven't already been processed
  l <- !file.exists(paths_to_RData)

  ## Subset netCDF files that actually exist at the URL ...
  ## to save time, only test files from the past 7 days
  ## (on April 20, 2021, files up to April 18, 2021 were available)
  check <- l & Dates > yesterday - 7
  l[check] <- vapply(urls[check], function(x) isTRUE(HEAD(x)$status_code == 200L), FALSE)
  if (!any(l)) {
    next
  }
  urls <- urls[l]
  paths_to_nc <- paths_to_nc[l]
  paths_to_RData <- paths_to_RData[l]

  ## Download
  download.file(urls, paths_to_nc, method = "libcurl")

  n <- length(paths_to_nc)
  for (j in seq_len(n)) {
    cat("Processing", sQuote(v_metoffice), "netCDF file", sprintf("%*d", nchar(n), j), "of", n, "...\n")
    nc <- nc_open(paths_to_nc[j])
    Z <- ncvar_get(nc, v_ncdf4)
    nc_close(nc)
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
      mc.cores = threads
    )
    f <- function(d) weighted.mean(d$v, d$population)
    x <- c(by(coords, coords$country_iso_alpha3, f))
    save(x, file = paths_to_RData[j])
    file.remove(paths_to_nc[j])
  }
}
