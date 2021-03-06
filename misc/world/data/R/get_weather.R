library("ncdf4")
library("httr")
library("parallel")
RNGkind("L'Ecuyer-CMRG")

## Data frame supplying latitude, longitude, population
## of major cities
load("coords.RData")
coords <- coords[order(coords$country_iso_alpha3), , drop = FALSE] # just in case

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
dirnames <- sprintf("weather/%s_mean", names(varnames))
for (d in dirnames[!dir.exists(dirnames)]) {
  dir.create(d, recursive = TRUE)
}

## Grid points common to all ncdf4 objects
path_to_grid <- "weather/grid.RData"
if (file.exists(path_to_grid)) {
  load(path_to_grid) # X Y dX dY
} else {
  url <- sprintf("%s/t1o5m_mean/global_daily_t1o5m_mean_20200101.nc", parent)
  path <- "weather/tmp.nc"
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
  save(X, Y, dX, dY, file = path_to_grid)
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

scale <- 1.5
n <- 1000L
num_simul_dl <- 7L
num_threads <- 2L

set.seed(960850250L)
options(timeout = 60 * num_simul_dl, warn = 1L)

for (i in seq_along(varnames)) {
  v_metoffice <- names(varnames)[i]
  v_ncdf4 <- varnames[[i]]

  ## Paths to output files
  paths_to_RData <- sprintf("weather/%s_mean/global_daily_%s_mean_%s.RData", v_metoffice, v_metoffice, Dates_Ymd)

  ## Subset output files that have not already been generated
  e <- file.exists(paths_to_RData)
  if (all(e)) {
    next
  }
  paths_to_RData <- paths_to_RData[!e]

  ## Web addresses for netCDF input files
  urls <- sprintf("%s/%s_mean/%s", parent, v_metoffice, sub("RData$", "nc", basename(paths_to_RData)))

  ## Test for availability of resource
  e <- vapply(urls, function(x) isTRUE(HEAD(x)$status_code == 200L), FALSE, USE.NAMES = FALSE)
  urls <- urls[e]
  paths_to_RData <- paths_to_RData[e]

  ## Paths to netCDF input files
  paths_to_nc <- sub("RData$", "nc", paths_to_RData)

  ## Subset netCDF input files that are not already in file system
  e <- file.exists(paths_to_nc)

  ## Download in segments
  d <- sum(!e)
  if (d > 0L) {
    urls_for_dl <- urls[!e]
    paths_to_nc_for_dl <- paths_to_nc[!e]
    for (j in split(seq_len(d), ceiling(seq_len(d) / num_simul_dl))) {
      download.file(urls_for_dl[j], paths_to_nc_for_dl[j], method = "libcurl")
    }
  }

  ## Process sequentially, deleting input once output is generated
  p <- length(paths_to_nc)
  f <- function(d) weighted.mean(d$v, d$population)
  for (j in seq_len(p)) {
    cat("Processing", sQuote(v_metoffice), "netCDF file",
        sprintf("%*d", nchar(p), j), "of", p, "...\n")
    nc <- try(nc_open(paths_to_nc[j]))
    if (inherits(nc, "try-error")) { # in case of corrupt initial download
      download.file(urls[j], paths_to_nc[j])
      nc <- nc_open(paths_to_nc[j])
    }
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
      mc.cores = num_threads
    )
    x <- c(by(coords, coords$country_iso_alpha3, f))
    save(x, file = paths_to_RData[j])
    file.remove(paths_to_nc[j])
  }
}
