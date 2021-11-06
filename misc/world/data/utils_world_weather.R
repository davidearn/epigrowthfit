#' Test availability of web resources concurrently
#' @description
#' An analogue of \code{file.exists} for URLs.
#' @param url
#'   A character vector listing URLs.
#' @param timeout
#'   A non-negative number giving the time in seconds to wait for results.
#' @param con_total,con_host
#'   Positive integers indicating the maximum number of concurrent
#'   HTTP requests, in total and for any given host.
#' @return
#' A logical vector indexing the elements of \code{url}
#' with response status code 200 (success).
#' @noRd
url.exists <- function(url,
                       timeout = getOption("timeout", 60),
                       con_total = getOption("con_total", 100L),
                       con_host = getOption("con_host", 6L)) {
  stopifnot(requireNamespace("curl"))
  res <- logical(length(url))
  names(res) <- url
  done <- function(req) {
    res[[req[["url"]]]] <<- identical(req[["status_code"]], 200L)
  }
  handle <- function() {
    curl::new_handle(nobody = TRUE)
  }
  pool <- curl::new_pool(total_con = con_total, host_con = con_host)
  for (i in seq_along(url)) {
    curl::curl_fetch_multi(url[i], done = done, handle = handle(), pool = pool)
  }
  curl::multi_run(timeout = timeout, poll = FALSE, pool = pool)
  res
}

#' Download web resources concurrently
#' @description
#' A wrapper function using \code{download.file(method = "libcurl")}
#' to download files in batches.
#' @param url
#'   A character vector listing URLs.
#' @param destfile
#'   A character vector listing paths where downloaded files should be saved.
#'   Tilde expansion is performed.
#' @param batchsize
#'   A positive integer indicating a maximum number of concurrent downloads.
#' @return
#' An integer vector of length \code{ceiling(length(url) / batchsize)}
#' containing the status code returned by each call to \code{download.file}.
batch.download.file <- function(url, destfile, batchsize = 6L) {
  stopifnot(capabilities("libcurl"))
  n <- length(url)
  if (batchsize >= n) {
    return(download.file(url, destfile, method = "libcurl"))
  }
  index <- seq_len(n)
  batch <- as.factor(as.integer(ceiling(index / batchsize)))
  l <- split(index, batch)
  res <- integer(length(l))
  for (i in seq_along(l)) {
    res[i] <- download.file(url[l[[i]]] , destfile[l[[i]]], method = "libcurl")
  }
  res
}

#' Uniform sampling in (longitude, latitude) space
#' @description
#' Samples uniformly from a box in (longitude, latitude) space.
#' Longitude and latitude are assumed to be [0, 360)- and [-90, 90]-valued.
#' @param n Sample size.
#' @param x0,x1,y0,y1 Limits of box.
#' @details
#' \code{x0 < 0}, \code{x1 > 360}, \code{y0 < -90}, and \code{y1 > 90}
#' are all tolerated, as wrapping is implemented in a sensible way.
#' @return
#' \code{cbind(longitude, latitude)}
#' @noRd
runif_ll <- function(n, x0 = 0, x1 = 360, y0 = -90, y1 = 90, ...) {
  x <- runif(n, x0, x1) %% 360
  y <- runif(n, y0, y1)
  k <- abs(y) > 90
  x[k] <- (x[k] + 180) %% 360
  y[k] <- sign(y[k]) * 180 - y[k]
  cbind(x, y, deparse.level = 0L)
}

#' Process Met Office weather data file
#' @description
#' For each of \code{N} cities, sample \code{n} points from a box in
#' (longitude, latitude) space and retrieve the value at each point
#' of a weather variable from a Met Office NetCDF file. Compute from
#' the resulting collection of values the city population-weighted
#' mean of the variable within countries.
#' @param path
#'   A character string indicating a path to a NetCDF (usually \code{.nc})
#'   file.
#' @param url
#'   A character string indicating a URL from which the NetCDF file can be
#'   downloaded in the event that a file does not exist at \code{path} or
#'   the existing file is corrupt.
#' @param varid
#'   A character string indicating a variable to read from the NetCDF file.
#' @param data
#'   A data frame with variables \code{country},
#'   \code{longitude}, \code{latitude}, and \code{population},
#'   supplying those details about \code{N = nrow(data)} cities.
#' @param n
#'   A positive integer indicating a sample size.
#' @param scale
#'   A positive number indicating a ratio of box side length to NetCDF grid
#'   spacing.
#' @return
#' A double vector corresponding elementwise to
#' \code{levels(as.factor(data[["country"]]))},
#' giving the city population-weighted mean of \code{varid} in each country.
#' @noRd
nc_process <- function(path, url = NULL, varid, data, n = 1000L, scale = 1) {
  stopifnot(requireNamespace("ncdf4"))
  cat("Processing NetCDF file:", path, "...\n")
  nc <- tryCatch(
    expr = {
      ncdf4::nc_open(path)
    },
    error = function(e) {
      if (is.null(url)) {
        stop(e)
      }
      download.file(url, path)
      eval(conditionCall(e))
    }
  )
  on.exit(ncdf4::nc_close(nc))
  z <- ncdf4::ncvar_get(nc, varid)
  x <- c(ncdf4::ncvar_get(nc, "longitude"))
  y <- c(ncdf4::ncvar_get(nc, "latitude"))
  dx <- 0.5 * scale * (x[2L] - x[1L])
  dy <- 0.5 * scale * (y[2L] - y[1L])

  ## FIXME: memory limits could be an issue when 'n * N' is , see 'help("Memory-limits")'

  R <- runif_ll(
    n = n * nrow(data),
    x0 = data[["longitude"]] - dx,
    x1 = data[["longitude"]] + dx,
    y0 = data[["latitude"]]  - dy,
    y1 = data[["latitude"]]  + dy
  )
  x[c(1L, length(x))] <- c(0, 360)
  y[c(1L, length(y))] <- c(-90, 90)
  i <- .bincode(R[, 1L], breaks = x, right = FALSE, include.lowest = FALSE)
  j <- .bincode(R[, 2L], breaks = y, right = FALSE, include.lowest = TRUE)

  fxw <- data.frame(
    f = data[["country"]],
    x = z[cbind(i, j)],
    w = data[["population"]]
  )
  c(by(fxw[c("x", "w")], fxw["f"], function(d) do.call(weighted.mean, d)))
}

#' Process Met Office weather data
#' @description
#' Downloads and processes NetCDF files from a Met Office server
#' and caches output for later reuse.
#' NetCDF files are deleted once output is generated.
#' @param path
#'   A character string containing the path to a root directory
#'   where downloads are (temporarily) saved and output is cached.
#' @param date
#'   A Date vector listing dates for which to download weather data.
#'   The Met Office provides data from January 1, 2020.
#' @param varid
#'   A character vector listing names of weather variables.
#' @param data,n,scale
#'   Arguments to \code{process_nc}.
#' @return
#' \code{NULL}, invisibly.
update_weather <- function(path,
                           date = seq(as.Date("2020-01-01"), Sys.Date() - 1, 1),
                           varid = c("temperature", "specific_humidity", "shortwave_radiation", "precipitation", "wind_speed"),
                           data,
                           n = 1000L,
                           scale = 1) {
  stopifnot(
    is.character(path),
    length(path) == 1L,
    dir.exists(path),
    inherits(date, "Date"),
    length(date) > 0L,
    !anyNA(date),
    is.data.frame(data),
    c("country", "longitude", "latitude", "population") %in% names(data),
    is.numeric(n),
    length(n) == 1L,
    n >= 1,
    is.numeric(scale),
    length(scale) == 1L,
    scale > 0
  )

  ymd <- format(date, "%Y%m%d")
  map <- cbind(
    metoffice = c("t1o5m", "sh", "sw", "precip", "windspeed"),
    netcdf = c("air_temperature", "specific_humidity", "m01s01i202", "precipitation_flux", "wind_speed")
  )
  rownames(map) <- eval(formals(sys.function())[["varid"]])
  varid <- unique(match.arg(varid, rownames(map), several.ok = TRUE))
  url <- "https://metdatasa.blob.core.windows.net/covid19-response/metoffice_global_daily"

  for (i in seq_along(varid)) {
    subdir <- file.path(path, varid[i])
    if (!dir.exists(subdir)) {
      dir.create(subdir)
    }

    file <- data.frame(
      ## Addresses of source files
      url = file.path(url, paste0(map[varid[i], "metoffice"], "_mean"), paste0("global_daily_", map[varid[i], "metoffice"], "_mean_", ymd, ".nc"), fsep = "/"),
      ## Paths to source files
      nc  = file.path(subdir, paste0(varid[i], "_", ymd, ".nc")),
      ## Paths to output files
      rds = file.path(subdir, paste0(varid[i], "_", ymd, ".rds")),
      stringsAsFactors = FALSE
    )

    ## Subset source files that have not already been processed
    e <- file.exists(file[["rds"]])
    if (all(e)) {
      next
    }
    file <- file[!e, , drop = FALSE]

    ## Subset source files that are actually available for download
    e <- url.exists(file[["url"]])
    if (!any(e)) {
      next
    }
    file <- file[e, , drop = FALSE]

    ## Identify source files already in file system
    e <- file.exists(file[["nc"]])

    ## Download missing source files in batches
    if (!all(e)) {
      batch.download.file(file[!e, "url"], file[!e, "nc"], batchsize = getOption("con_host", 6L))
    }

    ## FIXME:
    ## Possible disk space issue since all NetCDF files are temporarily
    ## but simultaneously stored on disk, at least while there is no cache.
    ## Fix would involve downloading and processing _in the same loop_.
    f <- function(k) {
      res <- nc_process(
        path = file[k, "nc"],
        url = file[k, "url"],
        varid = map[varid[i], "netcdf"],
        data = data,
        n = n,
        scale = scale
      )
      saveRDS(res, file = file[k, "rds"])
      file.remove(file[k, "nc"])
    }
    parallel::mclapply(seq_len(nrow(file)), f)
  }
  invisible(NULL)
}

.shape_weather <- function(path, date1 = seq(as.Date("2020-01-01"), Sys.Date() - 1, 1)) {
  lf <- list.files(path, pattern = "_\\d{8}\\.rds$", full.names = TRUE, recursive = FALSE)
  if (length(lf) == 0L) {
    return(NULL)
  }
  l <- lapply(lf, readRDS) # list of named double vectors
  date <- as.Date(sub("^.*_(\\d{8})\\.rds$", "\\1", lf), format = "%Y%m%d")
  value <- unlist(l, recursive = FALSE, use.names = TRUE)
  res <- data.frame(
    country = factor(names(value)),
    date = rep.int(date, lengths(l, use.names = FALSE)),
    value = value,
    row.names = NULL
  )
  m <- match(date1, date, 0L)
  i <- m == 0L
  if (!any(i)) {
    return(res)
  }
  s <- levels(res[["country"]])
  res_fill <- data.frame(
    country = gl(length(s), sum(i), labels = s),
    date = date1[i],
    value = NA_real_,
    row.names = NULL
  )
  rbind(res, res_fill)
}

shape_weather <- function(path, date1 = seq(as.Date("2020-01-01"), Sys.Date() - 1, 1)) {
  ld <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  l <- lapply(ld, .shape_weather, date1 = date1)
  res <- do.call(rbind, l)
  res[["varid"]] <- rep.int(gl(length(ld), 1L, labels = basename(ld)), vapply(l, NROW, 0L))
  s <- c("country", "varid", "date")
  o <- do.call(order, unname(res[s]))
  res <- res[o, c(s, "value"), drop = FALSE]
  row.names(res) <- NULL
  res
}
