ymd <- function(x, which = "ymd", drop = TRUE) {
  stopifnot(
    is.character(which),
    length(which) == 1L,
    !is.na(which),
    is.logical(drop),
    length(drop) == 1L,
    !is.na(drop)
  )
  X <- matrix(as.integer(unlist(strsplit(as.character(x), "-"), FALSE, FALSE)),
              nrow = length(x),
              ncol = 3L,
              byrow = TRUE,
              dimnames = list(NULL, c("y", "m", "d"))
  )
  j <- unique(match(strsplit(which, "")[[1L]], colnames(X), 0L))
  X[, j, drop = drop]
}

Dceiling <- function(x, to = c("day", "month", "year")) {
  stopifnot(inherits(x, "Date"))
  to <- match.arg(to)
  if (length(x) == 0L) {
    return(.Date(numeric(0L)))
  }
  x <- .Date(ceiling(unclass(x)))
  if (to == "day") {
    return(x)
  }
  X <- as.data.frame(ymd(x, drop = FALSE))
  if (to == "month") {
    X$m <- X$m + (X$d > 1L)
    X$y <- X$y + (i <- X$m == 13L)
    X$m[i] <- 1L
    as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # to == "year"
    X$y <- X$y + (X$m > 1L || X$d > 1L)
    as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
}

Dfloor <- function(x, to = c("day", "month", "year")) {
  stopifnot(inherits(x, "Date"))
  to <- match.arg(to)
  if (length(x) == 0L) {
    return(.Date(numeric(0L)))
  }
  x <- .Date(floor(unclass(x)))
  if (to == "day") {
    return(x)
  }
  X <- as.data.frame(ymd(x, drop = FALSE))
  if (to == "month") {
    as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # to == "year"
    as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
}

Daxis <- function(origin = .Date(0), minor = NULL, major = NULL,
                  show_minor = TRUE, show_major = TRUE) {
  usr <- par("usr")[1:2]
  Dusr <- origin + usr
  D0 <- min(Dceiling(Dusr[1L]), Dfloor(Dusr[2L]))
  D1 <- max(Dceiling(Dusr[1L]), Dfloor(Dusr[2L]), D0 + 1)
  t0 <- julian(D0, origin = origin)
  t1 <- julian(D1, origin = origin)
  w <- t1 - t0

  ## Determine tick coordinates and labels
  if (w <= 210) {
    ## Days
    by <- c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1L]
    minor_at_as_Date <- seq(D0, D1, by = by)
    minor_at <- julian(minor_at_as_Date, origin = D0)
    minor_labels <- ymd(minor_at_as_Date, "d")

    ## Months
    if (ymd(D0, "m") == ymd(D1, "m")) {
      major_at_as_Date <- D0
      major_at <- 0
    } else {
      major_at_as_Date <- seq(Dceiling(D0, "m"), D1, by = "m")
      major_at <- julian(major_at_as_Date, origin = D0)
      if (major_at[1L] > w / 8) {
        major_at_as_Date <- c(D0, major_at_as_Date)
        major_at <- c(0, major_at)
      }
    }
    major_labels <- months(major_at_as_Date, abbreviate = TRUE)
  }
  else if (w <= 3 * 365) {
    ## Months
    by <- c(1L, 2L, 3L)[w <= c(1, 2, 3) * 365][1L]
    minor_at_as_Date <- seq(Dceiling(D0, "m"), D1, by = paste(by, "m"))
    minor_at <- julian(minor_at_as_Date, origin = D0)
    minor_labels <- months(minor_at_as_Date, abbreviate = TRUE)

    ## Years
    if (ymd(D0, "y") == ymd(D1, "y")) {
      major_at_as_Date <- D0
      major_at <- 0
    } else {
      major_at_as_Date <- seq(Dceiling(D0, "y"), D1, by = "y")
      major_at <- julian(major_at_as_Date, origin = D0)
      if (major_at[1L] > w / 8) {
        major_at_as_Date <- c(D0, major_at_as_Date)
        major_at <- c(0, major_at)
      }
    }
    major_labels <- ymd(major_at_as_Date, "y")
  } else {
    ## Years
    by <- ceiling(ceiling(w / 365) / 7)
    minor_at_as_Date <- seq(Dceiling(D0, "year"), D1 + (by + 1) * 365, by = paste(by, "y"))
    minor_at <- julian(minor_at_as_Date, origin = D0)
    minor_labels <- ymd(minor_at_as_Date, "y")
    minor_at <- c(minor_at, (minor_at[-1L] + minor_at[-length(minor_at)]) / 2)
    length(minor_labels) <- length(minor_at)

    major_at <- numeric(0L)
    show_major <- FALSE
  }

  ## Minor axis
  if (show_minor) {
    args <- list(
      side = 1,
      at = t0 + minor_at,
      labels = minor_labels
    )
    do.call(axis, c(args, minor))
  }
  ## Major axis
  if (show_major) {
    args <- list(
      side = 1,
      at = t0 + major_at,
      labels = major_labels
    )
    do.call(axis, c(args, major))
  }
  list(minor = t0 + minor_at, major = t0 + major_at)
}

