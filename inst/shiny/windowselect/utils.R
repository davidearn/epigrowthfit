ymd <- function(x, which = "ymd", drop = TRUE) {
  stopifnot(
    inherits(x, "Date") || is.character(x),
    is.character(which),
    length(which) == 1L,
    !is.na(which),
    is.logical(drop),
    length(drop) == 1L,
    !is.na(drop)
  )
  X <- matrix(NA_integer_, nrow = length(x), ncol = 3L, dimnames = list(names(x), c("y", "m", "d")))
  if (length(x) > 0L) {
    ok <- is.finite(x)
    i <- as.integer(unlist(strsplit(as.character(x[ok]), "-")))
    X[ok, ] <- matrix(i, nrow = sum(ok), ncol = 3L, byrow = TRUE)
  }
  j <- match(strsplit(which, "")[[1L]], colnames(X), 0L)
  X[, j, drop = drop]
}

Dceiling <- function(x, to = c("day", "month", "year")) {
  stopifnot(inherits(x, "Date"))
  to <- match.arg(to)
  x <- .Date(ceiling(unclass(x)))
  if (to == "day" || !any(ok <- is.finite(x))) {
    return(x)
  }
  X <- as.data.frame(ymd(x[ok], drop = FALSE))
  if (to == "month") {
    X$m <- X$m + (X$d > 1L)
    X$y <- X$y + (i <- X$m == 13L)
    X$m[i] <- 1L
    x[ok] <- as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # to == "year"
    X$y <- X$y + (X$m > 1L || X$d > 1L)
    x[ok] <- as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
  x
}

Dfloor <- function(x, to = c("day", "month", "year")) {
  stopifnot(inherits(x, "Date"))
  to <- match.arg(to)
  x <- .Date(floor(unclass(x)))
  if (to == "day" || !any(ok <- is.finite(x))) {
    return(x)
  }
  X <- as.data.frame(ymd(x[ok], drop = FALSE))
  if (to == "month") {
    x[ok] <- as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # to == "year"
    x[ok] <- as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
  x
}

baxis <- function(side, a = NULL, b = NULL, at = NULL, labels = TRUE, ...) {
  stopifnot((side <- as.integer(side)) %in% 1:4)
  dots <- list(...)

  if (is.null(a) || is.null(b)) {
    gp <- par("usr", "xlog", "ylog")
    if (side %% 2L == 1L) {
      usr <- gp$usr[1:2]
      log <- gp$xlog
    } else {
      usr <- gp$usr[3:4]
      log <- gp$ylog
    }
    if (log) {
      usr <- 10^usr
    }
    if (is.null(a)) {
      a <- usr[1L]
    }
    if (is.null(b)) {
      b <- usr[2L]
    }
  }

  args <- list(side = side, at = c(a, b), labels = c("", ""))
  do.call(axis, c(args, replace(dots, "lwd.ticks", list(0))))

  args <- list(side = side, at = at, labels = labels)
  do.call(axis, c(args, replace(dots, c("lwd", "lwd.ticks"), list(0, dots[["lwd.ticks"]]))))
}

Daxis <- function(side, origin = .Date(0), minor = list(), major = list()) {
  stopifnot(
    (side <- as.integer(side)) %in% 1:4,
    inherits(origin, "Date"),
    length(origin) == 1L,
    is.finite(origin),
    is.list(minor) || is.null(minor),
    is.list(major) || is.null(major)
  )
  usr <- par("usr")[if (side %% 2L == 1L) 1:2 else 3:4]
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

    major_at <- double(0L)
    major <- NULL
  }

  ## Minor axis
  if (!is.null(minor)) {
    args <- list(
      side = side,
      at = t0 + minor_at,
      labels = minor_labels
    )
    do.call(baxis, c(args, minor))
  }
  ## Major axis
  if (!is.null(major)) {
    args <- list(
      side = side,
      at = t0 + major_at,
      labels = major_labels
    )
    do.call(baxis, c(args, major))
  }
  list(minor = t0 + minor_at, major = t0 + major_at)
}
