## Replacement for 'rle' that treats missing values like other constants
rle_patch <- function(x) {
  n <- length(x)
  if (n == 0L) {
    return(list(lengths = integer(0L), values = x))
  }
  l <- x[-n] != x[-1L]
  if (any(argna <- is.na(x))) {
    l[is.na(l)] <- FALSE
    l <- l | (argna[-n] & !argna[-1L]) | (!argna[-n] & argna[-1L])
  }
  i <- c(which(l), n)
  list(lengths = diff(c(0L, i)), values = x[i])
}

## Imputation by last observation carried forward
locf <- function(x, x0 = NULL, period = 1L) {
  if (period >= 2L) {
    index <- gl(period, 1L, length(x))
    split(x, index) <- lapply(split(x, index), locf, x0 = x0)
    return(x)
  }
  if (!anyNA(x)) {
    return(x)
  }
  rle_x <- rle_patch(x)
  y <- rle_x[["values"]]
  if (is.na(y[1L]) && !is.null(x0)) {
    y[1L] <- x0
  }
  if (anyNA(y[-1L])) {
    argna_y <- which(c(FALSE, is.na(y[-1L])))
    y[argna_y] <- y[argna_y - 1L]
  }
  rle_x[["values"]] <- y
  inverse.rle(rle_x)
}

## Replacement for 'mean.default' with a geometric option (but no 'trim')
Mean <- function(x, na.rm = FALSE, geom = FALSE, zero.rm = FALSE) {
  if (!geom) {
    return(mean(x, na.rm = na.rm))
  }
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero.rm) {
    return(exp(mean(x[x > 0], na.rm = na.rm)))
  }
  if (any(x == 0, na.rm = TRUE)) {
    return(0)
  }
  exp(mean(log(x), na.rm = na.rm))
}

## Mode
Mode <- function(x, na.rm = FALSE) {
  if (length(x) == 0L) {
    return(NaN)
  }
  argna_x <- is.na(x)
  y <- x[!argna_x]
  if (length(y) == 0L) {
    return(if (na.rm) NaN else NA)
  }
  u <- unique(y)
  f <- tabulate(match(y, u))
  max_f <- max(f)
  argmax_f <- f == max_f
  if (!na.rm && max_f <= max(c(0L, f[!argmax_f])) + sum(argna_x)) {
    return(NA)
  }
  u[argmax_f]
}

## Imputation by linear interpolation with a geometric option
linear <- function(x, x0 = NULL, x1 = NULL, period = 1L, geom = FALSE) {
  if (period >= 2L) {
    index <- gl(period, 1L, length(x))
    split(x, index) <- lapply(split(x, index), linear, x0 = x0, x1 = x1, geom = geom)
    return(x)
  }
  if (all(argna <- is.na(x))) {
    return(x)
  }
  if (argna[1L] && !is.null(x0)) {
    i <- seq_len(which.min(argna) - 1L)
    x[i] <- x0
    argna[i] <- FALSE
  }
  if (argna[n <- length(x)] && !is.null(x1)) {
    j <- n + 1L - seq_len(which.min(argna[n:1L]) - 1L)
    x[j] <- x1
    argna[j] <- FALSE
  }
  if (!any(argna)) {
    return(x)
  }
  l <- approx(
    x = seq_len(n),
    y = if (geom) log(x) else x,
    xout = which(argna),
    method = "linear",
    rule = 2L
  )
  x[argna] <- if (geom) exp(l[["y"]]) else l[["y"]]
  x
}

.get_summary <- function(data, date, start, end, func, ..., method, x0, x1, period, geom) {
  date1 <- seq(min(start, date), max(end, date), 1)
  data1 <- data[match(date1, date, NA_integer_), , drop = FALSE]
  if (method != "none") {
    data1[] <- switch(method,
      locf = lapply(data1, locf, x0 = x0, period = period),
      linear = lapply(data1, linear, x0 = x0, x1 = x1, period = period, geom = geom)
    )
  }
  subset <- Map(seq.int, match(start, date1), match(end, date1))

  f <- function(data, subset, func, ...) {
    data <- data[subset, , drop = FALSE]
    list(
      x = vapply(data, func, 0, ...),
      n = vapply(data, function(x) sum(!is.na(x)), 0L)
    )
  }

  l <- lapply(subset, f, data = data, func = func, ...)
  res <- sapply(c("x", "n"), function(el) matrix(unlist(lapply(l, `[[`, el)), nrow = length(l), ncol = length(data), byrow = TRUE, dimnames = list(NULL, names(data))), simplify = FALSE)
  names(res) <- toupper(names(res))
  res
}

#' @param data a data frame listing only numeric variables
#' @param date a Date vector of length \code{nrow(data)}
#' @param start,end Date vectors of equal length indicating left and right endpoints of fitting windows
#' @param index1 a factor grouping the rows of \code{data} and elements of \code{date}
#' @param index2 a factor grouping the elements of \code{start} and \code{end}
#' @param func a function computing a summary statistic from subsets of \code{data[[i]]}
#' @param ... optional arguments to \code{func}
#' @param method a character string indicating a method of imputation
#' @param x0,x1,period,geom optional arguments to imputation function
get_summary <- function(data, date, start, end, index1, index2, func, ...,
                        method = c("none", "locf", "linear"),
                        x0 = NULL, x1 = NULL, period = 1L, geom = FALSE) {
  index1 <- as.factor(index1)
  index2 <- as.factor(index2)
  stopifnot(
    is.data.frame(data),
    vapply(data, mode, "") == "numeric",
    inherits(date, "Date"),
    length(date) == nrow(data),
    is.finite(date),
    length(index1) == length(index)
    inherits(start, "Date"),
    is.finite(start),
    inherits(end, "Date"),
    length(end) == length(start),
    end > start,

    length(date) == nrow(data),
    is.finite(date),

  )
  method <- match.arg(method)
  start_index <- droplevels(start_index)
  d_index <- factor(d_index, levels = levels(start_index))
  Xl <- Map(get_summary_,
    d = split(d, d_index),
    d_Date = split(d_Date, d_index),
    start = split(start, start_index),
    func = list(func),
    ...,
    lag = lag,
    k = k,
    method = method,
    x0 = list(x0),
    x1 = list(x1),
    period = period,
    geom = geom
  )
  do.call(rbind, Xl)
}

