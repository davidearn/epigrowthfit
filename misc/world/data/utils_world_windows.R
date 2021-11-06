## Replacement for 'rle' that treats 'NA' and 'NaN' like other constants
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

## Impute by last observation carried forward
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
    return(mean.default(x, na.rm = na.rm))
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
    return(if (is.double(x)) NaN else x[NA_integer_])
  }
  argna_x <- is.na(x)
  y <- x[!argna_x]
  if (length(y) == 0L) {
    return(if (na.rm && is.double(x)) NaN else x[NA_integer_])
  }
  u <- unique(y)
  f <- tabulate(match(y, u))
  max_f <- max(f)
  argmax_f <- f == max_f
  if (!na.rm && max_f <= max(c(0L, f[!argmax_f])) + sum(argna_x)) {
    return(x[NA_integer_])
  }
  u[argmax_f]
}

## Impute by linear interpolation with a geometric option
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

.get_summary <- function(data, date, start, end, presummarize, summarize) {
  date1 <- seq(min(start, date), max(end, date), 1)
  data1 <- data[match(date1, date, NA_integer_), , drop = FALSE]
  data1[] <- lapply(data1, presummarize)
  subset <- Map(seq.int, match(start, date1), match(end, date1))

  f <- function(i) {
    d <- data1[i, , drop = FALSE]
    list(
      x = data.frame(lapply(d, summarize)),
      n = vapply(d, function(x) sum(!is.na(x)), 0L)
    )
  }

  l <- lapply(subset, f)
  list(
    X = do.call(rbind, lapply(l, `[[`, "x")),
    N = as.data.frame(matrix(unlist(lapply(l, `[[`, "n")), nrow = length(l), ncol = length(data), byrow = TRUE, dimnames = list(NULL, names(data))))
  )
}

#' @param data
#'   A data frame.
#' @param date
#'   A Date vector of length \code{nrow(data)}.
#' @param start,end
#'   Date vectors of equal length specifying summary windows.
#' @param f1
#'   A factor grouping the rows of \code{data} and elements of \code{date}.
#'   Levels not belonging to \code{levels(index2)} are dropped.
#' @param f2
#'   A factor grouping the elements of \code{start} and \code{end}.
#' @param presummarize
#'   A length-preserving function applied groupwise to \code{data[[i]]}
#'   before \code{summarize} is called, enabling, e.g., imputation or
#'   smoothing.
#' @param summarize
#'   A function returning a length-one summary statistic when passed
#'   any subset of \code{data[[i]]}.
get_summary <- function(data, date, start, end, f1, f2, presummarize = identity, summarize = mean) {
  f2 <- as.factor(f2)
  f1 <- factor(f1, levels = levels(f2))
  l <- Map(.get_summary,
    data = split(data, f1),
    date = split(date, f1),
    start = split(start, f2),
    end = split(end, f2),
    presummarize = list(presummarize),
    summarize = list(summarize)
  )
  i <- rep.int(NA_integer_, length(f2))
  res <- list(
    X = l[[1L]][["X"]][i, , drop = FALSE],
    N = l[[1L]][["N"]][i, , drop = FALSE]
  )
  split(res[["X"]], f2) <- lapply(l, `[[`, "X")
  split(res[["N"]], f2) <- lapply(l, `[[`, "N")
  res
}

## Extract index of last observation from an atomic vector
lo <- function(x, value = FALSE) {
  i <- length(x)
  while (i > 0L && is.na(x[i])) {
    i <- i - 1L
  }
  if (i == 0L) {
    i <- NA_integer_
  }
  if (value) x[i] else i
}

## Extract index of first observation from an atomic vector
fo <- function(x, value = FALSE, min = NULL) {
  n <- length(x)
  i <- 1L
  if (is.null(min)) {
    while (i <= n && is.na(x[i])) {
      i <- i + 1L
    }
  } else {
    while (i <= n && (is.na(x[i]) || x[i] < min)) {
      i <- i + 1L
    }
  }
  if (i == n + 1L) {
    i <- NA_integer_
  }
  if (value) x[i] else i
}
