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

geom_mean <- function(x, na.rm = FALSE, zero.rm = FALSE) {
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero.rm) {
    exp(mean(x[x > 0], na.rm = na.rm))
  } else {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  }
}

stat_mode <- function(x, na.rm = FALSE) {
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

locf <- function(x, x0 = NULL, period = 1L) {
  if (period >= 2L) {
    index <- gl(period, 1L, length(x))
    split(x, index, drop = TRUE) <-
      lapply(split(x, index, drop = TRUE), locf, x0 = x0)
    return(x)
  }
  if (!anyNA(x)) {
    return(x)
  }
  rle_x <- rle_patch(x)
  y <- rle_x$values
  if (is.na(y[1L]) && !is.null(x0)) {
    y[1L] <- x0
  }
  if (anyNA(y[-1L])) {
    argna_y <- which(c(FALSE, is.na(y[-1L])))
    y[argna_y] <- y[argna_y - 1L]
  }
  rle_x$values <- y
  inverse.rle(rle_x)
}

linear <- function(x, x0 = NULL, x1 = NULL, period = 1L, geom = FALSE) {
  if (period >= 2L) {
    index <- gl(period, 1L, length(x))
    split(x, index, drop = TRUE) <-
      lapply(split(x, index, drop = TRUE), linear, x0 = x0, x1 = x1, geom = geom)
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
    rule = 2L
  )
  x[argna] <- if (geom) exp(l$y) else l$y
  x
}

summarize <- function(d, i, func, ...) {
  di <- d[i, , drop = FALSE]
  f <- function(x) sum(!is.na(x))
  c(vapply(di, func, 0, ...), n = vapply(di, f, 0))
}

get_summary_ <- function(d, d_Date, start, func, ...,
                         lag, k, method, x0, x1, period, geom) {
  if (length(start) == 0L || nrow(d) == 0L) {
    X <- matrix(NA_real_,
      nrow = length(start),
      ncol = 2L * length(d),
      dimnames = list(NULL, c(names(d), paste0("n.", names(d))))
    )
    return(X)
  }
  if (k %% 2 == 0) {
    kL <- k / 2 - 1
    kR <- k / 2
  } else {
    kL <- kR <- (k - 1) / 2
  }
  first <- start - lag - kL
  last  <- start - lag + kR
  ab <- seq(min(first, d_Date), max(last, d_Date), by = 1)
  d <- d[match(ab, d_Date), , drop = FALSE]
  if (method != "none") {
    d[] <- switch(method,
      locf = lapply(d, locf, x0 = x0, period = period),
      linear = lapply(d, linear, x0 = x0, x1 = x1, period = period, geom = geom)
    )
  }
  l <- lapply(match(first, ab), seq.int, length.out = k)
  t(vapply(l, summarize, rep_len(0, 2L * length(d)), d = d, func = func, ...))
}

get_summary <- function(d, d_Date, d_index, start, start_index, func, ...,
                        lag = 14L, k = 14L,
                        method = c("none", "locf", "linear"),
                        x0 = NULL, x1 = NULL,
                        period = 1L, geom = FALSE) {
  stopifnot(lag >= 0, lag %% 1 == 0, k > 0, k %% 1 == 0)
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

