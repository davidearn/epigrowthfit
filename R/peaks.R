peaks <- function(x, q = 20, span = NULL, w = 5) {
  lo <- stats::loess(x ~ t,
    data = data.frame(t = seq_along(x), x),
    span = if (is.null(span)) q / length(x) else span
  )
  xbar <- lo$fitted
  w <- floor(w)
  pad <- rep(NA, w)
  xbar_padded <- c(pad, xbar, pad)
  bands <- stats::embed(xbar_padded, 2 * w + 1)

  ## Central element of each band
  bands_mid <- bands[, w+1]
  ## Remaining elements
  bands_sides <- bands[, -(w+1)]
  ## Maximum of remaining elements
  bands_sides_max <- apply(bands_sides, 1, max)

  list(x = x, xbar = xbar, peaks = which(bands_mid > bands_sides_max))
}
