#' \loadmathjax
#' Generation interval distribution
#'
#' @description
#' Generation interval density, distribution, and quantile functions.
#' Assumes:
#' * that the latent and infectious periods are integer-valued,
#' * that the latent period and infectious waiting time are independent, and
#' * that infectiousness is constant over the infectious period.
#'
#' @param x,q
#'   A numeric vector listing generation intervals.
#' @param p
#'   A numeric vector listing probabilities.
#' @param n
#'   A non-negative integer. The number of samples generated.
#' @param lat
#'   A numeric vector listing the probability \mjseqn{p_i} that the
#'   latent period is \mjseqn{i} (or the corresponding probability
#'   weight) for all \mjseqn{i \in \lbrace 1,\ldots,m \rbrace}.
#' @param inf
#'   A numeric vector listing the probability \mjseqn{q_i} that the
#'   infectious period is \mjseqn{i} (or the corresponding probability
#'   weight) for all \mjseqn{i \in \lbrace 1,\ldots,n \rbrace}.
#'
#' @return
#' A numeric vector with length equal to the length of the
#' first argument, or with length `n` in the case of `rgi()`.
#'
#' `dgi()` evaluates the density function \mjseqn{f_\text{gen}}.
#' `pgi()` evaluates the distribution function \mjseqn{F_\text{gen}}.
#' `qgi()` evaluates the quantile function, which is defined as the
#' left-continuous generalized inverse of `pgi()`.
#' `rgi()` samples from the distribution.
#'
#' @details
#' Let \mjseqn{\tau_\text{lat}} and \mjseqn{\tau_\text{inf}} be the
#' latent and infectious periods, \mjseqn{f_\text{lat}} the latent
#' period density function, and \mjseqn{F_\text{inf}} the infectious
#' period distribution function. Suppose the distributions
#' of \mjseqn{\tau_\text{lat}} and \mjseqn{\tau_\text{inf}}
#' are supported on \mjseqn{\lbrace 1,\ldots,m \rbrace} and
#' \mjseqn{\lbrace 1,\ldots,n \rbrace}, respectively, so that
#'
#' \mjtdeqn{\begin{array}{r@{}ll} f_\text{lat}(t) &= \sum_{i=1}^{m} p_i \delta(t - i)\,, \cr F_\text{inf}(t) & {}= \mathrm{1}USCORE{\lbrack 1,\infty)} \sum_{i=1}^{\min\lbrace n,\lfloor t \rfloor\rbrace} q_i\,, \end{array}}{\begin{align} f_\text{lat}(t) & {}= \sum_{i=1}^{m} p_i \delta(t - i)\,, \cr F_\text{inf}(t) &= \unicode{x1D7D9}USCORE{\lbrack 1,\infty)} \sum_{i=1}^{\min\lbrace n,\lfloor t \rfloor\rbrace} q_i\,, \end{align}}{LaTeX}
#'
#' where \mjseqn{p_i,q_i \in \lbrack 0,1 \rbrack},
#' \mjseqn{\sum_{i=1}^{m} p_i = \sum_{i=1}^{n} q_i = 1}
#' and \mjseqn{\delta} is the Dirac delta.
#'
#' Let \mjseqn{\tau_\text{wait}} be the infectious waiting time and
#' \mjseqn{f_\text{wait}} its density. Then the generation interval is
#' \mjseqn{\tau_\text{gen} = \tau_\text{lat} + \tau_\text{wait}}
#' and has density
#'
#' \mjsdeqn{f_\text{gen}(t) = (f_\text{lat} * f_\text{wait})(t) = \sum_{i=1}^{m} p_i f_\text{wait}(t - i)\,,}
#'
#' assuming independence of \mjseqn{\tau_\text{lat}} and
#' \mjseqn{\tau_\text{wait}}.
#'
#' Equation 5.7 of \insertCite{Sven07;textual}{epigrowthfit}
#' gives an expression for \mjseqn{f_\text{wait}} for the case
#' in which infectiousness is constant over the infectious period:
#'
#' \mjtdeqn{f_\text{wait} = \mathrm{1}USCORE{(0,\infty)} \frac{1 - F_\text{inf}}{\mathrm{E}\lbrack \tau_\text{inf} \rbrack}\,.}{f_\text{wait} = \unicode{x1D7D9}USCORE{(0,\infty)} \frac{1 - F_\text{inf}}{\unicode{x1D53C}\lbrack \tau_\text{inf} \rbrack}\,.}{LaTeX}
#'
#' In this case, it follows that \mjseqn{f_\text{gen}} is supported
#' on the interval \mjseqn{\lbrack 1,m+n)} and constant on the interval
#' \mjseqn{\lbrack i,i+1)} for all integers \mjseqn{i}. This allows
#' the distribution and quantile functions of \mjseqn{\tau_\text{gen}}
#' to be evaluated via linear interpolation.
#'
#' @references
#' \insertRef{Sven07}{epigrowthfit}
#'
#' @examples
#' lat <- inf <- rep(1, 10) # probability weights
#' m <- length(lat)
#' n <- length(inf)
#'
#' x <- seq(0, m + n + 1, by = 0.02)
#' fx <- dgi(x, lat, inf)
#' Fx <- pgi(x, lat, inf)
#'
#' plot(x, fx, type = "l", las = 1,
#'   xlab = "generation interval",
#'   ylab = "density function"
#' )
#' plot(x, Fx, type = "l", las = 1,
#'   xlab = "generation interval",
#'   ylab = "distribution function"
#' )
#'
#' @name gi
NULL

#' @rdname gi
#' @export
dgi <- function(x, lat, inf) {
  stop_if_not(
    is.numeric(x),
    m = "`x` must be numeric."
  )
  if (length(x) == 0L) {
    return(numeric(0L))
  }
  check_gi(lat, inf)
  lat <- lat / sum(lat)
  inf <- inf / sum(inf)
  m <- length(lat)
  n <- length(inf)

  d <- rep.int(NA_real_, length(x))
  d[x < 1 | x >= m + n] <- 0

  is_x_in_1_mpn <- (!is.na(x) & x >= 1 & x < m + n)
  if (any(is_x_in_1_mpn)) {
    ## Take advantage of the fact that the generation interval density
    ## is constant on intervals [i,i+1)
    x_floor <- floor(x[is_x_in_1_mpn])
    x_floor_unique <- unique(x_floor)
    k <- match(x_floor, x_floor_unique)

    ij_dim <- c(m, length(x_floor_unique))
    i <- .row(ij_dim)
    j <- .col(ij_dim)

    d[is_x_in_1_mpn] <- colSums(lat * dwait(x_floor_unique[j] - i, inf = inf))[k]
  }
  dim(d) <- dim(x)
  d
}

#' @rdname gi
#' @export
#' @importFrom stats approx
pgi <- function(q, lat, inf) {
  stop_if_not(
    is.numeric(q),
    m = "`q` must be numeric."
  )
  if (length(q) == 0L) {
    return(numeric(0L))
  }
  m <- length(lat)
  n <- length(inf)

  p <- rep.int(NA_real_, length(q))
  p[q <= 1] <- 0
  p[q >= m + n] <- 1
  is_q_in_1_mpn <- (!is.na(q) & q > 1 & q < m + n)
  if (any(is_q_in_1_mpn)) {
    p[is_q_in_1_mpn] <- approx(
      x = seq_len(m + n),
      y = c(0, cumsum(dgi(seq_len(m + n - 1L), lat = lat, inf = inf))),
      xout = q[is_q_in_1_mpn]
    )$y
  }
  dim(p) <- dim(q)
  p
}

#' @rdname gi
#' @export
#' @importFrom stats approx
qgi <- function(p, lat, inf) {
  stop_if_not(
    is.numeric(p),
    m = "`p` must be numeric."
  )
  if (length(p) == 0L) {
    return(numeric(0L))
  }
  m <- length(lat)
  n <- length(inf)

  q <- rep.int(NA_real_, length(p))
  q[p < 0 | p > 1] <- NaN
  q[p == 0] <- -Inf
  q[p == 1] <- m + n
  is_p_in_0_1 <- (!is.na(p) & p > 0 & p < 1)
  if (any(is_p_in_0_1)) {
    q[is_p_in_0_1] <- approx(
      x = c(0, cumsum(dgi(seq_len(m + n - 1L), lat = lat, inf = inf))),
      y = seq_len(m + n),
      xout = p[is_p_in_0_1]
    )$y
  }
  dim(q) <- dim(p)
  q
}

#' @rdname gi
#' @export
#' @importFrom stats runif
rgi <- function(n, lat, inf) {
  stop_if_not(
    is.numeric(n),
    length(n) == 1L,
    is.finite(n),
    n >= 0,
    m = "`n` must be a non-negative number."
  )
  if (n == 0L) {
    return(numeric(0L))
  }
  check_gi(lat, inf)
  lat <- lat / sum(lat)
  inf <- inf / sum(inf)

  ## Latent period is integer-valued,
  ## equal to  `i` with probability `lat[i]`
  rlat <- sample(seq_along(lat), size = n, replace = TRUE, prob = lat)

  ## Infectious waiting time is real-valued, with
  ## distribution supported on [0,length(inf)) and
  ## with constant density on subintervals [i,i+1)
  rwait_floor <- sample(seq_along(inf) - 1L, size = n, replace = TRUE,
                        prob = dwait(seq_along(inf) - 1L, inf = inf))
  rwait <- runif(n, min = rwait_floor, max = rwait_floor + 1L)

  ## Sum of latent period and infectious waiting time
  ## yields generation interval
  rlat + rwait
}
