#' \loadmathjax
#' The generation interval distribution
#'
#' @description
#' Evaluate the generation interval density and distribution functions
#' implied by discrete distributions of the latent and infectious periods.
#' These functions are valid under the assumptions (i) that the latent
#' period and infectious waiting time are independent and (ii) that
#' infectiousness is constant over the infectious period.
#'
#' @param x A numeric vector listing values for the generation interval
#'   as a number of days.
#' @param lat A numeric vector listing the probability
#'   \mjseqn{p_i} that the latent period is \mjseqn{i} days,
#'   for all \mjseqn{i \in \lbrace 1,\ldots,m \rbrace}.
#'   Replaced with `lat / sum(lat)` in the event that
#'   `sum(lat) != 1`.
#' @param inf A numeric vector listing the probability
#'   \mjseqn{q_i} that the infectious period is \mjseqn{i} days,
#'   for all \mjseqn{i \in \lbrace 1,\ldots,n \rbrace}.
#'   Replaced with `inf / sum(inf)` in the event
#'   that `sum(inf) != 1`.
#'
#' @return
#' A numeric vector with length `length(x)`.
#' For `dgi()`, the `i`th element is the generation interval
#' density function (see \mjseqn{f_\text{gen}} in Details)
#' evaluated at `x[i]`.
#' For `pgi()`, the `i`th element is the generation interval
#' distribution function (see \mjseqn{F_\text{gen}} in Details)
#' evaluated at `x[i]`.
#'
#' @details
#' Let \mjseqn{\tau_\text{lat}} and \mjseqn{\tau_\text{inf}} be
#' the latent and infectious periods, let \mjseqn{f_\text{lat}}
#' and \mjseqn{f_\text{inf}} be their density functions, and let
#' \mjseqn{F_\text{lat}} and \mjseqn{F_\text{inf}} be their
#' distribution functions. Suppose \mjseqn{f_\text{lat}} and
#' \mjseqn{f_\text{inf}} are supported on
#' \mjseqn{\lbrace 1,\ldots,m \rbrace} days
#' and \mjseqn{\lbrace 1,\ldots,n \rbrace} days, respectively,
#' so that
#'
#' \mjtdeqn{\begin{array}{r@{}ll} f_\text{lat}(t) &= \sum_{i=1}^{m} p_i \delta(t - i)\,, \cr f_\text{inf}(t) & {}= \sum_{i=1}^{n} q_i \delta(t - i)\,, \end{array}}{\begin{align} f_\text{lat}(t) & {}= \sum_{i=1}^{m} p_i \delta(t - i)\,, \cr f_\text{inf}(t) &= \sum_{i=1}^{n} q_i \delta(t - i)\,, \end{align}}{LaTeX}
#'
#' where \mjseqn{p_i,q_i \in \lbrack 0,1 \rbrack},
#' \mjseqn{\sum_{i=1}^{m} p_i = \sum_{i=1}^{n} q_i = 1}
#' and \mjseqn{\delta} is the Dirac delta.
#'
#' Let \mjseqn{\tau_\text{wait}} be the infectious waiting time and
#' \mjseqn{f_\text{wait}} its density. Then the generation interval
#' is \mjseqn{\tau_\text{gen} = \tau_\text{lat} + \tau_\text{wait}}
#' and has density function
#' \mjseqn{f_\text{gen} = f_\text{lat} * f_\text{wait}},
#' assuming independence of \mjseqn{\tau_\text{lat}} and
#' \mjseqn{\tau_\text{wait}}.
#'
#' Equation 5.7 of \insertCite{Sven07;textual}{epigrowthfit}
#' gives an expression for \mjseqn{f_\text{wait}}, assuming
#' that infectiousness is constant over the infectious period:
#'
#' \mjtdeqn{f_\text{wait} = \mathrm{1}USCORE{(0,\infty)} \frac{1 - F_\text{inf}}{\mathrm{E}\lbrack \tau_\text{inf} \rbrack}\,.}{f_\text{wait} = \unicode{x1D7D9}USCORE{(0,\infty)} \frac{1 - F_\text{inf}}{\unicode{x1D53C}\lbrack \tau_\text{inf} \rbrack}\,.}{LaTeX}
#'
#' In this case, one can show that
#'
#' \mjtdeqn{f_\text{gen}(t) = \left\lbrace \begin{array}{l@{\qquad}ll} 0\,, & t \in (-\infty,1)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} p_1\,, & t \in \lbrack 1,2)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \big(p_{\lfloor t \rfloor} + \sum_{i=1}^{\lfloor t \rfloor-1} p_i \big(1 - \sum_{j=1}^{\min\lbrace n,\lfloor t \rfloor-i \rbrace} q_j\big)\big)\,, & t \in \lbrack 2,m+1)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \sum_{i=1}^{m} p_i \big(1 - \sum_{j=1}^{\min\lbrace n,\lfloor t \rfloor-i \rbrace} q_j\big)\,, & t \in \lbrack m+1,m+n)\,, \cr 0\,, & t \in \lbrack m+n,\infty)\,. \end{array} \right.}{f_\text{gen}(t) = \begin{cases} 0\,, & t \in (-\infty,1)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} p_1\,, & t \in \lbrack 1,2)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \big(p_{\lfloor t \rfloor} + \sum_{i=1}^{\lfloor t \rfloor-1} p_i \big(1 - \sum_{j=1}^{\min\lbrace n,\lfloor t \rfloor-i \rbrace} q_j\big)\big)\,, & t \in \lbrack 2,m+1)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \sum_{i=1}^{m} p_i \big(1 - \sum_{j=1}^{\min\lbrace n,\lfloor t \rfloor-i \rbrace} q_j\big)\,, & t \in \lbrack m+1,m+n)\,, \cr 0\,, & t \in \lbrack m+n,\infty)\,. \end{cases}}{LaTeX}
#'
#' Note that \mjseqn{f_\text{gen}} is supported on the interval
#' \mjseqn{(1,m+n)} and constant on the interval \mjseqn{(i,i+1)}
#' for all integers \mjseqn{i}. Hence the probability that
#' \mjseqn{\tau_\text{gen} \in (i,i+1\rbrack} is simply
#'
#' \mjsdeqn{g_i = \int_{i}^{i+1} f_\text{gen}(s)\,\text{d}s = \int_{i}^{i+1} f_\text{gen}(i+1/2)\,\text{d}s = f_\text{gen}(i+1/2)\,.}
#'
#' Then the distribution function of \mjseqn{\tau_\text{gen}} is
#'
#' \mjtdeqn{F_\text{gen}(t) = \int_{-\infty}^{t} f_\text{gen}(s)\,\text{d}s = \left\lbrace \begin{array}{l@{\qquad}ll} 0\,, & t \in (-\infty,1\rbrack\,, \cr (t - 1) g_1\,, & t \in (1,2\rbrack\,, \cr (t - i) g_i + \sum_{j=1}^{i-1} g_j\,, & t \in (i,i+1\rbrack\,,\quad i \in \lbrace 2,\ldots,m+n-1 \rbrace\,, \cr 1\,, & t \in (m+n,\infty)\,. \end{array} \right.}{\begin{align} F_\text{gen}(t) &= \int_{-\infty}^{t} f_\text{gen}(s)\,\text{d}s \cr &= \begin{cases} 0\,, & t \in (-\infty,1\rbrack\,, \cr (t - 1) g_1\,, & t \in (1,2\rbrack\,, \cr (t - i) g_i + \sum_{j=1}^{i-1} g_j\,, & t \in (i,i+1\rbrack\,,\quad i \in \lbrace 2,\ldots,m+n-1 \rbrace\,, \cr 1\,, & t \in (m+n,\infty)\,. \end{cases} \end{align}}{LaTeX}
#'
#' @references
#' \insertRef{Sven07}{epigrowthfit}
#'
#' @seealso [compute_R0()]
#'   for computing the basic reproduction number implied by
#'   an initial growth rate and a generation interval distribution
#'
#' @examples
#' data(plague_latent_period)
#' lat <- plague_latent_period$relfreq
#' m <- length(lat)
#'
#' data(plague_infectious_period)
#' inf <- plague_infectious_period$relfreq
#' n <- length(inf)
#'
#' x <- seq(0, m+n+1, by = 0.02)
#' fx <- dgi(x, lat, inf)
#' Fx <- pgi(x, lat, inf)
#'
#' plot(x, fx, cex = 0.1, las = 1,
#'   xlab = "number of days",
#'   ylab = "density function"
#' )
#' plot(x, Fx, type = "l", las = 1,
#'   xlab = "number of days",
#'   ylab = "distribution function"
#' )
#'
#' @name generation-interval
NULL

#' @rdname generation-interval
#' @export
dgi <- function(x, lat, inf) {
  check(x,
    what = "numeric",
    len = c(1, Inf),
    "`x` must be numeric and have nonzero length."
  )
  for (a in c("lat", "inf")) {
    a_val <- get(a, inherits = FALSE)
    check(a_val,
      what = "numeric",
      len = c(1, Inf),
      sprintf("`%s` must be numeric and have nonzero length.", a)
    )
    check(a_val,
      val = c(0, Inf),
      yes = function(x) all(is.finite(x)),
      sprintf("`%s` must not contain missing, infinite, or negative values.", a)
    )
    check(a_val,
      no = function(x) all(x == 0),
      sprintf("`%s` must have at least one positive element.", a)
    )
  }

  if (length(x) > 1) {
    return(sapply(x, dgi, lat = lat, inf = inf))
  }

  lat <- lat / sum(lat)
  inf <- inf / sum(inf)
  m <- length(lat)
  n <- length(inf)
  re <- 1 / sum((1:n) * inf)

  if (is.na(x)) {
    NA
  } else if (x < 1) {
    0
  } else if (x < 2) {
    re * lat[1]
  } else if (x < m + 1) {
    s <- 0
    for (i in 1:(floor(x)-1)) {
      s <- s + lat[i] * (1 - sum(inf[1:min(n, floor(x)-i)]))
    }
    re * (lat[floor(x)] + s)
  } else if (x < m + n) {
    s <- 0
    for (i in 1:m) {
      s <- s + lat[i] * (1 - sum(inf[1:min(n, floor(x)-i)]))
    }
    re * s
  } else {
    0
  }
}

#' @rdname generation-interval
#' @export
pgi <- function(x, lat, inf) {
  check(x,
    what = "numeric",
    len = c(1, Inf),
    "`x` must be numeric and have nonzero length."
  )

  if (length(x) > 1) {
    return(sapply(x, pgi, lat = lat, inf = inf))
  }

  m <- length(lat)
  n <- length(inf)

  if (is.na(x)) {
    NA
  } else if (x <= 1) {
    0
  } else if (x <= 2) {
    (x - 1) * dgi(1.5, lat, inf)
  } else if (x <= m + n) {
    i <- ceiling(x) - 1
    j <- 1:(i-1)
    (x - i) * dgi(i + 0.5, lat, inf) + sum(dgi(j + 0.5, lat, inf))
  } else {
    1
  }
}
