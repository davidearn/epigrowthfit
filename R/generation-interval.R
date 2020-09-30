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
#'   as numbers of days.
#' @param latent A numeric vector listing the probability
#'   \mjseqn{p_i} that the latent period is \mjseqn{i} days,
#'   for all \mjseqn{i \in \lbrace 1,\ldots,m \rbrace}.
#'   Replaced with `latent / sum(latent)` in the event that
#'   `sum(latent) != 1`.
#' @param infectious A numeric vector listing the probability
#'   \mjseqn{q_i} that the infectious period is \mjseqn{i} days,
#'   for all \mjseqn{i \in \lbrace 1,\ldots,n \rbrace}.
#'   Replaced with `infectious / sum(infectious)` in the event
#'   that `sum(infectious) != 1`.
#'
#' @return
#' A numeric vector `y` with length `length(x)`.
#' For `dgi()`, `y[i]` is the generation interval density function
#' (see \mjseqn{f_Z} in Details) evaluated at `x[i]`.
#' For `pgi()`, `y[i]` is the generation interval distribution function
#' (see \mjseqn{F_Z} in Details) evaluated at `x[i]`.
#'
#' @details
#' Let \mjseqn{X} and \mjseqn{Y} be the latent and infectious
#' periods, let \mjseqn{f_X} and \mjseqn{f_Y} be their density
#' functions, and let \mjseqn{F_X} and \mjseqn{F_Y} be their
#' distribution functions. Suppose \mjseqn{f_X} and \mjseqn{f_Y}
#' are supported on \mjseqn{\lbrace 1,\ldots,m \rbrace} days
#' and \mjseqn{\lbrace 1,\ldots,n \rbrace} days, respectively,
#' so that
#'
#' \mjsdeqn{\begin{align} f_X(t) &= \sum_{i=1}^{m} p_i \delta(t - i)\,, \cr f_Y(t) &= \sum_{i=1}^{n} q_i \delta(t - i)\,, \end{align}}
#'
#' where \mjseqn{p_i,q_i \in \lbrack 0,1 \rbrack},
#' \mjseqn{\sum_{i=1}^{m} p_i = \sum_{i=1}^{n} q_i = 1}
#' and \mjseqn{\delta} is the Dirac delta.
#'
#' Let \mjseqn{W} be the infectious waiting time and \mjseqn{f_W}
#' its density. Then the generation interval is \mjseqn{Z = X + W}.
#' and has density function \mjseqn{f_Z = f_X * f_W}, assuming
#' independence of \mjseqn{X} and \mjseqn{W}.
#'
#' Equation 5.7 of \insertCite{Sven07;textual}{epigrowthfit}
#' gives an expression for \mjseqn{f_W}, assuming that
#' infectiousness is constant over the infectious period:
#'
#' \mjsdeqn{f_W = \unicode{x1D7D9}_{(0,\infty)} \frac{1 - F_Y}{\unicode{x1D53C}\lbrack Y \rbrack}\,.}
#'
#' In this case, one can show that
#'
#' \mjsdeqn{f_Z(t) = \begin{cases} 0\,, & t \in (-\infty,1\rbrack\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \sum_{i=1}^{\min\lbrace m,\lceil t \rceil-1 \rbrace} p_i\,, & t \in (1,2)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \sum_{i=1}^{\min\lbrace m, \lceil t \rceil-1 \rbrace} p_i \big(1 - \sum_{j=1}^{\min\lbrace n,\lfloor t \rfloor-i \rbrace} q_j\big)\,, & t \in \big(\lbrack 2,m+1 \rbrack \cap \unicode{x2124}\big) \cup (m+1,\infty)\,, \cr \big(\sum_{i=1}^{n} i q_i\big)^{-1} \big(p_{\lceil t \rceil-1} + \sum_{i=1}^{\lfloor t \rfloor-1} p_i \big(1 - \sum_{j=1}^{\min\lbrace n,\lfloor t \rfloor-i \rbrace} q_j\big)\big)\,, & t \in \lbrack 2,m+1 \rbrack \setminus \unicode{x2124}\,. \end{cases}}
#'
#' Note that \mjseqn{f_Z} is supported on the interval \mjseqn{(1,m+n\rbrack}
#' and constant on the interval \mjseqn{(i,i+1)} for all integers \mjseqn{i}.
#' Hence the probability that \mjseqn{Z \in (i,i+1\rbrack} is simply
#'
#' \mjsdeqn{z_i = \int_{i}^{i+1} f_Z(t)\,\text{d}t = \int_{i}^{i+1} f_Z(i+\tfrac{1}{2})\,\text{d}t = f_Z(i+\tfrac{1}{2})\,.}
#'
#' Then the distribution function of \mjseqn{Z} is
#'
#' \mjseqn{\begin{align} F_Z(t) &= \int_{-\infty}^{t} f_Z(s)\,\text{d}s \cr &= \begin{cases} 0\,, & t \in (-\infty,1\rbrack\,, \cr (t - 1) z_1\,, & t \in (1,2\rbrack\,, \cr (t - i) z_i + \sum_{j=1}^{i-1} z_j\,, & t \in (i,i+1\rbrack\,,\quad i \in \lbrace 2,\ldots,m+n-1 \rbrace\,, \cr 1\,, & t \in (m+n,\infty)\,. \end{cases} \end{align}}
#'
#' @references
#' \insertRef{Sven07}{epigrowthfit}
#'
#' @seealso [compute_R0()] for computing the basic reproduction number
#'   implied by an initial growth rate and a generation interval distribution
#'
#' @examples
#' data(plague_latent_period)
#' latent <- plague_latent_period$relfreq
#' m <- length(latent)
#'
#' data(plague_infectious_period)
#' infectious <- plague_infectious_period$relfreq
#' n <- length(infectious)
#'
#' x <- seq(0, m+n+1, by = 0.5)
#' fx <- dgi(x, latent, infectious)
#' Fx <- pgi(x, latent, infectious)
#' plot(x, fx, las = 1, xlab = "number of days", ylab = "density function")
#' plot(x, Fx, las = 1, xlab = "number of days", ylab = "distribution function")
#'
#' @name generation-interval
NULL

#' @rdname generation-interval
#' @export
dgi <- function(x, latent, infectious) {
  if (missing(x)) {
    stop("Missing argument `x`.")
  } else if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be numeric and have nonzero length.")
  }
  if (missing(latent)) {
    stop("Missing argument `latent`.")
  } else if (!is.numeric(latent) || length(latent) == 0) {
    stop("`latent` must be numeric and have nonzero length.")
  } else if (!all(is.finite(latent)) || any(latent < 0)) {
    stop("`latent` must not contain missing, infinite, or negative values.")
  } else if (all(latent == 0)) {
    stop("`latent` must have at least one positive element.")
  }
  if (!is.numeric(infectious) || length(infectious) == 0) {
    stop("`infectious` must be numeric and have nonzero length.")
  } else if (!all(is.finite(infectious)) || any(infectious < 0)) {
    stop("`infectious` must not contain missing, infinite, or negative values.")
  } else if (all(infectious == 0)) {
    stop("`infectious` must have at least one positive element.")
  }

  latent <- latent / sum(latent)
  infectious <- infectious / sum(infectious)
  m <- length(latent)
  n <- length(infectious)
  re <- 1 / sum((1:n) * infectious)

  if (length(x) == 1) {
    if (!is.finite(x)) {
      NA
    } else if (x <= 1) {
      0
    } else if (x < 2) {
      re * sum(latent[1:min(m, ceiling(x)-1)])
    } else if (x %% 1 == 0 || x > m + 1) {
      s <- 0
      for (i in 1:min(m, ceiling(x)-1)) {
        s <- s + latent[i] * (1 - sum(infectious[1:min(n, floor(x)-i)]))
      }
      re * s
    } else {
      s <- 0
      for (i in 1:(floor(x)-1)) {
        s <- s + latent[i] * (1 - sum(infectious[1:min(n, floor(x)-i)]))
      }
      re * (latent[ceiling(x)-1] + s)
    }
  } else {
    sapply(x, dgi, latent = latent, infectious = infectious)
  }
}

#' @rdname generation-interval
#' @export
pgi <- function(x, latent, infectious) {
  if (missing(x)) {
    stop("Missing argument `x`.")
  } else if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be numeric and have nonzero length.")
  }
  if (missing(latent)) {
    stop("Missing argument `latent`.")
  } else if (!is.numeric(latent) || length(latent) == 0) {
    stop("`latent` must be numeric and have nonzero length.")
  } else if (!all(is.finite(latent)) || any(latent < 0)) {
    stop("`latent` must not contain missing, infinite, or negative values.")
  } else if (all(latent == 0)) {
    stop("`latent` must have at least one positive element.")
  }
  if (!is.numeric(infectious) || length(infectious) == 0) {
    stop("`infectious` must be numeric and have nonzero length.")
  } else if (!all(is.finite(infectious)) || any(infectious < 0)) {
    stop("`infectious` must not contain missing, infinite, or negative values.")
  } else if (all(infectious == 0)) {
    stop("`infectious` must have at least one positive element.")
  }

  m <- length(latent)
  n <- length(infectious)

  if (length(x) == 1) {
    if (!is.finite(x)) {
      NA
    } else if (x <= 1) {
      0
    } else if (x <= 2) {
      (x - 1) * dgi(1.5, latent, infectious)
    } else if (x <= m + n) {
      i <- ceiling(x) - 1
      j <- 1:(i-1)
      (x - i) * dgi(i + 0.5, latent, infectious) + sum(dgi(j + 0.5, latent, infectious))
    } else {
      1
    }
  } else {
    sapply(x, pgi, latent = latent, infectious = infectious)
  }
}
