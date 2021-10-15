#' Prior distributions
#'
#' Functions used by \code{\link{egf}} to specify prior distributions
#' of bottom level mixed effects model parameters.
#'
#' @param mu
#'   A numeric vector listing means.
#' @param sigma
#'   A positive numeric vector listing standard deviations.
#' @param eta
#'   A positive numeric vector listing values for the shape parameter,
#'   with 1 corresponding to a uniform distribution over the space of
#'   symmetric positive definite matrices with unit diagonal elements.
#'   Lesser (greater) values concentrate the probability density around
#'   such matrices whose determinant is nearer to 0 (1).
#' @param df
#'   A numeric vector listing degrees of freedom.
#'   \code{df} must be greater than \code{nrow(scale) - 1}.
#'   (If either \code{df} or \code{scale} has length greater than 1,
#'   then this condition is checked pairwise after recycling.)
#' @param scale
#'   A list of symmetric positive definite \link{numeric} matrices,
#'   or a matrix to be placed in a list of length 1.
#' @param tol
#'   A non-negative number specifying a tolerance for non-positive definiteness
#'   of \code{scale}. All eigenvalues of \code{scale} must exceed
#'   \code{-tol * rho}, where \code{rho} is the spectral radius of \code{scale}.
#'   (However, regardless of \code{tol}, \code{\link{diag}(scale)} must
#'   be positive, as standard deviations are processed on the log scale.)
#'
#' @return
#' A list inheriting from class \code{"egf_prior"}, with elements:
#' \item{family}{
#'   A character string naming a family of distributions.
#' }
#' \item{parameters}{
#'   A named list of numeric vectors specifying parameter values.
#' }
#'
#' @examples
#' Normal(mu = 0, sigma = 1)
#' Normal(mu = -5:5, sigma = c(0.1, 1))
#'
#' LKJ(eta = 2)
#'
#' U <- matrix(rnorm(9L), 3L, 3L)
#' UTU <- t(U) %*% U
#' UUT <- U %*% t(U)
#' Wishart(df = 6, scale = UTU)
#' InverseWishart(df = 6, scale = list(UTU, UUT))
#'
#' @name egf_prior
NULL

#' @rdname egf_prior
#' @export
Normal <- function(mu = 0, sigma = 1) {
  stopifnot(
    is.numeric(mu),
    length(mu) > 0L,
    is.finite(mu),
    is.numeric(sigma),
    length(sigma) > 0L,
    is.finite(sigma),
    sigma > 0
  )
  res <- list(
    family = "norm",
    parameters = list(
      mu = as.numeric(mu),
      sigma = as.numeric(sigma)
    )
  )
  class(res) <- "egf_prior"
  res
}

#' @rdname egf_prior
#' @export
LKJ <- function(eta = 1) {
  stopifnot(
    is.numeric(eta),
    length(eta) > 0L,
    is.finite(eta),
    eta > 0
  )
  res <- list(
    family = "lkj",
    parameters = list(eta = as.numeric(eta))
  )
  class(res) <- "egf_prior"
  res
}

#' @rdname egf_prior
#' @export
Wishart <- function(df, scale, tol = 1e-06) {
  stopifnot(
    is.numeric(tol),
    length(tol) == 1L,
    is.finite(tol),
    tol >= 0
  )
  if (is.matrix(scale)) {
    scale <- list(scale)
  } else {
    stopifnot(is.list(scale))
  }
  for (i in seq_along(scale)) {
    stopifnot(
      is.numeric(scale[[i]]),
      length(scale[[i]]) > 0L,
      is.finite(scale[[i]]),
      isSymmetric.matrix(scale[[i]]),
      (e <- eigen(scale[[i]], symmetric = TRUE, only.values = TRUE)$values) > -tol * abs(e[1L]),
      diag(scale[[i]]) > 0
    )
  }
  stopifnot(
    is.numeric(df),
    length(df) > 0L,
    is.finite(df),
    rep.int(df, length(scale)) > rep.int(vapply(scale, nrow, 0L), length(df)) - 1L
  )
  scale[] <- lapply(scale, cov2theta)
  res <- list(
    family = "wishart",
    parameters = list(df = as.numeric(df), scale = unname(scale))
  )
  class(res) <- "egf_prior"
  res
}

#' @rdname egf_prior
#' @export
InverseWishart <- function(df, scale, tol = 1e-06) {
  res <- Wishart(df = df, scale = scale, tol = tol)
  res$family <- "invwishart"
  res
}
