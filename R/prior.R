#' Prior distributions
#'
#' Functions used by \code{\link{egf}} to specify prior distributions
#' of bottom level mixed effects model parameters.
#'
#' @param mu
#'   A \link{numeric} vector listing means.
#' @param sigma
#'   A positive \link{numeric} vector listing standard deviations.
#' @param eta
#'   A positive \link{numeric} vector listing values for the shape
#'   parameter, with \code{1} corresponding to a uniform distribution
#'   over the space of symmetric positive definite matrices with
#'   unit diagonal elements. Lesser (greater) values concentrate the
#'   probability density around such matrices whose determinant is
#'   nearer to 0 (1).
#' @param df
#'   A \link{numeric} vector listing degrees of freedom.
#'   \code{df} must be greater than \code{\link{nrow}(scale) - 1}.
#'   (If either \code{df} or \code{scale} has length greater than 1,
#'   then this condition is checked pairwise after recycling.)
#' @param scale
#'   A \link{list} of symmetric positive definite \link{numeric} matrices,
#'   or a \link{matrix} to be coerced to a list of length 1.
#' @param tol
#'   A positive number specifying a tolerance for non-positive definiteness
#'   of \code{scale}. All eigenvalues must exceed \code{-tol * rho},
#'   where \code{rho} is the spectral radius of \code{scale}.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_prior"},
#' with elements:
#' \item{family}{
#'   A \link{character} string naming a family of distributions.
#' }
#' \item{parameters}{
#'   A named \link{list} of \link{numeric} vectors specifying parameter values.
#' }
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
  class(res) <- c("egf_prior", "list")
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
  class(res) <- c("egf_prior", "list")
  res
}

#' @rdname egf_prior
#' @export
Wishart <- function(df, scale, tol = 1e-06) {
  stopifnot(
    is.numeric(tol),
    length(tol) == 1L,
    is.finite(tol),
    tol > 0
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
      isSymmetric.matrix(scale[[i]]),
      (e <- eigen(scale[[i]], symmetric = TRUE, only.values = TRUE)$values) > -tol * max(abs(e))
    )
  }
  stopifnot(
    is.numeric(df),
    length(df) > 0L,
    is.finite(df),
    rep.int(df, length(scale)) > rep.int(vapply(scale, nrow, 0L), length(df)) - 1L
  )
  scale <- lapply(scale, function(S) {
    log_sd <- 0.5 * log(diag(S))
    R <- chol(S)
    R[] <- R * rep(1 / diag(R), each = nrow(R))
    chol <- R[upper.tri(R)]
    c(log_sd, chol)
  })
  res <- list(
    family = "wishart",
    parameters = list(df = as.numeric(df), scale = unname(scale))
  )
  class(res) <- c("egf_prior", "list")
  res
}

#' @rdname egf_prior
#' @export
InverseWishart <- function(df, scale, tol = 1e-06) {
  replace(Wishart(df = df, scale = scale, tol = tol), "family", list("invwishart"))
}
