#' Prior distributions
#'
#' Functions used by \code{\link{egf}} to specify prior distributions
#' of bottom level mixed effects model parameters.
#'
#' @param mu
#'   A \link{numeric} vector listing means.
#' @param sigma
#'   A positive \link{numeric} vector listing standard deviations.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_prior"},
#' with elements:
#' \item{family}{
#'   A \link{character} string naming a family of distributions.
#' }
#' \item{parameters}{
#'   A named \link{list} of \link{numeric} vectors listing parameter values.
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
