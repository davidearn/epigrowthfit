#' Prior distributions
#'
#' Functions used by \code{\link{egf}} to specify prior distributions
#' of nonlinear and dispersion model parameters.
#'
#' @param mu
#'   Mean.
#' @param sigma
#'   Standard deviation. Must be positive.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_prior"},
#' with elements:
#' \item{distribution}{
#'   A \link{character} string naming the prior distribution.
#' }
#' \item{parameters}{
#'   A named \link{numeric} vector listing parameter values.
#' }
#'
#' @name egf_prior
NULL

#' @rdname egf_prior
#' @export
Normal <- function(mu = 0, sigma = 1) {
  stop_if_not_number_in_interval(mu, -Inf, Inf, "()")
  stop_if_not_number_in_interval(sigma, 0, Inf, "()")
  out <- list(distribution = "Normal", parameters = c(mu = mu, sigma = sigma))
  class(out) <- c("egf_prior", "list")
  out
}
