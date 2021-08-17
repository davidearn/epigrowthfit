#' Prior distributions
#'
#' Functions used by \code{\link{egf}} to specify prior distributions
#' of top level nonlinear model parameters.
#'
#' @param mu
#'   Mean.
#' @param sigma
#'   Standard deviation. Must be positive.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_prior"},
#' with elements:
#' \item{family}{
#'   A \link{character} string naming a family of distributions.
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
  res <- list(family = "norm", parameters = c(mu = mu, sigma = sigma))
  class(res) <- c("egf_prior", "list")
  res
}
