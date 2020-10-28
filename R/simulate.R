#' Simulate models of epidemic growth
#'
#' @description
#' A method for simulating incidence models specified
#' by objects of class "egf".
#'
#' @param object An "egf" object.
#' @param nsim A positive integer specifying a number of simulations.
#' @param seed An integer specifying a seed for RNG, otherwise `NULL`.
#' @param time A numeric vector listing increasing time points,
#'   expressed as numbers of days since `object$date[1]`
#'   for "egf_init" objects and since `object$init$date[1]`
#'   for "egf" objects. Missing values are not tolerated and
#'   `length(time) >= 2` is required.
#' @param ... Unused optional arguments.
#'
#' @return
#' An "egf_sim" object. A list with elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`refdate`}{Matches `object$init$date[1]`.}
#'   \item{`int_inc`}{A matrix with `length(time)-1` rows and `nsim`
#'     columns, such that `int_inc[i, j]` is the number of cases
#'     observed between `time[i]` and `time[i+1]` in simulation `j`
#'     of `nsim`. Row vector `int_inc[i, ]` is sampled from a
#'     Poisson or negative binomial distribution (depending on
#'     `object$init$distr`) with mean `predict(object, time)$int_inc[i]`.
#'     The negative binomial dispersion parameter is taken from
#'     `object$theta_fit[["nbdisp"]]`.
#'   }
#'   \item{`cum_inc`}{A matrix with `length(time)` rows and `nsim`
#'     columns, such that `cum_inc[i, j]` is the number of cases
#'     observed up to `time[i]` in simulation `j`. Column vector
#'     `cum_inc[, j]` is computed as `c0 + cumsum(c(0, int_inc[, j]))`,
#'     where `c0 = predict(object, time)$cum_inc[1]` is the expected
#'     value of cumulative incidence at `time[1]` conditional on
#'     parameter vector `object$theta_fit`.
#'   }
#'   \item{`object`}{Matches argument.}
#' }
#'
#' @seealso [egf()], [plot.egf_sim()]
#' @name simulate.egf
NULL

#' @rdname simulate.egf
#' @export
#' @importFrom stats rpois rnbinom
simulate.egf <- function(object, nsim = 1, seed = NULL,
                         time = object$init$time, ...) {
  if (!is.numeric(time) || length(time) < 2) {
    stop("`time` must be numeric and have length 2 or greater.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }
  if (!is.numeric(nsim) || length(nsim) != 1 || !isTRUE(nsim >= 1)) {
    stop("`nsim` must be a positive integer.")
  }
  if (!is.null(seed) && (!is.numeric(seed) || !is.finite(seed[1]))) {
    stop("`seed` must be `NULL` or an integer.")
  }

  ## Predicted curves
  cum_inc <- object$eval_cum_inc(time)
  int_inc <- diff(cum_inc)

  ## Simulated curves
  if (object$init$distr == "pois") {
    set.seed(seed)
    sim <- replicate(nsim, rpois(int_inc, lambda = int_inc))
  } else if (object$init$distr == "nbinom") {
    nbdisp <- object$theta_fit[["nbdisp"]]
    set.seed(seed)
    sim <- replicate(nsim, rnbinom(int_inc, mu = int_inc, size = nbdisp))
  }

  out <- list(
    time = time,
    refdate = object$init$date[1],
    cum_inc = cum_inc[1] + apply(rbind(0, sim), 2, cumsum),
    int_inc = sim,
    object = object
  )
  structure(out, class = c("egf_sim", "list"))
}
