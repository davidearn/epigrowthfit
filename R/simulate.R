#' Simulate models of epidemic growth
#'
#' @description
#' A method for simulating incidence models specified
#' by objects of class "egf".
#'
#' @param object
#'   An "egf" object.
#' @param nsim
#'   A positive integer specifying a number of simulations.
#' @param seed
#'   An integer specifying a seed for RNG, otherwise `NULL`.
#' @param time
#'   A numeric vector of length 2 or greater listing increasing
#'   time points. Times must be expressed as numbers of days since
#'   `with(object$init, date[first])`.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' An "egf_sim" object. A list with elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`refdate`}{Matches `with(object$init, date[first])`.}
#'   \item{`int_inc`}{
#'     A matrix with `length(time)` rows and `nsim` columns.
#'     `int_inc[1, j]` is `NA` for all `j`. For `i > 1`, `int_inc[i, j]`
#'     is the number of cases observed between `time[i]` and `time[i-1]`
#'     in simulation `j` of `nsim`, sampled from a Poisson or negative
#'     binomial distribution (depending on `object$init$distr`) with
#'     mean `diff(object$eval_cum_inc(time[c(i-1, i)]))`. If the latter,
#'     then the negative binomial dispersion parameter is taken from
#'     `object$theta_fit[["nbdisp"]]`.
#'   }
#'   \item{`cum_inc`}{
#'     A matrix with `length(time)` rows and `nsim` columns.
#'     `cum_inc[i, j]` is the number of cases observed up to `time[i]`
#'     in simulation `j`. Column vector `cum_inc[, j]` is computed as
#'     `object$eval_cum_inc(time[1]) + cumsum(c(0, int_inc[-1, j]))`.
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
                         time = with(object$init, time[first:last]),
                         ...) {
  check(time,
    what = "numeric",
    len = c(2, Inf),
    "`time` must be numeric and have length 2 or greater."
  )
  check(time,
    no = anyNA,
    "`time` must not have missing values."
  )
  check(time,
    yes = function(x) all(diff(x) > 0),
    "`time` must be increasing."
  )
  check(time,
    val = with(object$init, time[c(first, last)]),
    action = "warn",
    "There are elements of `time` outside of the fitting window."
  )
  check(nsim,
    what = "numeric",
    len = 1,
    val = c(1, Inf),
    yes = is.finite,
    "`nsim` must be a positive integer."
  )
  if (!is.null(seed)) {
    check(seed,
      what = "numeric",
      yes = function(x) is.finite(x[1]),
      "`seed` must be an integer or `NULL`."
    )
  }

  ## Predicted curves
  cum_inc_pred <- object$eval_cum_inc(time)
  int_inc_pred <- diff(cum_inc_pred)

  ## Simulated curves
  if (object$init$distr == "pois") {
    set.seed(seed)
    int_inc_sim <- replicate(nsim,
      rpois(n = int_inc_pred, lambda = int_inc_pred)
    )
  } else if (object$init$distr == "nbinom") {
    nbdisp <- object$theta_fit[["nbdisp"]]
    set.seed(seed)
    int_inc_sim <- replicate(nsim,
      rnbinom(n = int_inc_pred, mu = int_inc_pred, size = nbdisp)
    )
  }
  cum_inc_sim <- cum_inc_pred[1] + apply(rbind(0, int_inc_sim), 2, cumsum)

  out <- list(
    time = time,
    refdate = with(object$init, date[first]),
    cum_inc = cum_inc_sim,
    int_inc = rbind(NA, int_inc_sim),
    object = object
  )
  structure(out, class = c("egf_sim", "list"))
}
