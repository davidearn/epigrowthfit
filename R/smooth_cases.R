#' Fit cubic splines to interval incidence
#'
#' @description
#' Fits cubic splines of varying smoothness to an interval incidence
#' time series using [stats::smooth.spline()], and looks for local
#' extrema by checking where the first derivative changes sign.
#' Intended for use in conjunction with [plot.smooth_cases()]
#' and [egf_init()], to make the process of selecting fitting windows
#' more robust to noise.
#'
#' @param log A logical scalar. If `TRUE`, then `log10(1+cases)` is
#'   smoothed instead of `cases`.
#' @param spar A numeric vector with elements in the interval (0,1\].
#'   Degree of smoothing increases with `spar` (smoothing parameter).
#'   See [stats::smooth.spline()].
#' @param ... Optional arguments to [stats::smooth.spline()].
#' @inheritParams egf_init
#'
#' @return
#' A "smooth_cases" object. A list containing copies of arguments
#' `date`, `cases`, `log`, and `spar`, with these additional elements:
#'
#' \describe{
#'   \item{`time`}{Time as a number of days since `date[1]`.
#'     Equal to `as.numeric(date - date[1])`.
#'   }
#'   \item{`ss`}{A list of `length(spar)` "smoothing.spline" objects,
#'     giving the list output of [stats::smooth.spline()] for each
#'     element of `spar`.
#'   }
#'   \item{`peaks`}{An integer vector indexing peaks in `cases` or
#'     `log10(1+cases)`.
#'   }
#'   \item{`troughs`}{An integer vector indexing troughs in `cases` or
#'     `log10(1+cases)`.
#'   }
#' }
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' x <- smooth_cases(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   log = TRUE,
#'   spar = seq(0.55, 0.9, by = 0.05)
#' )
#' plot(x)
#'
#' @seealso [plot.smooth_cases()]
#' @export
#' @import stats
smooth_cases <- function(date, cases, log = TRUE,
                         spar = seq(0.55, 0.9, by = 0.05), ...) {
  check(date,
    what = "Date",
    len = c(5, Inf),
    "`date` must be of class \"Date\" and have ",
    "length ", 5, " or greater."
  )
  check(date,
    no = anyNA,
    "`date` must not have missing values."
  )
  check(date,
    yes = function(x) all(diff(x) > 0),
    "`date` must be increasing."
  )
  check(cases,
    what = "numeric",
    len = length(date) - 1,
    "`cases` must be numeric and have length `length(date)-1`."
  )
  check(cases,
    val = c(0, Inf),
    yes = function(x) all(is.finite(x)),
    "`cases` must not contain missing, infinite, or negative values."
  )
  check(log,
    what = "logical",
    len = 1,
    opt = c(TRUE, FALSE),
    "`log` must be TRUE or FALSE."
  )
  check(spar,
    what = "numeric",
    val = c(0, 1),
    rel = c(">", "<="),
    "`spar` must be numeric with elements in (0,1]."
  )

  time <- as.numeric(date - date[1])
  ss <- lapply(spar, function(x) {
    smooth.spline(
      x = time[-1],
      y = if (log) log10(cases + 1) else cases,
      w = diff(time),
      spar = x,
      ...
    )
  })
  dss <- lapply(ss, function(x) predict(x, deriv = 1)$y)
  peaks <- lapply(dss, function(x) 1L + which(diff(sign(x)) < 0))
  troughs <- lapply(dss, function(x) 1L + which(diff(sign(x)) > 0))

  l <- as.list(environment())
  l <- l[c("date", "time", "cases", "log", "spar",
           "ss", "peaks", "troughs")]
  structure(l, class = c("smooth_cases", "list"))
}
