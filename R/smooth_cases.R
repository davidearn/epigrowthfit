#' Fit cubic splines to interval incidence
#'
#' @description
#' Fits a cubic spline to an interval incidence time series
#' using [stats::smooth.spline()], then looks for local extrema
#' by checking when the first derivative changes sign. Intended
#' for use in conjunction with [egf_init()] to facilitate
#' fitting window selection.
#'
#' @param object
#'   An "egf_init" or "egf" object, in which case `date` and `cases`
#'   are ignored, or `NULL`, in which case `date` and `cases` are
#'   required.
#' @param date
#'   A Date vector of length 5 or greater listing increasing time
#'   points. Should start at or before the date of the first observed
#'   case in an epidemic wave.
#' @param cases
#'   A numeric vector of length `length(date)`. For `i > 1`, `cases[i]`
#'   must specify the number of cases observed between `date[i-1]` and
#'   `date[i]`. `cases[1]` is ignored. Missing values are tolerated
#'   unless `sum(!is.na(cases[-1])) < 4`.
#' @param log
#'   A logical scalar. If `TRUE`, then `log10(1+cases)` is smoothed
#'   instead of `cases`.
#' @param spar
#'   A number in the interval (0,1\]. Degree of smoothing increases
#'   with `spar` (smoothing parameter). See [stats::smooth.spline()].
#' @param ...
#'   Optional arguments to [stats::smooth.spline()].
#'
#' @return
#' The default method returns a "smooth_cases" object, which is
#' a list containing copies of arguments `date`, `cases`, `log`,
#' and `spar`, with these additional elements:
#'
#' \describe{
#'   \item{`ss`}{
#'     The list output of [stats::smooth.spline()].
#'   }
#'   \item{`peaks`, `troughs`}{
#'     Integer vectors indexing dates of peaks and troughs,
#'     respectively, in the cubic spline specified by `ss`.
#'     Dates are obtained as `date[peaks]` and `date[troughs]`.
#'   }
#' }
#'
#' The method for class "egf_init" applies the default method
#' to `date = object$date` and `cases = object$cases`.
#'
#' The method for class "egf" applies the default method
#' to `date = object$init$date` and `cases = object$init$cases`.
#'
#' @details
#' A cubic spline is fit to `x = as.numeric(date[-1]-date[1])` and
#' `y = cases[-1]` or `y = log(1+cases[-1])`. The points are weighted
#' by `diff(date)` to account for possibly unequal spacing.
#'
#' @examples
#' data(canadacovid)
#' ontario <- subset(canadacovid, province == "ON")
#' x <- smooth_cases(
#'   date = ontario$date,
#'   cases = ontario$new_confirmed,
#'   log = TRUE,
#'   spar = 0.7
#' )
#' plot(x)
#'
#' @seealso [plot.smooth_cases()]
#' @export
#' @import stats
smooth_cases <- function(object = NULL, date, cases,
                         log = TRUE, spar = 0.5, ...) {
  check(log,
    what = "logical",
    len = 1,
    no = is.na,
    "`log` must be TRUE or FALSE."
  )
  check(spar,
    what = "numeric",
    len = 1,
    val = c(0, 1),
    rel = c(">", "<="),
    "`spar` must be a number in (0,1]."
  )

  UseMethod("smooth_cases", object)
}

#' @rdname smooth_cases
#' @export
#' @importFrom stats smooth.spline predict
smooth_cases.default <- function(object, date, cases, log, spar, ...) {
  check(date,
    what = "Date",
    len = c(5, Inf),
    "`date` must be of class \"Date\" and have length 5 or greater."
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
    len = length(date),
    "`cases` must be numeric and have length `length(date)`."
  )
  check(cases[-1],
    val = c(0, Inf),
    rel = c(">=", "<"),
    "Elements of `cases` must be finite and non-negative."
  )
  check(cases[-1],
    no = function(x) sum(!is.na(x)) < 4,
    "`cases[-1]` must have at least 4 elements that are not `NA`."
  )
  cases[1] <- NA

  data <- data.frame(
    time = days(date, since = date[1]),
    cases = cases,
    log_cases = log10(1 + cases),
    weight = c(NA, ddiff(date))
  )
  data <- na.omit(data)
  varname <- if (log) "log_cases" else "cases"
  ss <- smooth.spline(
    x = data$time,
    y = data[[varname]],
    w = data$weight,
    spar = spar,
    ...
  )
  dss <- predict(ss, data$time, deriv = 1)$y
  peaks <- 1L + which(diff(sign(dss)) < 0)
  troughs <- 1L + which(diff(sign(dss)) > 0)

  l <- list(
    date = date,
    cases = cases,
    log = log,
    spar = spar,
    ss = ss,
    peaks = peaks,
    troughs = troughs
  )
  structure(l, class = c("smooth_cases", "list"))
}

#' @rdname smooth_cases
#' @export
smooth_cases.egf_init <- function(object, date, cases, log, spar, ...) {
  date <- object$date
  cases <- object$cases
  NextMethod("smooth_cases", object)
}

#' @rdname smooth_cases
#' @export
smooth_cases.egf <- function(object, date, cases, log, spar, ...) {
  date <- object$init$date
  cases <- object$init$cases
  NextMethod("smooth_cases", object)
}
