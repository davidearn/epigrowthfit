#' Fit cubic splines to interval incidence
#'
#' @description
#' Fits a cubic spline to an interval incidence time series using
#' [stats::smooth.spline()], then looks for local extrema by checking
#' when the first derivative changes sign. Intended for use in
#' conjunction with [egf_init()] to facilitate fitting window selection.
#'
#' @param formula
#'   A formula of the form `y ~ x` used to locate an interval incidence
#'   time series in `data`. Alternatively, an "egf_init" or "egf" object.
#' @param log
#'   A logical scalar. If `TRUE`, then `log10(1+cases)` is smoothed
#'   instead of `cases`.
#' @param spar
#'   A number in the interval (0,1\]. Degree of smoothing increases
#'   with `spar` (smoothing parameter). See [stats::smooth.spline()].
#' @param ...
#'   Optional arguments to [stats::smooth.spline()].
#' @inheritParams egf_init
#'
#' @return
#' The default method returns a list of class "smooth_cases"
#' containing copies of arguments `log` and `spar` and these
#' additional elements:
#'
#' \describe{
#'   \item{`data`}{
#'     A data frame with variables
#'     `date`,
#'     `time = as.integer(date - date[1])`,
#'     `difftime = c(NA, diff(time))`,
#'     `cases`, and
#'     `log1p_cases = log10(1+cases)`.
#'     `date` and `cases` are obtained from `data` using `formula`.
#'     Names in `formula` are discarded.
#'   }
#'   \item{`ss`}{
#'     The list output of [stats::smooth.spline()].
#'   }
#'   \item{`peaks`, `troughs`}{
#'     Integer vectors indexing dates of peaks and troughs,
#'     respectively, in the cubic spline specified by `ss`.
#'     Dates are obtained as `date[peaks]` and `date[troughs]`.
#'   }
#'   \item{`call`}{
#'     The call to `smooth_cases()`, allowing the output
#'     to be updated using [stats::update()].
#'   }
#' }
#'
#' The methods for classes "egf_init" and "egf" apply the default
#' method with `formula = cases ~ date` and `data = formula$data`.
#'
#' @details
#' A cubic spline is fit to `x = as.integer(date[-1]-date[1])` and
#' `y = cases[-1]` or `y = log10(1+cases[-1])`, depending on `log`.
#' The points are weighted by `diff(date)` to account for possibly
#' unequal spacing.
#'
#' @examples
#' data(canadacovid)
#' ontario <- subset(canadacovid, province == "ON")
#' x <- smooth_cases(new_confirmed ~ date,
#'   data = ontario,
#'   log = TRUE,
#'   spar = 0.7
#' )
#' plot(x)
#' v <- c("2020-03-01", "2020-03-28", "2020-09-01", "2020-09-26")
#' dline(v, lty = 2, col = "#CCCCCC")
#'
#' @seealso [plot.smooth_cases()]
#' @export
#' @import stats
smooth_cases <- function(formula = cases ~ date, ...) {
  UseMethod("smooth_cases", formula)
}

#' @rdname smooth_cases
#' @export
#' @importFrom stats smooth.spline predict
smooth_cases.default <- function(formula = cases ~ date,
                                 data = data.frame(date, cases),
                                 date,
                                 cases,
                                 log = TRUE,
                                 spar = 0.5,
                                 dfmt = "%Y-%m-%d",
                                 ...) {
  check(formula,
    what = "formula",
    len = 3,
    yes = function(x) is.name(x[[2]]) && is.name(x[[3]]),
    "`formula` must be a formula of the form `y ~ x`."
  )
  check(data,
    what = c("data.frame", "list", "environment"),
    "`data` must be a data frame, list, or environment."
  )
  dn <- all.vars(formula[[3]])
  cn <- all.vars(formula[[2]])
  found <- c(dn, cn) %in% names(data)
  check(!found,
    no = any,
    "`formula` variables not found in `data`:\n",
    paste(c(dn, cn)[!found], collapse = ", ")
  )
  date <- data[[dn]]
  if (is.character(date)) {
    date <- try(as.Date(date, tryFormats = dfmt), silent = TRUE)
  }
  check(date,
    what = "Date",
    sprintf("`%s` must be of class \"Date\" or so coercible with\n`as.Date(%s, tryFormats = dfmt)`.", dn, dn)
  )
  check(date,
    len = c(5, Inf),
    sprintf("`%s` must have length 5 or greater.", dn)
  )
  check(date,
    no = anyNA,
    sprintf("`%s` must not have missing values.", dn)
  )
  check(date,
    yes = function(x) all(diff(x) > 0),
    sprintf("`%s` must be increasing.", dn)
  )
  cases <- data[[cn]]
  check(cases,
    what = "numeric",
    len = length(date),
    sprintf("`%s` must be numeric and have length `length(%s)`.", cn, dn)
  )
  check(cases[-1],
    val = c(0, Inf),
    rel = c(">=", "<"),
    sprintf("Elements of `%s` must be finite and non-negative.", cn)
  )
  check(cases[-1],
    no = function(x) sum(!is.na(x)) < 4,
    sprintf("`%s[-1]` must have at least 4 elements that are not `NA`.", cn)
  )
  cases[1] <- NA
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

  data <- data.frame(
    date = date,
    time = days(date, since = date[1]),
    difftime = c(NA, ddiff(date)),
    cases = cases,
    log1p_cases = log10(1 + cases)
  )
  data_na_omit <- na.omit(data)
  varname <- if (log) "log1p_cases" else "cases"
  ss <- smooth.spline(
    x = data_na_omit$time,
    y = data_na_omit[[varname]],
    w = data_na_omit$difftime,
    spar = spar,
    ...
  )
  dss <- predict(ss, data$time, deriv = 1)$y
  peaks <- 1L + which(diff(sign(dss)) < 0)
  troughs <- 1L + which(diff(sign(dss)) > 0)

  l <- list(
    data = data,
    log = log,
    spar = spar,
    ss = ss,
    peaks = peaks,
    troughs = troughs,
    call = match.call()
  )
  structure(l, class = c("smooth_cases", "list"))
}

#' @rdname smooth_cases
#' @export
smooth_cases.egf_init <- function(formula, log = TRUE, spar = 0.5, ...) {
  formula <- cases ~ date
  data <- formula$data
  NextMethod("smooth_cases", formula)
}

#' @rdname smooth_cases
#' @export
smooth_cases.egf <- function(formula, log = TRUE, spar = 0.5, ...) {
  formula <- cases ~ date
  data <- formula$data
  NextMethod("smooth_cases", formula)
}
