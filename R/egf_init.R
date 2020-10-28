#' \loadmathjax
#' Define a fitting window and initial parameter estimates
#'
#' @description
#' Given an interval incidence time series, constructs an
#' "egf_init" object specifying a fitting window and initial
#' estimates of model parameters. Attempts are made to define
#' a window and estimates that are reasonable in the absence
#' of any user input, but the default behaviour is not robust,
#' and it is likely that some optional arguments must be set
#' explicitly.
#'
#' @param date A Date vector listing increasing time points.
#'   Should start at or before the date of the first observed case
#'   in an epidemic.
#' @param cases A numeric vector with length `length(date)-1`.
#'   `cases[i]` is the number of cases observed between `date[i]`
#'   and `date[i+1]`. Here, "cases" can mean infections,
#'   reported infections, or reported deaths (either disease deaths or
#'   deaths due to multiple causes including the disease of interest).
#'   Missing values are not tolerated.
#' @param curve One of `"exponential"`, `"logistic"`, and `"richards"`,
#'   indicating a model of expected cumulative incidence. See Details 1.
#' @param distr One of `"pois"` and `"nbinom"`, indicating a
#'   model of observed interval incidence. See Details 1.
#' @param include_baseline A logical scalar. If `TRUE`, then
#'   the model of expected cumulative incidence will include
#'   a linear baseline. Assign `TRUE` if `cases` counts deaths
#'   due to multiple causes and `FALSE` otherwise. See Details 1.
#' @param theta_init A named numeric vector specifying positive
#'   initial estimates of relevant model parameters:
#'
#'   \describe{
#'     \item{`r`}{\mjseqn{\lbrace\,r\,\rbrace}
#'       Initial exponential growth rate expressed per day.
#'     }
#'     \item{`c0`}{\mjseqn{\lbrace\,c_0\,\rbrace}
#'       Expected cumulative incidence on `date[first]`.
#'       Here, "cumulative incidence" refers to the number of cases
#'       observed since the start of the epidemic wave (i.e., since
#'       `date[first]`), which is not necessarily the number of
#'       cases observed since the start of the epidemic (i.e., since
#'       `date[1]`). Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{\mjseqn{\lbrace\,K\,\rbrace}
#'       Expected epidemic final size. This is the expected number
#'       of cases observed over the full course of the epidemic wave.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{\mjseqn{\lbrace\,t_\textrm{half}\,\rbrace}
#'       Time at which the epidemic wave is expected to attain
#'       half its final size, expressed as a number of days since
#'       `date[first]`.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`p`}{\mjseqn{\lbrace\,p\,\rbrace}
#'       Richards shape parameter. Used only if `curve = "richards"`.
#'     }
#'     \item{`b`}{\mjseqn{\lbrace\,b\,\rbrace}
#'       Baseline linear growth rate expressed per day.
#'       Used only if `include_baseline = TRUE`.
#'     }
#'     \item{`nbdisp`}{\mjseqn{\lbrace\,k\,\rbrace}
#'       Negative binomial dispersion parameter.
#'       Used only if `distr = "nbinom"`.
#'     }
#'   }
#'
#'   `theta_init` can be `NULL` or a vector specifying a subset
#'   of the relevant parameters. Unspecified parameters are
#'   set internally. See Value.
#' @param peak An integer in `seq_along(cases)` indexing a
#'   peak in `cases`. If `cases` describes just one epidemic
#'   wave, then the default value is typically acceptable.
#'   If `cases` spans multiple waves, then `peak` should be
#'   set explicitly, as its default value may not locate the
#'   peak in the focal wave.
#' @param last An integer in `seq_along(cases)`. The last
#'   observation in the fitting window will be `cases[last]`.
#'   See Details 2.
#' @param first An integer in `seq_along(cases)`, or `NULL`.
#'   If non-`NULL`, then the first observation in the fitting
#'   window will be `cases[first]`. If `NULL`, then `first`
#'   is chosen internally according to `wlen`, and otherwise
#'   according to `min_wlen` and `max_wlen`. See Details 2.
#' @param wlen An integer indicating the number of observations
#'   in the fitting window. Must be at least the number of
#'   parameters in the model being fit. See Details 2.
#' @param min_wlen,max_wlen Integers indicating the minimum
#'   and maximum number of observations in the fitting window.
#'   Along with `last`, these determine the range of integers
#'   considered for `first` if `first` and `wlen` are `NULL`.
#'   If `cases` describes just one epidemic wave, then the
#'   default values are typically acceptable. If `cases` spans
#'   multiple epidemic waves, then `max_wlen` should be set
#'   explicitly. See Details 2.
#'
#' @return
#' An "egf_init" object. A list containing copies of arguments
#' `date`, `cases`, `last`, `curve`, `distr`, and `include_baseline`,
#' with these additional elements:
#'
#' \describe{
#'   \item{`time`}{Time as a number of days since `date[1]`.
#'     Equal to `as.numeric(date - date[1])`.
#'   }
#'   \item{`first`}{A copy of the argument if set explicitly.
#'     Otherwise, the result of an internal selection algorithm.
#'     See Details 2.
#'   }
#'   \item{`theta_init`}{A named numeric vector whose elements are
#'     the subset of `r`, `c0`, `K`, `thalf`, `p`, `b`, and
#'     `nbdisp` relevant to `curve`, `distr`, and `include_baseline`.
#'     Values from the argument are retained if they are positive
#'     numbers and discarded otherwise. Parameters not specified
#'     in the argument are assigned their default value:
#'     \describe{
#'       \item{`r`, `c0`}{`beta1` and `exp(beta0)`, respectively,
#'         where `beta1` and `beta0` are the slope and intercept
#'         of a linear model fit to the first half of
#'         `log(cumsum(cases[first:last]) + 0.1)`, where
#'         the corresponding time points are taken to be
#'         `time[(first:last)+1] - time[first]`.
#'       }
#'       \item{`K`}{`2 * sum(cases[first:peak])`}
#'       \item{`thalf`}{`time[peak+1] - time[first]`}
#'       \item{`p`}{1}
#'       \item{`nbdisp`}{1}
#'       \item{`b`}{1}
#'     }
#'     Can be extracted with `coef(object, log = FALSE)`.
#'   }
#'   \item{`log_theta_init`}{Log-transformed `theta_init`. Identical to
#'     `log(theta_init)`, but with `"log_"` prepended to the names.
#'     Can be extracted with `coef(object, log = TRUE)`.
#'   }
#'   \item{`eval_cum_inc`}{A closure with numeric arguments
#'     `time` and `theta` (default is `theta_init`) evaluating
#'     expected cumulative incidence at `time` days using
#'     parameter vector `theta`. Elements of `theta` must be
#'     named as in `theta_init`. `predict(object, time)` wraps
#'     `eval_cum_inc(time, theta = theta_init)` and provides
#'     additional useful information.
#'   }
#'   \item{`call`}{The call to `egf_init()`, allowing the output
#'     to be updated using [`update()`][stats::update()].
#'   }
#' }
#'
#' @details
#' ## 1. Models
#'
#' A full description of the models of expected cumulative
#' incidence and observed interval incidence specified by
#' arguments `curve`, `distr`, and `include_baseline`, can
#' be found in the package vignette, accessible with
#' `vignette("epigrowthfit-vignette")`.
#'
#' ## 2. Fitting window selection
#'
#' The "fitting window" is the subset of `cases` used
#' when fitting model parameters to the data, starting
#' at index `first` and ending at index `last`. If
#' `cases` contains data from multiple epidemic waves,
#' then this subset should not include data from more
#' than one wave. Regardless of the model being fit,
#' the fitting window should start at the point in the
#' focal wave at which `cases` begins growing roughly
#' exponentially. This is precisely the point at which
#' `log(cases)` becomes roughly linear. Where the
#' fitting window should end depends on the cumulative
#' incidence curve being fit to the data.
#' If `curve = "exponential"`, then the window should
#' end at the point in the focal wave at which `cases`
#' stops growing exponentially. This is the point at
#' which `log(cases)` *stops* being linear.
#' If `curve = "logistic"` or `curve = "richards"`, then
#' the fitting window should end near the peak in `cases`
#' during the focal wave. If the wave is incomplete and
#' the peak has not yet occurred (e.g., when fitting in
#' real time, during a growing epidemic), then the window
#' should simply end at `cases[length(cases)]`.
#'
#' The default behaviour of `egf_init()` tries to make
#' usage with `curve = "logistic"` and `curve = "richards"`
#' (which should be preferred in practice over
#' `curve = "exponential"`) as simple as possible
#' by reducing the need for user input. Specifically, if
#' `last` is not set explicitly, then it is set by default
#' to `min(length(cases), peak+1)`, ensuring that the
#' fitting window ends near the peak in `cases` during
#' the focal wave. If `first` is not set explicitly but
#' `wlen` is, then `first` is set to `max(1, last-wlen+1)`.
#' It neither `first` nor `wlen` is set explicitly, then
#' `egf_init()` attempts to choose `first` such that
#' `cases[first]` occurs at the "base" of the focal wave,
#' an approximation of the point at which `cases` begins
#' growing exponentially. This selection is made as follows:
#'
#' First, constraints on `first` are determined according
#' to `last`, `min_wlen` and `max_wlen`. The length
#' `wlen = last-first+1` of the fitting window is
#' constrained to be at least `min_wlen` at most `max_wlen`.
#' This means that `first` must satisfy
#' `first >= last-max_wlen+1` and `first <= last-min_wlen+1`.
#' The default values of `min_wlen` and `max_wlen` make the
#' range of valid integers as large as possible. This is
#' typically acceptable when `cases` describes just one
#' epidemic wave. However, when `cases` spans multiple waves,
#' it is necessary to set `max_wlen` so that the minimum
#' allowed value of `first` indexes a point before exponential
#' growth begins in the focal wave, but after any previous waves.
#'
#' Second, `m = min(cases[(last-max_wlen+1):(last-min_wlen+1)])`
#' is computed. Provided `max_wlen` is set appropriately, this
#' value can be described as the minimum of `cases` in the trough
#' between the focal wave and the previous wave (or the start of
#' the time series, if the focal wave is the first or only wave).
#'
#' Finally, `first` is set to the greatest integer `i` in
#' `(last-max_wlen+1):(last-min_wlen+1)` such that `cases[i] = m`.
#'
#' If `m = 0`, which typically arises when the focal wave is the
#' first wave in the epidemic, then the selected index `first`
#' is such that `cases[first] = 0`. In this situation, `egf_init()`
#' enforces `cases[first] > 0` by assigning `first <- first + 1`,
#' provided that this assignment does not cause `first` to exceed
#' the maximum allowed value, namely `last-min_wlen+1`.
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' x <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   curve = "logistic",
#'   distr = "nbinom"
#' )
#' print(x)
#' coef(x, log = FALSE)
#' coef(x, log = TRUE)
#' time_obs <- x$time
#' time_pred <- seq(min(time_obs), max(time_obs), by = median(diff(time_obs)))
#' pred <- predict(x, time = time_pred)
#' plot(x, inc = "interval")
#' plot(x, inc = "cumulative")
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [egf()], [methods for class "egf_init"][egf_init-methods]
#' @export
#' @import stats
egf_init <- function(date,
                     cases,
                     curve = "logistic",
                     distr = "nbinom",
                     include_baseline = FALSE,
                     theta_init = NULL,
                     peak = min_wlen - 1 + which.max(cases[min_wlen:length(cases)]),
                     last = min(length(cases), peak + 1),
                     first = NULL,
                     wlen = NULL,
                     min_wlen = 6,
                     max_wlen = last) {
  ### VALIDATE INPUT #################################################

  ## FIXME: Move cumbersome checks into a `validate_args()` function
  fails <- function(x, what = NULL, len = NULL, rel = "==", opt = NULL) {
    if (what == "numeric") {
      what <- c("numeric", "integer")
    }
    passes <-
      (is.null(what) || inherits(x, what)) &&
      (is.null(len) || match.fun(rel)(length(x), len)) &&
      (is.null(opt) || x %in% opt)
    !passes
  }
  if (fails(curve, "character", 1, opt = c("exponential", "logistic", "richards"))) {
    stop("`curve` must be one of \"exponential\", \"logistic\", \"richards\".")
  }
  if (fails(distr, "character", 1, opt = c("nbinom", "pois"))) {
    stop("`distr` must be one of \"nbinom\", \"pois\".")
  }
  if (fails(include_baseline, "logical", 1, opt = c(TRUE, FALSE))) {
    stop("`include_baseline` must be one of TRUE, FALSE.")
  }
  npar <- switch(curve, exponential = 2, logistic = 3, richards = 4) +
    (distr == "nbinom") + include_baseline
  if (fails(cases, "numeric", npar, ">=")) {
    stop("`cases` must be numeric and have length ", npar, " or greater.")
  } else if (!all(is.finite(cases)) || !all(cases >= 0)) {
    stop("`cases` must not contain missing, infinite, or negative values.")
  }
  if (fails(date, "Date", length(cases) + 1)) {
    stop("`date` must be of class \"Date\" and have length `length(cases)+1`.")
  } else if (anyNA(date)) {
    stop("`date` must not contain missing values.")
  } else if (!all(diff(date) > 0)) {
    stop("`date` must be increasing.")
  }
  if (fails(peak, "numeric", 1, opt = 1:length(cases))) {
    stop("`peak` must be an element of `1:length(cases)`.")
  }
  if (fails(last, "numeric", 1, opt = npar:length(cases))) {
    stop("`last` must be an element of `", npar, ":length(cases)`.")
  }
  if (is.null(first)) {
    if (is.null(wlen)) {
      if (fails(min_wlen, "numeric", 1, opt = npar:last)) {
        stop("`min_wlen` must be an element of `", npar, ":last`.")
      }
      if (fails(max_wlen, "numeric", 1, opt = min_wlen:last)) {
        stop("`max_wlen` must be an element of `min_wlen:last`.")
      }
    } else if (fails(wlen, "numeric", 1, opt = npar:last)) {
      stop("`wlen` must be an element of `", npar, ":last`.")
    }
  } else if (fails(first, "numeric", 1, opt = 1:(last-npar+1))) {
    stop("`first` must be `NULL` or an element of `",
         "1:(last-", npar - 1, ")`.")
  }
  if (!is.null(theta_init) && (!is.numeric(theta_init) || is.null(names(theta_init)))) {
    warning("`theta_init` is not a named numeric vector, setting `theta_init = NULL`.",
            call. = FALSE)
    theta_init <- NULL
  }


  ### FITTING WINDOW #################################################

  if (is.null(first)) {
    if (is.null(wlen)) {
      min_first <- last - max_wlen + 1
      max_first <- last - min_wlen + 1
      first <- max_first - which.min(cases[max_first:min_first]) + 1
      if (cases[first] == 0 && first < max_first) {
        first <- first + 1
      }
    } else {
      first <- last - wlen + 1
    }
  }


  ### INITIAL PARAMETER ESTIMATES ####################################

  ## Parameters that `theta_init` must specify
  par <- switch(curve,
    exponential = c("r", "c0"),
    logistic    = c("r", "K", "thalf"),
    richards    = c("r", "K", "thalf", "p")
  )
  if (distr == "nbinom") {
    par <- c(par, "nbdisp")
  }
  if (include_baseline) {
    par <- c(par, "b")
  }

  ## Initialize `theta_init` if necessary
  if (is.null(theta_init)) {
    theta_init <- numeric(0)
  }

  ## Dispense with unwanted elements of `theta_init`
  l <- length(theta_init)
  theta_init <- theta_init[names(theta_init) %in% par]
  theta_init <- theta_init[is.finite(theta_init) & theta_init > 0]
  if (length(theta_init) < l) {
    warning("Discarding extraneous and improperly defined elements of `theta_init`.",
            call. = FALSE)
  }

  ## Parameters that `theta_init` does not specify
  par_missing <- setdiff(par, names(theta_init))

  ## Fit a linear model to log cumulative incidence
  ## in the first half of the fitting window
  m <- max(2, floor((last - first + 1) / 2))
  time <- as.numeric(date - date[1])
  lm_data <- data.frame(
    x = time[(first:last)+1] - time[first],
    y = log(cumsum(cases[first:last]) + 0.1)
  )
  lm_coef <- coef(lm(y ~ x, data = lm_data, subset = 1:m))

  ## Default values of all parameters
  val <- c(
    r      = lm_coef[[2]],
    c0     = exp(lm_coef[[1]]),
    K      = 2 * sum(cases[first:peak]),
    thalf  = time[peak+1] - time[first],
    p      = 1,
    nbdisp = 1,
    b      = 1
  )

  ## Take from `val` what is missing in `theta_init`
  theta_init[par_missing] <- val[par_missing]
  theta_init <- theta_init[par]

  ## Define a closure that evaluates expected
  ## cumulative incidence at desired time points
  eval_cum_inc <- function(time, theta = theta_init) {
    ## Wave baseline
    c1 <- if (first > 1) sum(cases[1:(first-1)]) else 0
    ## Cumulative incidence above wave baseline
    c2 <- eval_model(time - environment(eval_cum_inc)$time[first],
      curve = curve,
      include_baseline = include_baseline,
      theta = theta
    )
    c1 + c2
  }


  out <- list(
    date = date,
    time = time,
    cases = cases,
    first = first,
    last = last,
    curve = curve,
    distr = distr,
    include_baseline = include_baseline,
    theta_init = theta_init,
    log_theta_init = setNames(log(theta_init), paste0("log_", names(theta_init))),
    eval_cum_inc = eval_cum_inc,
    call = match.call()
  )
  structure(out, class = c("egf_init", "list"))
}
