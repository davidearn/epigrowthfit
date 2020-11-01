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
#'   in an epidemic wave.
#' @param cases A numeric vector of length `length(date)-1`.
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
#'     \item{`r`}{Initial exponential growth rate, expressed per day.}
#'     \item{`c0`}{Expected cumulative incidence on `date[first]`.
#'       Here, "cumulative incidence" refers to the number of cases
#'       observed since the start of the epidemic wave (i.e., since
#'       `date[first]`), which is not necessarily the number of
#'       cases observed since the start of the epidemic (i.e., since
#'       `date[1]`). Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{Expected epidemic final size. This is the expected
#'       number of cases observed over the full course of the epidemic
#'       wave. Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{Time at which the epidemic wave is expected
#'       to attain half its final size, expressed as a number of days
#'       since `date[first]`.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`p`}{Richards shape parameter.
#'       Used only if `curve = "richards"`.
#'     }
#'     \item{`b`}{Baseline linear growth rate, expressed per day.
#'       Used only if `include_baseline = TRUE`.
#'     }
#'     \item{`nbdisp`}{Negative binomial dispersion parameter.
#'       Used only if `distr = "nbinom"`.
#'     }
#'   }
#'
#'   `theta_init` can be `NULL` or a vector specifying a subset
#'   of the relevant parameters. Unspecified parameters are
#'   set internally. See Value.
#' @param peak,first,last Integers in `seq_along(cases)` indexing
#'   a peak in `cases` (`peak`) and endpoints of the fitting window
#'   (`first`, `last`). Alternatively, Dates between `min(date)`
#'   and `max(date)` or characters coercible to such a Date via
#'   `as.Date(x)` (e.g., "YYYY-MM-DD"). In this case, coercion
#'   from Date to index is done by `which.min(abs(date[-1] - x))`.
#'   If `NULL` (default), then an index is chosen internally.
#'   The default behaviour generates reasonable fitting windows
#'   only in special cases, and so should be accepted with caution.
#'   See Details 2.
#' @param min_first,max_first Bounds on `first` used when `first`
#'   is chosen internally. Like `first`, can be integer, Date,
#'   or character. If `NULL` (default), then the least strict
#'   bounds are used internally. See Details 2.
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
#'   \item{`first`, `last`}{Integers in `seq_along(cases)` such
#'     that the first and last observations in the fitting window
#'     are `cases[first]` and `cases[last]`. (This means that the
#'     fitting window includes cases observed between `date[first]`
#'     and `date[last+1]`.)
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
#'         `log(cumsum(cases[first:last]) + 1)`, where
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
#'   \item{`log_theta_init`}{Log-transformed `theta_init`. Identical
#'     to `log(theta_init)`, but with `"log_"` prepended to the names.
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
#'     to be updated using [stats::update()].
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
#' (which should be preferred in practice over `curve = "exponential"`)
#' as simple as possible by reducing the need for user input.
#'
#' If `peak = NULL`, then `peak` is set to `which.max(cases)`
#' internally. This default should be accepted only if `cases`
#' gives data for exactly one epidemic wave and only if `cases`
#' is close to smooth. Otherwise, `peak` should be set explicitly.
#' [smooth_cases()] can be used to find candidate values for `peak`.
#'
#' If `last = NULL`, then `last` is set to
#' `max(length(cases), peak+1)` internally.
#' If `peak < npar-1`, where `npar` is the number of
#' model parameters, then `peak` is replaced with
#' `npar-1+which.max(cases[npar:length(cases)])`
#' for the purpose of defining `last`, ensuring that
#' `last` is at least `npar`. If `peak` is specified
#' appropriately, then this default ensures that the
#' fitting window ends at the peak in the focal wave.
#'
#' If `first = NULL`, then `first` is set to
#' `max_first-which.min(cases[max_first:min_first])+1` internally.
#' This is the greatest index `i` in `min_first:max_first` satisfying
#' `cases[i] = min(cases[min_first:max_first])`. This default should
#' be accepted only if `min_first` points to the peak of the wave
#' before the focal wave (or to the start of the epidemic, if the
#' focal wave is the first wave) and only if `cases` is close to
#' smooth. If these conditions hold, then the chosen `first` will
#' point to the "base" of the focal wave, an approximation of the
#' time when `cases` begins to grow exponentially. Otherwise,
#' `first` should be set explicitly. [smooth_cases()] can
#' (and should) be used to find candidate values for `first`.
#'
#' If `min_first = NULL`, then `min_wlen` is set to 1 internally.
#' If `max_first = NULL`, then `max_wlen` is set to `last-npar+1`
#' internally. These defaults allow the fitting window to contain
#' anywhere between `npar` and `last` observations. The default
#' for `min_first` should be accepted only if the focal wave is
#' the first wave. Otherwise, it should be set explicitly to point
#' to the peak of the wave before the focal wave.
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' sc <- smooth_cases(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   log = TRUE,
#'   spar = seq(0.55, 0.9, by = 0.05)
#' )
#' plot(sc)
#' x <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   curve = "logistic",
#'   distr = "nbinom",
#'   peak = 57 # based on `spar = 0.7`
#' )
#' print(x)
#' coef(x, log = FALSE)
#' coef(x, log = TRUE)
#' predict(x)
#' plot(x, inc = "interval")
#' plot(x, inc = "cumulative")
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [smooth_cases()], [egf()], [plot.egf_init()]
#' @export
#' @import stats
egf_init <- function(date,
                     cases,
                     curve = "logistic",
                     distr = "nbinom",
                     include_baseline = FALSE,
                     theta_init = NULL,
                     peak = NULL,
                     last = NULL,
                     first = NULL,
                     min_first = NULL,
                     max_first = NULL) {
  ### VALIDATE MODEL ###################################################

  check(curve,
    what = "character",
    len = 1,
    opt = c("exponential", "logistic", "richards"),
    "`curve` must be one of ",
    "\"exponential\", \"logistic\", \"richards\"."
  )
  check(distr,
    what = "character",
    len = 1,
    opt = c("nbinom", "pois"),
    "`distr` must be one of \"nbinom\", \"pois\"."
  )
  check(include_baseline,
    what = "logical",
    len = 1,
    opt = c(TRUE, FALSE),
    "`include_baseline` must be TRUE or FALSE."
  )


  ### VALIDATE DATA ####################################################

  npar <- switch(curve, exponential = 2, logistic = 3, richards = 4) +
    (distr == "nbinom") + include_baseline

  check(date,
    what = "Date",
    len = c(npar + 1, Inf),
    "`date` must be of class \"Date\" and have ",
    "length ", npar + 1, "or greater."
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


  ### SELECT/VALIDATE FITTING WINDOW ###################################

  ## Converts numeric/Date/character value of `peak`, etc.
  ## to an integer indexing `cases`
  ndc_to_index <- function(name) {
    x <- get(name, envir = parent.frame(), inherits = FALSE)
    check(x,
      what = c("numeric", "Date", "character"),
      len = 1,
      "`", name, "` must have class ",
      "\"numeric\", \"Date\", or \"character\" and length 1."
    )
    if (is.numeric(x)) {
      check(x,
        val = c(1, length(cases)),
        "Numeric `", name, "` must not be less than 1 ",
        "or greater than `length(cases)`."
      )
      x <- as.integer(x)
    } else if (inherits(x, "Date")) {
      check(x,
        opt = seq(min(date), max(date), by = 1),
        "Date `", name, "` must not be earlier than `min(date)` ",
        "or later than `max(date)`."
      )
      x <- which.min(abs(date[-1] - x))
    } else if (is.character(x)) {
      x <- try(as.Date(x), silent = TRUE)
      check(x,
        not = "try-error",
        opt = seq(min(date), max(date), by = 1),
        "Character `", name, "` must be coercible ",
        "to Date between `min(date)` and `max(date)`."
      )
      x <- which.min(abs(date[-1] - x))
    }
    x
  }

  ## Default `peak` is (smallest) index of maximum of `cases`
  if (is.null(peak)) {
    peak <- which.max(cases)
  } else {
    peak <- ndc_to_index("peak")
  }

  ## Default `last` is `peak+1`, but must be careful to
  ## enforce `last >= npar` and `last <= length(cases)`
  if (is.null(last)) {
    if (peak < npar - 1) {
      peak_excl_start <- npar - 1 + which.max(cases[npar:length(cases)])
      last <- min(length(cases), peak_excl_start + 1)
    } else {
      last <- min(length(cases), peak + 1)
    }
  } else {
    last <- max(npar, ndc_to_index("last"))
  }

  ## Default `first` is (greatest) index of minimum of
  ## `cases[min_first:max_first]`, but must be careful to
  ## enforce `1 <= min_first <= max_first <= last-npar+1`
  if (is.null(first)) {
    if (is.null(min_first)) {
      min_first <- 1
    } else {
      min_first <- min(last - npar + 1, ndc_to_index("min_first"))
    }
    if (is.null(max_first)) {
      max_first <- last - npar + 1
    } else {
      max_first <- max(min_first, ndc_to_index("max_first"))
    }
    first <- max_first - which.min(cases[max_first:min_first]) + 1
  } else {
    first <- min(last - npar + 1, ndc_to_index("first"))
  }


  ### VALIDATE/SELECT INITIAL PARAMETER ESTIMATES ######################

  if (is.null(theta_init)) {
    theta_init <- numeric(0)
  } else {
    check(theta_init,
      what = "numeric",
      no = function(x) is.null(names(x)),
      "`theta_init` must be a named numeric vector or `NULL`."
    )
  }

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

  ## Dispense with unwanted elements of `theta_init`
  theta_init_strict <- theta_init[
    names(theta_init) %in% par &
    is.finite(theta_init) &
    theta_init > 0
  ]
  if (length(theta_init_strict) < length(theta_init)) {
    rm_names <- setdiff(names(theta_init), names(theta_init_strict))
    warning("Discarding user-specified elements of `theta_init`:\n",
            paste(rm_names, collapse = ", "),
            call. = FALSE)
  }
  theta_init <- theta_init_strict

  ## Parameters that `theta_init` does not specify
  par_missing <- setdiff(par, names(theta_init))

  ## Fit a linear model to log cumulative incidence
  ## in the first half of the fitting window
  m <- max(2, floor((last - first + 1) / 2))
  time <- as.numeric(date - date[1])
  lm_data <- data.frame(
    x = time[(first:last)+1] - time[first],
    y = log(cumsum(cases[first:last]) + 1)
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


