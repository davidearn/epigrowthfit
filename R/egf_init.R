#' Define a fitting window and initial parameter estimates
#'
#' @description
#' Defines a fitting window and initial estimates of model parameters
#' given an interval incidence time series. Used to initialize [egf()].
#'
#' @param date
#'   A Date vector listing increasing time points. Should start at
#'   or before the date of the first observed case in an epidemic wave.
#' @param cases
#'   A numeric vector of length `length(date)`. For `i > 1`, `cases[i]`
#'   must specify the number of cases observed between `date[i-1]` and
#'   `date[i]`. `cases[1]` is ignored. Missing values are tolerated but
#'   affect the behaviour of [plot.egf_init()]. See Details 4.
#' @param curve
#'   One of `"exponential"`, `"logistic"`, and `"richards"`, indicating
#'   a model of expected cumulative incidence. See Details 1.
#' @param distr
#'   One of `"pois"` and `"nbinom"`, indicating a model of observed
#'   interval incidence. See Details 1.
#' @param include_baseline
#'   A logical scalar. If `TRUE`, then the model of expected cumulative
#'   incidence will include a linear baseline. Assign `TRUE` if `cases`
#'   counts deaths due to multiple causes and `FALSE` otherwise.
#'   See Details 1.
#' @param theta_init
#'   A named numeric vector specifying positive initial estimates of
#'   relevant model parameters:
#'
#'   \describe{
#'     \item{`r`}{Initial exponential growth rate expressed per day.}
#'     \item{`c0`}{
#'       Expected cumulative incidence on `date[first]`, where
#'       "cumulative incidence" refers to the number of cases observed
#'       since `date[first]`. Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{
#'       Expected epidemic final size. This is the expected number of
#'       cases observed over the full course of the epidemic wave.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{
#'       Time at which the epidemic wave is expected to attain half its
#'       final size expressed as a number of days since `date[first]`.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`p`}{
#'       Richards shape parameter.
#'       Used only if `curve = "richards"`.
#'     }
#'     \item{`b`}{
#'       Baseline linear growth rate expressed per day.
#'       Used only if `include_baseline = TRUE`.
#'     }
#'     \item{`nbdisp`}{
#'       Negative binomial dispersion parameter.
#'       Used only if `distr = "nbinom"`.
#'     }
#'   }
#'
#'   `theta_init` can be `NULL` or a vector specifying a subset of the
#'   relevant parameters. Unspecified parameters are set internally.
#'   See Value.
#' @param peak,first,last
#'   Integers in `seq_along(date)` indexing the time of a peak in
#'   `cases` (`peak`) and endpoints of the fitting window (`first`,
#'   `last`). Alternatively, Dates between `min(date)` and `max(date)`
#'   or characters coercible to such a Date via `as.Date()` (e.g.,
#'   "YYYY-MM-DD"). In this case, coercion from Date to index is done
#'   by `which.min(abs(date - x))`. If `NULL` (default), then an index
#'   is chosen internally. See Details 2 and 3.
#' @param min_first,max_first
#'   Bounds on `first` used when `first` is chosen internally. Like
#'   `first`, can be integer, Date, or character. If `NULL` (default),
#'   then the least strict bounds are used internally. See Details 2
#'   and 3.
#'
#' @return
#' An "egf_init" object. A list containing copies of arguments
#' `date`, `cases`, `last`, `curve`, `distr`, and `include_baseline`,
#' with these additional elements:
#'
#' \describe{
#'   \item{`first`, `last`}{
#'     Integers in `seq_along(date)` such that the endpoints of the
#'     fitting window are `date[first]` and `date[last]`.
#'   }
#'   \item{`time`}{
#'     An integer vector giving time as a number of days since
#'     `date[first]`. Equal to `as.integer(date - date[first])`.
#'   }
#'   \item{`theta_init`}{
#'     A named numeric vector whose elements are the subset of `r`,
#'     `c0`, `K`, `thalf`, `p`, `b`, and `nbdisp` relevant to `curve`,
#'     `distr`, and `include_baseline`. Values from the argument are
#'     retained if they are positive numbers and discarded otherwise.
#'     Parameters not specified in the argument are assigned their
#'     default value:
#'     \describe{
#'       \item{`K`}{`2 * sum(cases[(first+1):peak])`}
#'       \item{`thalf`}{`time[peak]`}
#'       \item{`p`}{1}
#'       \item{`nbdisp`}{1}
#'       \item{`b`}{1}
#'     }
#'     Can be extracted with `coef(object, log = FALSE)`.
#'   }
#'   \item{`log_theta_init`}{
#'     Log-transformed `theta_init`. Identical to `log(theta_init)`,
#'     but with `"log_"` prepended to the names. Can be extracted
#'     with `coef(object, log = TRUE)`.
#'   }
#'   \item{`eval_cum_inc`}{
#'     A closure with numeric arguments `time` and `theta`
#'     (default is `theta_init`) evaluating expected cumulative
#'     incidence at `time` days since `date[first]` conditional
#'     on parameter vector `theta`. Elements of `theta` must be
#'     named as in `theta_init`. `predict(object, time)` can be
#'     used instead of `object$eval_cum_inc(time)`.
#'   }
#'   \item{`call`}{
#'     The call to `egf_init()`, allowing the output to be updated
#'     using [stats::update()].
#'   }
#' }
#'
#' @details
#' ## 1. Models
#'
#' A full description of the models of expected cumulative
#' incidence and observed interval incidence specified by
#' arguments `curve`, `distr`, and `include_baseline` can
#' be found in the package vignette, accessible with
#' `vignette("epigrowthfit-vignette")`.
#'
#' The number `npar` of model parameters is given by
#' `switch(curve, exponential = 2, logistic = 3, richards = 4)`,
#' plus one if `distr = "nbinom"`,
#' plus one if `include_baseline = TRUE`.
#'
#' ## 2. Fitting window selection
#'
#' The "fitting window" is the time interval with endpoints
#' `date[first]` and `date[last]`. The cases observed during
#' this interval are counted in `cases[(first+1):last]`.
#' Only these `last-first` elements of `cases` are used when
#' fitting model parameters to the data, hence `last-first`
#' is constrained to be at least `npar`.
#'
#' If `cases` spans multiple epidemic waves, then the fitting
#' window must not include data from more than one wave. The
#' fitting window should start when `cases` (restricted to the
#' focal wave), begins growing roughly exponentially (linearly
#' on a logarithmic scale). When it should end depends on the
#' model of expected cumulative incidence being fit to the data.
#' If `curve = "exponential"`, then the window should end when
#' `cases` (restricted to the focal wave) stops growing
#' exponentially. If `curve %in% c("logistic", "richards")`,
#' then it should end at the time of the peak in `cases`
#' (during the focal wave) or, if the peak has not yet occurred,
#' at `date[length(date)]`.
#'
#' ## 3. Default behaviour
#'
#' `egf_init()` tries to make fitting window selection with
#' `curve %in% c("logistic", "richards")` as simple as possible
#' by reducing the need for user input. However, behaviour with
#' `peak`, `last`, `first`, `min_first`, and `max_first` all
#' set to `NULL` (see below) relies on several assumptions.
#' It can be expected to define a reasonable fitting window if
#' (i) `cases` gives data for exactly one epidemic wave,
#' (ii) `cases` is close to equally spaced, and
#' (iii) `cases` is close to smooth,
#' but if any of these conditions fail to hold,
#' then `peak`, `last`, and `first` may need to be set explicitly.
#' [smooth_cases()] simplifies this task by fitting cubic splines
#' to the data and providing the times of peaks and troughs in
#' the fitted curves as indices of `date`. The time of the trough
#' before the focal wave and the time of the peak in the focal wave
#' can be used as (approximate) markers for when the fitting window
#' should start and end.
#'
#' If `peak = NULL`, then `peak` is set to `which.max(cases)`
#' internally.
#'
#' If `last = NULL`, then `last` is set to `peak` internally,
#' unless `peak < npar+1`,
#' in which case `last` is set to `npar+which.max(cases[-(1:npar)])`
#' in order to enforce `last >= npar+1`.
#'
#' If `first = NULL`, then `first` is set to
#' `max_first-which.min(cases[max_first:min_first])+1` internally.
#' This is the greatest index `i` in `min_first:max_first` satisfying
#' `cases[i] = min(cases[min_first:max_first])`. This default assumes
#' that `min_first` points to the peak of the wave before the focal
#' wave or to the start of the epidemic if the focal wave is the first
#' wave. In this case, the chosen index `first` will point to the
#' "base" of the focal wave, a crude approximation of the time when
#' `cases` begins to grow exponentially.
#'
#' If `min_first = NULL`, then `min_first` is set to 1 internally.
#' If `max_first = NULL`, then `max_first` is set to `last-npar`
#' internally. These are the least strict bounds on `first` in that
#' they allow the fitting window to contain as few as `npar` and as
#' many as `last` elements of `cases`. The default `min_first = 1`
#' should be accepted only if the focal wave is the first wave and
#' otherwise set explicitly (see above paragraph).
#'
#' ## 4. Missing data
#'
#' Missing values in `cases` are tolerated insofar that they do not
#' prevent fitting of a model to the data by `egf(object)` unless
#' `sum(!is.na(cases[(first+1):last])) < npar`. However, they do
#' prevent calculation of cumulative incidence since `date[1]` and
#' so cause `plot(object, inc = "cumulative")` to throw an error.
#'
#' @examples
#' data(canadacovid)
#' ontario <- subset(canadacovid, province == "ON")
#' sc <- smooth_cases(
#'   date = ontario$date,
#'   cases = ontario$new_confirmed,
#'   log = TRUE,
#'   spar = 0.7
#' )
#' plot(sc)
#' v <- c("2020-03-06", "2020-09-01", "2020-10-01")
#' dline(v, lty = 2, col = "#CCCCCC")
#' x1 <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmed,
#'   curve = "logistic",
#'   distr = "nbinom",
#'   peak = 76,
#'   first = 32
#' )
#' x2 <- update(x1, peak = 241, first = 211)
#' print(x1)
#' coef(x1, log = FALSE)
#' coef(x1, log = TRUE)
#' predict(x1)
#' plot(x1, inc = "interval")
#' plot(x2, inc = "interval", add = TRUE)
#' plot(x1, inc = "cumulative")
#' plot(x2, inc = "cumulative", add = TRUE)
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [egf()], [smooth_cases()], [plot.egf_init()]
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
    opt = c("pois", "nbinom"),
    "`distr` must be one of \"pois\", \"nbinom\"."
  )
  check(include_baseline,
    what = "logical",
    len = 1,
    no = is.na,
    "`include_baseline` must be TRUE or FALSE."
  )

  ## Model parameters
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
  npar <- length(par)


  ### VALIDATE DATA ####################################################

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
    len = length(date),
    "`cases` must be numeric and have length `length(date)`."
  )
  check(cases[-1],
    val = c(0, Inf),
    rel = c(">=", "<"),
    "Elements of `cases` must be finite and non-negative."
  )
  cases[1] <- NA


  ### SELECT/VALIDATE FITTING WINDOW ###################################

  ## Converts numeric/Date/character value of `peak`, etc.
  ## to integer indexing `date`
  ndc_to_index <- function(name) {
    x <- get(name, envir = parent.frame(), inherits = FALSE)
    check(x,
      what = c("numeric", "Date", "character"),
      len = 1,
      "`", name, "` must have class \"numeric\", \"Date\",\n",
      "or \"character\" and length 1."
    )
    if (is.numeric(x)) {
      check(x,
        opt = seq_along(date),
        "Numeric `", name, "` must be an integer in `seq_along(date)."
      )
      x <- as.integer(x)
    } else if (inherits(x, "Date")) {
      check(x,
        opt = seq(min(date), max(date), by = 1),
        "Date `", name, "` must not be earlier than `min(date)`\n",
        "or later than `max(date)`."
      )
      x <- which.min(abs(date - x))
    } else if (is.character(x)) {
      x <- try(as.Date(x), silent = TRUE)
      check(x,
        not = "try-error",
        opt = seq(min(date), max(date), by = 1),
        "Character `", name, "` must be coercible ",
        "to Date between `min(date)` and `max(date)`."
      )
      x <- which.min(abs(date - x))
    }
    x
  }

  ## Default `peak` is (smallest) index of maximum of `cases`
  if (is.null(peak)) {
    peak <- which.max(cases)
  } else {
    peak <- ndc_to_index("peak")
  }

  ## Default `last` is `peak`, but must take care to enforce
  ## `last >= npar+1`
  if (is.null(last)) {
    if (peak >= npar + 1) {
      last <- peak
    } else {
      last <- npar + which.max(cases[-(1:npar)])
    }
  } else {
    last <- max(npar + 1, ndc_to_index("last"))
  }

  ## Default `first` is (greatest) index of minimum of
  ## `cases[min_first:max_first]`, but must take care to
  ## enforce `1 <= min_first <= max_first <= last-npar`
  if (is.null(first)) {
    if (is.null(min_first)) {
      min_first <- 1
    } else {
      min_first <- min(last - npar, ndc_to_index("min_first"))
    }
    if (is.null(max_first)) {
      max_first <- last - npar
    } else {
      max_first <- max(min_first, ndc_to_index("max_first"))
    }
    first <- max_first - which.min(cases[max_first:min_first]) + 1
  } else {
    first <- min(last - npar, ndc_to_index("first"))
  }

  ## Fitting window must contain at least `npar` observations
  ## after missing values are excluded
  if (sum(!is.na(cases[(first+1):last])) < npar) {
    warning("Length of `cases[(first+1):last]`, excluding NA,\n",
            "is less than the number of model parameters.",
            call. = FALSE)
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
  time <- days(date, since = date[first])
  h <- max(2, floor((last - first) / 2))
  x <- time[first+(1:h)]
  y <- log(1 + cumsum(cases[first+(1:h)]))
  lm_coef <- try(silent = TRUE, expr = {
    coef(lm(y ~ x, data = data.frame(x, y), na.action = na.omit))
  })
  lm_fail_flag <- inherits(lm_coef, "try-error")

  ## Default values of all parameters
  val <- c(
    r      = if (lm_fail_flag) 0.1 else lm_coef[[2]],
    c0     = if (lm_fail_flag) 1 else exp(lm_coef[[1]]),
    K      = 2 * sum(cases[(first+1):peak]),
    thalf  = time[peak],
    p      = 1,
    nbdisp = 1,
    b      = 1
  )

  ## Take from `val` what is missing in `theta_init`
  theta_init[par_missing] <- val[par_missing]
  theta_init <- theta_init[par]

  ## Define a closure that evaluates expected cumulative incidence
  ## at times in days since the start of the fitting window
  eval_cum_inc <- function(time, theta = theta_init) {
    eval_model(time,
      curve = curve,
      include_baseline = include_baseline,
      theta = theta
    )
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
