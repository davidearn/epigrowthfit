#' \loadmathjax
#' Define a fitting window and initial parameter estimates
#'
#' @description
#' Given an interval incidence time series, constructs an
#' "egf_init" object specifying a fitting window and initial
#' estimates of model parameters. Attempts are made to define
#' a window and estimates that are reasonable in the absence
#' of any user input.
#'
#' @param date A Date vector listing increasing time points.
#'   Should start at or before the date of the first observed case.
#' @param cases A numeric vector with length `length(date)-1`.
#'   `cases[i]` is the number of cases observed between `date[i]`
#'   and `date[i+1]`. Here, "cases" can mean infections, reported
#'   infections, or reported deaths (either disease deaths or deaths
#'   due to multiple causes including the disease of interest).
#'   Missing values are not tolerated.
#' @param curve One of `"exponential"`, `"logistic"`, and `"richards"`,
#'   indicating a model of expected cumulative incidence. See Details 2.
#' @param distr One of `"pois"` and `"nbinom"`, indicating a
#'   model of observed interval incidence. See Details 2.
#' @param include_baseline A logical scalar. If `TRUE`, then
#'   the model of expected cumulative incidence will include
#'   a linear baseline \mjseqn{b t}. Assign `TRUE` if `cases`
#'   counts deaths due to multiple causes. See Details 2.
#' @param theta0 A named numeric vector specifying positive
#'   initial estimates of relevant model parameters:
#'
#'   \describe{
#'     \item{`r`}{\mjseqn{\lbrace\,r\,\rbrace}
#'       Initial exponential growth rate expressed per day.
#'     }
#'     \item{`c0`}{\mjseqn{\lbrace\,c_0\,\rbrace}
#'       Expected cumulative incidence on `date[1]`. This is
#'       the expected number of cases observed up to `date[1]`.
#'       Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{\mjseqn{\lbrace\,K\,\rbrace}
#'       Expected epidemic final size. This is the expected number
#'       of cases observed over the full course of the epidemic.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{\mjseqn{\lbrace\,t_\text{half}\,\rbrace}
#'       Time at which the epidemic is expected to attain half its
#'       final size, expressed as a (possibly non-integer) number
#'       of days since `date[1]`.
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
#'   `theta0` can be `NULL` or a vector specifying a subset
#'   of the relevant parameters. Unspecified parameters are
#'   set internally. See Value.
#' @param min_wlen An integer scalar. The minimum number
#'   of observations in the fitting window. Must be at
#'   least the number of parameters being fit.
#' @param peak An integer in `seq_along(cases)` indexing
#'   the peak in `cases`. The index of the last observation
#'   in the fitting window will be `peak` or `peak+1`.
#'   The default is like `which.max(cases)` but careful
#'   to prevent `peak < min_wlen`. See Details 1.
#' @param first An integer in `seq_along(cases)` indexing
#'   the first observation in the fitting window, or `NULL`.
#'   Non-`NULL` values may not be respected, depending on
#'   other arguments. See Details 1.
#' @param first_level A numeric scalar or `NULL`. Can be
#'   used to define `first` if `first = NULL`. See Details 1.
#' @param skip_zero A logical scalar. If `TRUE`, then an
#'   attempt is made when defining `first` to ensure that
#'   `cases[first] > 0`. See Details 1.
#'
#' @return
#' An "egf_init" object. A list with elements:
#'
#' \describe{
#'   \item{`date`}{Matches argument.}
#'   \item{`time`}{Time as a number of days since `date[1]`.
#'     Equal to `as.numeric(date - date[1])`.
#'   }
#'   \item{`cases`}{Matches argument.}
#'   \item{`first`}{An integer indexing the start of the fitting window,
#'     so that the first observation is `cases[first]`. May not match
#'     the argument of the same name. See Details 1.
#'   }
#'   \item{`last`}{An integer indexing the end of the fitting window,
#'     so that the last observation is `cases[last]`. Equal to `peak`
#'     if `curve = "exponential"` and `max(length(cases), peak + 1)`
#'     otherwise.
#'   }
#'   \item{`curve`}{Matches argument.}
#'   \item{`distr`}{Matches argument.}
#'   \item{`include_baseline`}{Matches argument.}
#'   \item{`theta0`}{A named numeric vector whose elements
#'     are the subset of `r`, `c0`, `K`, `thalf`, `p`, `b`,
#'     and `nbdisp` relevant to `curve`, `include_baseline`,
#'     and `distr`. Values from the argument of the same
#'     name are retained if they are positive and discarded
#'     otherwise. Parameters not specified in the argument
#'     are assigned their default value:
#'     \describe{
#'       \item{`r`, `c0`}{`beta1` and `exp(beta0)`, respectively,
#'         where `beta1` and `beta0` are the slope and intercept
#'         of a linear model fit to `log(cumsum(cases) + 0.1)`
#'         within the first half of the fitting window.
#'       }
#'       \item{`K`}{`sum(cases)`}
#'       \item{`thalf`}{`time[peak+1]`}
#'       \item{`p`}{1}
#'       \item{`nbdisp`}{1}
#'       \item{`b`}{1}
#'     }
#'   }
#'   \item{`log_theta0`}{Log-transformed `theta0`. Identical to
#'     `log(theta0)`, but with `"log_"` prepended to the names.
#'   }
#'   \item{`cum_inc`}{A closure with numeric arguments `time`
#'     and `theta` (default is `theta0`), evaluating expected
#'     cumulative incidence at `time` days using parameter
#'     vector `theta`. Elements must be named as in `theta0`.
#'   }
#'   \item{`call`}{The call to `egf_init()`, making the output
#'     reproducible with `eval(call)`.
#'   }
#' }
#'
#' @details
#' ## 1. Fitting window selection
#'
#' The "fitting window" is the subset of `cases` used
#' when fitting model parameters to the data. It ends
#' at index `last`, where
#' `last = peak` if `curve = "exponential"` and
#' `last = max(length(cases), peak + 1)` otherwise.
#'
#' The length `wlen = last-first+1` of the fitting window
#' is constrained to be at least `min_wlen`. This implies
#' two constraints on the arguments: `peak >= min_wlen` if
#' `curve = "exponential"` and `peak >= min_wlen-1` otherwise,
#' and `first <= peak-min_wlen+1`. An error will be thrown
#' if arguments `peak` and `first` (unless `first = NULL`)
#' do not conform to these constraints.
#'
#' If `first = NULL` and `first_level = NULL`, then `first`
#' is set internally to the greatest index `i` in `1:(peak-min_wlen+1)`
#' for which `cases[i] = min(cases[1:(peak-min_wlen+1)])`.
#' If `first = NULL` and `first_level != NULL`, then `first`
#' is set internally to one plus the greatest index `i` in
#' `1:(peak-min_wlen)` for which `cases[i] < first_level * max(cases)`.
#' If no such index exists, then `first` is set to 1.
#' Finally, if `first != NULL`, then the supplied value is
#' respected (unless `skip_zeros = TRUE`) and `first_level`
#' is ignored.
#'
#' If `skip_zeros = TRUE` and `first` (as defined above) is
#' such that `cases[first] = 0` and `first < peak-min_wlen+1`,
#' then `first` is increased by one until `cases[first] > 0`
#' or `first = peak-min_wlen+1`.
#'
#' ## 2. Models
#'
#' A full description of the models of expected cumulative
#' incidence and observed interval incidence, specified by
#' `curve`, `distr`, and `include_baseline`, can be found
#' in the package vignette, accessible with
#' `vignette("epigrowthfit-vignette")`.
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' x <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   curve = "richards",
#'   distr = "nbinom"
#' )
#' print(x)
#' coef(x, log = FALSE)
#' coef(x, log = TRUE)
#' time_obs <- x$time
#' time_pred <- seq(min(time_obs), max(time_obs), by = median(diff(time_obs)))
#' predict(x, time = time_pred)
#' plot(x, inc = "cumulative")
#' plot(x, inc = "interval")
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
                     theta0 = NULL,
                     min_wlen = 6,
                     peak = min_wlen - 1 + which.max(cases[min_wlen:length(cases)]),
                     first = NULL,
                     first_level = NULL,
                     skip_zero = TRUE) {
  ### CHECKS ON INPUT ################################################

  if (!is.character(curve) || length(curve) != 1 ||
        !curve %in% c("exponential", "logistic", "richards")) {
    stop("`curve` must be an element of ",
         "`c(\"exponential\", \"logistic\", \"richards\")`.")
  }
  if (!is.character(distr) || length(distr) != 1 ||
        !distr %in% c("nbinom", "pois")) {
    stop("`distr` must be an element of `c(\"nbinom\", \"pois\")`.")
  }
  if (!is.logical(include_baseline) || length(include_baseline) != 1 ||
        is.na(include_baseline)) {
    stop("`include_baseline` must be `TRUE` or `FALSE`.")
  }
  npar <- switch(curve, exponential = 2, logistic = 3, richards = 4) +
    (distr == "nbinom") + include_baseline
  pe <- 1 * (curve %in% c("logistic", "richards"))
  if (missing(cases)) {
    stop("Missing argument `cases`.")
  } else if (!is.numeric(cases) || length(cases) < npar) {
    stop("`cases` must be numeric and have length no less than ", npar, ".")
  } else if (!all(is.finite(cases)) || !all(cases >= 0)) {
    stop("`cases` must not contain missing, infinite, or negative values.")
  }
  if (missing(date)) {
    stop("Missing argument `date`.")
  } else if (!inherits(date, "Date") || length(date) != length(cases)+1) {
    stop("`date` must be of class \"Date\" and have length `length(cases)+1`.")
  } else if (any(is.na(date))) {
    stop("`date` must not contain missing values.")
  } else if (!all(diff(date) > 0)) {
    stop("`date` must be increasing.")
  }
  if (!is.numeric(min_wlen) || length(min_wlen) != 1 ||
        !min_wlen %in% npar:length(cases)) {
    stop("`min_wlen` must be at least ",
         "the number of parameters in the model (", npar, ").")
  }
  min_peak <- min_wlen - pe
  if (!is.numeric(peak) || length(peak) != 1 ||
        !peak %in% min_peak:length(cases)) {
    stop("`peak` must be an element of `", min_peak, ":length(cases)`.")
  }
  last <- min(length(cases), peak + pe)
  max_first <- last - min_wlen + 1
  if (!is.null(first)) {
    if (!is.numeric(first) || length(first) != 1 ||
          !first %in% seq_len(max_first)) {
      stop("`first` must be `NULL` or an element of `1:", max_first, "`.")
    }
  }
  if (!is.null(first_level)) {
    if (!is.numeric(first_level) || length(first_level) != 1 ||
          !is.finite(first_level) || !is.null(first)) {
      warning("`first` is non-`NULL`. Ignoring `first_level`.", call. = FALSE)
      first_level <- NULL
    }
  }
  if (!is.logical(skip_zero) || length(skip_zero) != 1 || is.na(skip_zero)) {
    warning("Setting `skip_zero = TRUE`.", call. = FALSE)
    skip_zero <- TRUE
  }
  if (!is.null(theta0) && (!is.numeric(theta0) || is.null(names(theta0)))) {
    warning("`theta0` is not a named numeric vector, setting `theta0 = NULL`.",
            call. = FALSE)
    theta0 <- NULL
  }


  ### FITTING WINDOW #################################################

  if (is.null(first)) {
    message("Starting with `first = NULL`.")
    if (max_first == 1) {
      first <- 1
      message("Setting `first = ", first, "`.")
    }
    if (is.null(first_level)) {
      ## Set `first` equal to the index of the minimum of
      ## `cases[1:max_first]`. If there is more than one
      ## such index, then choose the greatest.
      first <- max_first - which.min(cases[max_first:1]) + 1
      message("Setting `first = ", first, "`.")
    } else {
      ## Set `first` equal to one plus the index of the observation in
      ## `cases[1:(max_first-1)]` less than `first_level * max(cases)`.
      ## If there is more than one such index, then choose the greatest.
      ## If there is no such index, then set `first` equal to 1.
      is_below_level <- cases[1:(max_first-1)] < first_level * max(cases)
      first <- if (any(is_below_level)) 1 + max(which(is_below_level)) else 1
      message("Setting `first = ", first, "`.")
    }
  }

  ## Enforce `cases[first] > 0` if desired and possible
  if (cases[first] == 0 && skip_zero && first < max_first) {
    is_nz <- cases[(first+1):max_first] > 0
    if (any(is_nz)) {
      new_first <- first + min(which(is_nz))
      message("`cases[first] = 0` with `first = ", first,
              "`, setting `first = ", new_first, "`.")
      first <- new_first
    } else {
      warning("`cases[i] = 0` for all `i` in `first:max_first`, ",
              "setting `first = max_first`.",
              call. = FALSE)
      first <- max_first
    }
  }


  ### PARAMETER ESTIMATES ############################################

  ## Parameters that `theta0` must specify
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

  ## Initialize `theta0` if necessary
  if (is.null(theta0)) {
    theta0 <- numeric(0)
  }

  ## Dispense with unwanted elements
  l <- length(theta0)
  theta0 <- theta0[names(theta0) %in% par]
  theta0 <- theta0[is.finite(theta0) & theta0 > 0]
  if (length(theta0) < l) {
    warning("Discarding improperly named or valued elements of `theta0`.",
            call. = FALSE)
  }

  ## Parameters that `theta0` does not specify
  par_missing <- setdiff(par, names(theta0))

  ## Fit a linear model to `log(cumsum(cases) + 0.1)`
  ## in the first half of the fitting window
  m <- max(2, floor((last - first + 1) / 2))
  time <- as.numeric(date - date[1])
  lm_data <- data.frame(x = time[-1], y = cumsum(cases))
  lm_data <- lm_data[first:(first+m), ]
  lm_coef <- coef(lm(log(y + 0.1) ~ x, data = lm_data))

  ## Default values of all parameters
  val <- c(
    r      = lm_coef[[2]],
    x0     = exp(lm_coef[[1]]),
    K      = sum(cases),
    thalf  = time[peak+1],
    p      = 1,
    nbdisp = 1,
    b      = 1
  )

  ## Take from `val` what is missing in `theta0`
  theta0[par_missing] <- val[par_missing]
  theta0 <- theta0[par]

  ## Define a closure that evaluates expected cumulative incidence
  ## at desired time points
  cum_inc <- function(time, theta = theta0) {
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
    theta0 = theta0,
    log_theta0 = setNames(log(theta0), paste0("log_", names(theta0))),
    cum_inc = cum_inc,
    call = match.call()
  )
  structure(out, class = c("egf_init", "list"))
}
