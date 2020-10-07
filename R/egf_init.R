#' \loadmathjax
#' Define fitting windows and initial parameter estimates
#'
#' @description
#' Given an incidence or mortality time series, constructs an
#' "egf_init" object specifying a fitting window and initial
#' estimates of model parameters. Pains are taken to define a
#' window and estimates that are reasonable in the absence of
#' any user input. For model details, see [egf()].
#'
#' @param date A Date vector listing increasing time points.
#' @param cases A numeric vector with length `length(date)-1`.
#'   `cases[i]` is the number of cases observed between
#'   `date[i]` and `date[i+1]`. Here, "cases" can mean
#'   (reported) infections or (reported) deaths. Missing
#'   values are not tolerated.
#' @param curve One of `"exponential"`, `"logistic"`, and `"richards"`,
#'   indicating a phenomenological model for cumulative incidence.
#' @param include_baseline A logical scalar. If `TRUE`, then the
#'   cumulative incidence model will include a linear term \mjseqn{b t}.
#'   Assign `TRUE` only if `cases` counts deaths due to all causes and
#'   you want to model growth above a baseline mortality rate \mjseqn{b}.
#' @param distr One of `"pois"` and `"nbinom"`, indicating an
#'   observation model for interval incidence.
#' @param theta0 A named numeric vector specifying positive
#'   initial estimates of model parameters:
#'
#'   \describe{
#'     \item{`r`}{\mjseqn{\lbrace\,r\,\rbrace}
#'       Initial (exponential) growth rate expressed per day.
#'     }
#'     \item{`x0`}{\mjseqn{\lbrace\,x_0\,\rbrace}
#'       Expected cumulative incidence on `date[first]` (see `first` below).
#'       This is the expectation of the number of cases observed
#'       up to `date[first]`. Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{\mjseqn{\lbrace\,K\,\rbrace}
#'       Expected epidemic final size. This is the expectation of the
#'       number of cases observed over the full course of the epidemic.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{\mjseqn{\lbrace\,t_\text{half}\,\rbrace}
#'       Expected time at which the epidemic attains half its
#'       final size, expressed as a (possibly non-integer) number
#'       of days since `date[1]`.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`p`}{\mjseqn{\lbrace\,p\,\rbrace}
#'       Richards shape parameter. Used only if `curve = "richards"`.
#'     }
#'     \item{`b`}{\mjseqn{\lbrace\,b\,\rbrace}
#'       Baseline (linear) growth rate expressed per day.
#'       Used only if `include_baseline = TRUE`.
#'     }
#'     \item{`nbdisp`}{\mjseqn{\lbrace\,k\,\rbrace}
#'       Negative binomial dispersion parameter.
#'       Used only if `distr = "nbinom"`.
#'     }
#'   }
#'
#'   `theta0` can be `NULL` or a vector specifying a subset
#'   of the relevant parameters. Absent parameters are set
#'   internally. See Value.
#' @param min_wlen An integer scalar. The minimum number
#'   of observations in the fitting window. Must be at
#'   least the number of parameters being fit.
#' @param peak An integer in `seq_along(cases)` indexing
#'   the peak in `cases`. The index of the last observation
#'   in the fitting window will be `peak` or `peak+1`.
#'   The default is like `which.max(cases)` but careful
#'   to prevent `peak < min_wlen`. See Details 4.
#' @param first An integer in `seq_along(cases)` indexing
#'   the first observation in the fitting window, or `NULL`.
#'   Non-`NULL` values may not be respected, depending on
#'   other arguments. See Details 4.
#' @param first_level A numeric scalar or `NULL`. Can be
#'   used to define `first` if `first = NULL`. See Details 4.
#' @param skip_zero A logical scalar. If `TRUE`, then an
#'   attempt is made when defining `first` to ensure that
#'   `cases[first] > 0`. See Details 4.
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
#'     the argument of the same name. See Details 4.
#'   }
#'   \item{`last`}{An integer indexing the end of the fitting window,
#'     so that the last observation is `cases[last]`. Equal to `peak`
#'     if `curve = "exponential"` and `max(length(cases), peak + 1)`
#'     otherwise.
#'   }
#'   \item{`curve`}{Matches argument.}
#'   \item{`include_baseline`}{Matches argument.}
#'   \item{`distr`}{Matches argument.}
#'   \item{`theta0`}{A named numeric vector whose elements
#'     are the subset of `r`, `x0`, `K`, `thalf`, `p`, `b`,
#'     and `nbdisp` relevant to `curve`, `include_baseline`,
#'     and `distr`. Values from the argument of the same
#'     name are retained if they are positive and discarded
#'     otherwise. Parameters not specified in the argument
#'     are assigned their default value below:
#'
#'     \begin{describe}{
#'       \item{`r`, `x0`}{`beta1` and `exp(beta0)`, respectively,
#'         where `beta1` and `beta0` are the slope and intercept
#'         of a linear fit to `log(cumsum(cases) + 0.1)` within
#'         the first half of the fitting window. Here, "intercept"
#'         means the value of the fit at `time[first]`.
#'       }
#'       \item{`K`}{`sum(cases)`}
#'       \item{`thalf`}{`time[peak+1]`}
#'       \item{`p`}{1}
#'       \item{`b`}{1}
#'       \item{`nbdisp`}{1}
#'     }
#'
#'   }
#'   \item{`log_theta0`}{A named numeric vector, identical to
#'     `log(theta0)` with `"log_"` prepended to the names.
#'   }
#'   \item{`call`}{The call to `egf_init()`, making the output
#'     reproducible with `eval(call)`.
#'   }
#' }
#'
#' @details
#' ## 1. Phenomenological models
#' Let \mjseqn{x(t)} be the expected number of cases observed
#' up to time \mjseqn{t} (i.e., expected cumulative incidence),
#' and let \mjseqn{x_0 = x(0) > 0}. Ignoring any baseline growth
#' (see Details 2), \mjseqn{x(t)} is modeled as an exponential,
#' logistic, or Richards curve.
#'
#' ### Exponential model
#' If \mjseqn{x(t)} follows
#'
#' \mjsdeqn{x'(t) = r x(t)\,,\qquad r > 0\,,}
#'
#' then \mjseqn{x(t)} grows exponentially as
#'
#' \mjsdeqn{x(t) = x_0 e^{r t}\,.}
#'
#' Hence the exponential model for \mjseqn{x(t)} requires
#' fitting two parameters to observed data: the exponential
#' growth rate \mjseqn{r} and initial value \mjseqn{x_0}.
#'
#' The exponential model ignores depletion of susceptible
#' individuals and implies continuous exponential growth
#' in \mjseqn{x(t)}. Hence it will only agree with epidemic
#' data during the (typically short) initial exponential
#' growth phase. \insertCite{Ma+14;textual}{epigrowthfit}
#' show that estimates of \mjseqn{r} are sensitive to the
#' choice of fitting window, and that more robust fits to
#' the data are likely to obtained with the logistic and
#' Richards models, at negligible cost.
#'
#' ### Logistic model
#' If \mjseqn{x(t)} follows
#'
#' \mjsdeqn{x'(t) = r x(t) \bigg(1 - \frac{x(t)}{K}\bigg)\,,\qquad r, K > 0\,,}
#'
#' and if \mjseqn{x_0 \in (0,K)}, then \mjseqn{x(t)} grows as
#'
#' \mjsdeqn{x(t) = \frac{K}{1 + \big(\frac{K}{x_0} - 1\big) e^{-r t}}}
#'
#' and increases to \mjseqn{K} as \mjseqn{t \to \infty}.
#' The logistic model can be reparametrized as
#'
#' \mjsdeqn{x(t) = \frac{K}{1 + e^{-r (t - t_\text{half})}}\,,}
#'
#' where \mjseqn{t_\text{half}} is the time at which
#' cumulative incidence attains half its final size,
#' satisfying \mjseqn{x(t_\text{half}) = \frac{K}{2}}.
#'
#' The reparametrized logistic model requires fitting
#' \mjseqn{r}, \mjseqn{K}, and \mjseqn{t_\text{half}}
#' to observed data.
#'
#' ### Richards
#' If \mjsdeqn{x(t)} follows
#'
#' \mjsdeqn{x'(t) = r x(t) \bigg(1 - \bigg(\frac{x(t)}{K}\bigg)^p\bigg)\,,\qquad r, K, p > 0\,,}
#'
#' and if \mjseqn{x_0 \in (0,K)}, then \mjseqn{x(t)} grows as
#'
#' \mjsdeqn{x(t) = \frac{K}{\big\lbrack 1 + \big(\big(\frac{K}{x_0}\big)^p - 1\big) e^{-r p t} \big\rbrack^{1/p}}}
#'
#' and increases to \mjseqn{K} as \mjseqn{t \to \infty}.
#' The Richards model can be reparametrized as
#'
#' \mjsdeqn{x(t) = \frac{K}{\big\lbrack 1 + (2^p - 1) * e^{-r p (t - t_\text{half}}) \big\rbrack^{1/p}}\,,}
#'
#' where, as with the logistic model, \mjseqn{t_\text{half}}
#' satisfies \mjseqn{x(t_\text{half}) = \frac{K}{2}}.
#'
#' The reparametrized logistic model requires fitting
#' \mjseqn{r}, \mjseqn{K}, \mjseqn{t_\text{half}}, and \mjseqn{p}
#' to observed data.
#'
#' ## 2. Baseline growth
#' For many historical epidemics, the observed data are counts
#' deaths due to all causes, not only the disease of interest.
#' Growth in disease mortality over time can still be understood
#' phenomenologically, provided that baseline mortality (deaths
#' unrelated to the epidemic) and disease mortality are modeled
#' separately. To account for baseline mortality, it is assumed
#' that deaths occur at a constant rate \mjseqn{b > 0} in the
#' absence of an epidemic. Then, for example, the logistic model
#' (see Details 1) becomes
#'
#' \mjsdeqn{x(t) = b t + \frac{K}{1 + \big(\frac{K}{x_0} - 1\big) e^{-r t}}\,,}
#'
#' where \mjseqn{x(t)} is to be interpreted as expected
#' cumulative mortality instead of expected cumulative incidence.
#' Accounting for baseline mortality requires that \mjseqn{b} is
#' fit in addition to the other model parameters.
#'
#' ## 3. Observation model
#' Let \mjseqn{Y(t_1,t_2)} be the number of cases observed between times
#' \mjseqn{t_1} and \mjseqn{t_2 > t_1} (i.e., observed interval incidence).
#' \mjseqn{Y(t_1,t_2)} is modeled as either a Poisson-distributed random
#' variable with mean \mjseqn{x(t_2) - x(t_1)},
#'
#' \mjsdeqn{Y(t_1,t_2) \sim \mathrm{Poisson}\big(x(t_2) - x(t_1)\big)\,,}
#'
#' or as negative binomial-distributed random variable with mean
#' \mjseqn{x(t_2) - x(t_1)} and dispersion \mjseqn{k > 0},
#'
#' \mjsdeqn{Y(t_1,t_2) \sim \mathrm{NegativeBinomial}\big(x(t_2) - x(t_1), k\big)\,.}
#'
#' The negative binomial observation model requires that \mjseqn{k} is fit
#' in addition to the other model parameters.
#'
#' ## 4. Fitting window selection
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
#' plot(x)
#' plot(x, inc = "cumulative")
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [methods for class "egf_init"][egf_init-methods], [egf()]
#' @export
#' @import stats
egf_init <- function(date,
                     cases,
                     curve = "logistic",
                     include_baseline = FALSE,
                     distr = "nbinom",
                     theta0 = NULL,
                     min_wlen = 6,
                     peak = min_wlen - 1 + which.max(cases[min_wlen:length(cases)]),
                     first = NULL,
                     first_level = NULL,
                     skip_zero = TRUE) {
  ### CHECKS ON INPUT -----------------------------------------------

  if (!is.character(curve) || length(curve) != 1 ||
        !curve %in% c("exponential", "logistic", "richards")) {
    stop("`curve` must be an element of ",
         "`c(\"exponential\", \"logistic\", \"richards\")`.")
  }
  if (!is.logical(include_baseline) || length(include_baseline) != 1 ||
        is.na(include_baseline)) {
    stop("`include_baseline` must be `TRUE` or `FALSE`.")
  }
  if (!is.character(distr) || length(distr) != 1 ||
        !distr %in% c("nbinom", "pois")) {
    stop("`distr` must be an element of `c(\"nbinom\", \"pois\")`.")
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
         "the number of parameters being fit (", npar, ").")
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
  if (!is.null(theta0) && !is.numeric(theta0)) {
    warning("Setting `theta0 = NULL`.", call. = FALSE)
    theta0 <- NULL
  }


  ### FITTING WINDOW ################################################

  if (max_first == 1) {
    first <- 1
  }
  if (is.null(first)) {
    if (is.null(first_level)) {
      ## Set `first` equal to the index of the minimum of
      ## `cases[1:max_first]`. If there is more than one
      ## such index, then choose the greatest.
      first <- max_first - which.min(cases[max_first:1]) + 1
    } else {
      ## Set `first` equal to one plus the index of the observation in
      ## `cases[1:(max_first-1)]` less than `first_level * max(cases)`.
      ## If there is more than one such index, then choose the greatest.
      ## If there is no such index, then set `first` equal to 1.
      is_below_level <- cases[1:(max_first-1)] < first_level * max(cases)
      first <- if (any(is_below_level)) 1 + max(which(is_below_level)) else 1
    }
  }

  ## Enforce `cases[first] > 0` if desired and possible
  if (cases[first] == 0 && skip_zero && first < max_first) {
    is_nz <- cases[(first+1):max_first] > 0
    if (any(is_nz)) {
      new_first <- first + min(which(is_nz))
      message("`cases[first] = 0` with `first = ", first,
              "`, setting `first <- ", new_first, "`.")
      first <- new_first
    } else {
      warning("`cases[i] = 0` for all `i` in `first:max_first`, ",
              "setting `first <- max_first`.",
              call. = FALSE)
      first <- max_first
    }
  }


  ### PARAMETER ESTIMATES ###########################################

  ## Parameters that `theta0` must specify
  par <- switch(curve,
    exponential = c("r", "x0"),
    logistic    = c("r", "K", "thalf"),
    richards    = c("r", "K", "thalf", "p")
  )
  if (include_baseline) {
    par <- c(par, "b")
  }
  if (distr == "nbinom") {
    par <- c(par, "nbdisp")
  }

  ## Initialize `theta0` if necessary
  if (is.null(theta0) || is.null(names(theta0))) {
    theta0 <- numeric(0)
  }

  ## Dispense with unwanted elements
  l <- length(theta0)
  theta0 <- theta0[names(theta0) %in% par]
  theta0 <- theta0[is.finite(theta0) & theta0 > 0]
  if (length(theta0) < l) {
    warning("Discarding improperly named or valued elements of `theta0`.")
  }

  ## Parameters that `theta0` doesn't specify
  par_missing <- setdiff(par, names(theta0))

  ## Fit a linear model to `log(cumsum(cases) + 0.1)`
  ## in the first half of the fitting window
  m <- max(2, floor((last - first + 1) / 2))
  time <- as.numeric(date - date[1])
  lm_data <- data.frame(time = time[-1] - time[first], cases)
  lm_data <- lm_data[first:(first+m), ]
  lm_coef <- coef(lm(log(cumsum(cases) + 0.1) ~ time, data = lm_data))

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


  out <- list(
    date   = date,
    time   = time,
    cases  = cases,
    first  = first,
    last   = last,
    curve  = curve,
    include_baseline = include_baseline,
    distr  = distr,
    theta0 = theta0,
    log_theta0 = setNames(log(theta0), paste0("log_", names(theta0))),
    call   = match.call()
  )
  structure(out, class = c("egf_init", "list"))
}
