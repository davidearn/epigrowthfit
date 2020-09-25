#' Define fitting windows and initial parameter estimates
#'
#' @description
#' Creates a list specifying a fitting window and initial
#' parameter estimates. Pains are taken to define a window
#' and estimates that are reasonable in the absence of any
#' user input.
#'
#' @param time A numeric vector listing (increasing) time points.
#' @param cases A numeric vector with length `length(time)`.
#'   `cases[i]` is the number of cases observed at time `time[i]`.
#'   Here, "cases" can be (reported) infections or (reported) deaths.
#'   Missing values are not tolerated.
#' @param model One of `"exponential"`, `"logistic"`, and `"richards"`,
#'   indicating a phenomenological model.
#' @param error One of `"poisson"` and `"nbinom"`,
#'   indicating an observation model.
#' @param theta0 A list with numeric scalar elements specifying
#'   initial estimates of these parameters:
#'
#'   \describe{
#'     \item{`r`}{Initial growth rate expressed in reciprocal units of `time`.}
#'     \item{`x0`}{Expectation of the number of cases observed
#'       at time `time[1]`. Used only if `model = "exponential"`.
#'     }
#'     \item{`K`}{Expectation of the number of cases observed over the full
#'       course of the epidemic (i.e., the expected final epidemic size).
#'       Used only if `model %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{Expectation of the time at which the epidemic attains
#'       half of its final size. Used only if
#'       `model %in% c("logistic", "richards")`
#'     }
#'     \item{`p`}{Shape parameter. Used only if `model = "richards"`.}
#'     \item{`nbdisp`}{Dispersion parameter. Used only if `error = "nbinom"`.}
#'   }
#'
#'   `theta0` can be `NULL` or a list specifying only a subset of
#'   the relevant parameters. In this case, the absent parameters
#'   are set internally. See Value.
#' @param r_if_leq0 A numeric scalar. Assigned to `r`
#'   if `!r %in% names(theta0)` and `beta1 <= 0`.
#' @param min_wlen An integer scalar. The minimum number
#'   of observations in the fitting window. Ignored if
#'   less than the number of model parameters. See Details.
#' @param peak An integer in `seq_along(cases)` indexing
#'   the peak in `cases`. The index of the last observation
#'   in the fitting window will be `peak` or `peak+1`.
#'   See Details.
#' @param first An integer in `seq_along(cases)` indexing
#'   the first observation in the fitting window, or `NULL`.
#'   See Details.
#' @param first_level A numeric scalar or `NULL`. Can be
#'   used to define `first` if `first = NULL`. See Details.
#' @param skip_zero A logical scalar. If `TRUE`, then an
#'   attempt is made when defining `first` to ensure that
#'   `cases[first] > 0`. See Details.
#' @param debug_lm A logical scalar. If `TRUE`, then
#'   a diagnostic plot of the linear model (see Value)
#'   is generated in addition to the usual output.
#'
#' @return
#' A list with integer elements `first`, `last`, `peak`
#' and `wlen` and a list element `theta0`.
#'
#' `first:last` specifies the chosen fitting window and
#' `wlen = last-first+1` the window length. `peak` will
#' match its value in the function call. `last` will be
#' equal to  `peak` if  `model = "exponential"` and
#' equal to `max(length(cases), peak + 1)` otherwise.
#' `first` can differ from its value in the function
#' call depending on the values of other arguments (see Details).
#'
#' `theta0` specifies the chosen initial parameter estimates.
#' `names(theta0)` is the subset of
#' `c("r", "x0", "K", "thalf", "p", "nbdisp")` relevant for
#' the indicated `model` and `error`. Parameter values supplied
#' in the function call are retained. Parameters not specified
#' in the function call are defined as follows:
#'
#' \describe{
#'   \item{`r`}{The slope of a linear least squares fit to
#'     `log(cases + 0.1)` within the first half of the chosen
#'     fitting window. If the slope is negative or zero, then
#'     `r` is assigned the value of argument `r_if_leq0`.
#'   }
#'   \item{`x0`}{The exponential of the intercept of the
#'     linear fit described above.
#'   }
#'   \item{`K`}{`sum(cases)`}
#'   \item{`thalf`}{`times[peak]`}
#'   \item{`p`}{1.0001}
#'   \item{`nbdisp`}{1}
#' }
#'
#' @details
#' # Fitting window selection
#'
#' The fitting window ends at index `last = peak`, where `last = peak`
#' if `model = "exponential"` and `last = max(length(cases), peak + 1)`
#' otherwise.
#'
#' The length `wlen = last-first+1` of the fitting window is constrained
#' to be at least `m = max(min_wlen, npar)`, where `npar` is the number
#' of parameters being fit:
#' 2 for `model = "exponential"`,
#' 3 for `model = "logistic"`, and
#' 4 for `model = "richards"`, plus
#' 0 for `error = "poisson"` and
#' 1 for `error = "nbinom"`.
#' This implies two constraints on the arguments: `peak >= m`
#' if `model = "exponential"` and `peak >= m-1` otherwise,
#' and `first <= peak-m+1`. An error will be thrown
#' if arguments `peak` and `first` (unless `first = NULL`)
#' do not conform to these constraints.
#'
#' If `first = NULL` and `first_level = NULL`, then `first`
#' is set internally to the greatest index `i` in `1:(peak-m+1)`
#' for which `cases[i] = min(cases[1:(peak-m+1)])`.
#' If `first = NULL` and `first_level != NULL`, then `first`
#' is set internally to one plus the greatest index `i` in
#' `1:(peak-m)` for which `cases[i] < first_level * max(cases)`.
#' (If no such index exists, then `first` is set to 1.)
#' Finally, if `first != NULL`, then the supplied value is
#' respected (unless `skip_zeros = TRUE`) and `first_level`
#' is ignored.
#'
#' If `skip_zeros = TRUE` and `first` as defined above
#' is such that `cases[first] = 0` and `first < peak-m+1`,
#' then `first` is increased by one until `cases[first] > 0`
#' or `first = peak-m+1`.
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' get_init(
#'   time  = ontario$time,
#'   cases = ontario$new_confirmations,
#'   model = "richards",
#'   error = "nbinom"
#' )
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @import graphics
#' @import stats
#' @export
get_init <- function(time,
                     cases,
                     model = "exponential",
                     error = "poisson",
                     theta0 = NULL,
                     r_if_leq0 = 0.1,
                     min_wlen = switch(model, exponential = 2, logistic = 3, richards = 4) + (error == "nbinom"),
                     peak = min_wlen - 1 + which.max(cases[min_wlen:length(cases)]),
                     first = NULL,
                     first_level = NULL,
                     skip_zero = TRUE,
                     debug_lm = FALSE) {
  ### CHECKS ON INPUT -----------------------------------------------

  if (!is.character(model) || length(model) != 1 ||
        !model %in% c("exponential", "logistic", "richards")) {
    stop("`model` must be an element of ",
         "`c(\"exponential\", \"logistic\", \"richards\")`.")
  }
  if (!is.character(error) || length(error) != 1 ||
        !error %in% c("nbinom", "poisson")) {
    stop("`error` must be an element of `c(\"nbinom\", \"poisson\")`.")
  }
  npar <- switch(model, exponential = 2, logistic = 3, richards = 4) +
    (error == "nbinom")
  pe <- 1 * (model %in% c("logistic", "richards"))
  if (missing(time)) {
    stop("Missing argument `time`.")
  } else if (!is.numeric(time) || length(time) < npar)  {
    stop("`time` must be numeric and have length no less than ", npar, ".")
  } else if (!all(is.finite(time))) {
    stop("`time` must not contain missing or infinite values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }
  if (missing(cases)) {
    stop("Missing argument `cases`.")
  } else if (!is.numeric(cases) || length(cases) != length(time)) {
    stop("`cases` must be numeric and have length equal to `length(time)`.")
  } else if (!all(is.finite(cases)) || !all(cases >= 0)) {
    stop("`cases` must not contain missing, infinite, or negative values.")
  }
  if (!is.numeric(min_wlen) || length(min_wlen) != 1 ||
        !min_wlen %in% npar:length(cases)) {
    warning("`min_wlen` set to the number of parameters: ", npar, ".", call. = FALSE)
    min_wlen <- npar
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
    warning("`skip_zero` set to `TRUE`.", call. = FALSE)
    skip_zero <- TRUE
  }
  if (!is.null(theta0) && !is.list(theta0)) {
    warning("`theta0` set to `NULL`.", call. = FALSE)
    theta0 <- NULL
  }
  if (!is.numeric(r_if_leq0) || length(r_if_leq0) != 1 ||
        !isTRUE(r_if_leq0 > 0)) {
    warning("`r_if_leq0` set to 0.1.", call. = FALSE)
    r_if_leq0 <- 0.1
  }
  if (!is.logical(debug_lm) || length(debug_lm) != 1 || is.na(debug_lm)) {
    warning("`debug_lm` set to `FALSE`.", call. = FALSE)
    debug_lm <- FALSE
  }

  ### FITTING WINDOW ------------------------------------------------

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
      warning("`cases[first] = 0` with `first = ", first,
              "`, setting `first <- ", new_first, "`.",
              call. = FALSE)
      first <- new_first
    } else {
      warning("`cases[i] = 0` for all `i` in `first:max_first`, ",
              "setting `first <- max_first`.",
              call. = FALSE)
      first <- max_first
    }
  }

  ## Fitting window length
  wlen <- last - first + 1

  ### PARAMETER ESTIMATES -------------------------------------------

  if (is.null(theta0)) {
    theta0 <- list()
  }

  ## Parameters that `theta0` must specify
  par_needed <- switch(model,
    exponential = c("r", "x0"),
    logistic    = c("r", "thalf", "K"),
    richards    = c("r", "thalf", "K", "p")
  )
  if (error == "nbinom") {
    par_needed <- c(par_needed, "nbdisp")
  }
  ## Parameters that `theta0` doesn't specify
  par_missing <- par_needed[!par_needed %in% names(theta0)]

  ## Fit a linear model to `log(cases + 0.1)` in the first half
  ## of the fitting window. This accommodates zeros in `cases`.
  m <- max(2, floor(wlen / 2))
  lm_data <- data.frame(time, cases)[first:(first+m), ]
  lm_coef <- coef(lm(log(cases + 0.1) ~ I(time - time[1]), data = lm_data))
  if (isTRUE(lm_coef[[2]] <= 0)) {
    warning("Estimated `r <= 0`, ",
            "setting `r` equal to `r_if_leq0` (", r_if_leq0, ").",
            call. = FALSE)
  }

  ## Values for parameters that `theta0` doesn't specify
  par_vals <- list(
    r      = if (isTRUE(lm_coef[[2]] <= 0)) r_if_leq0 else lm_coef[[2]],
    K      = sum(cases),
    x0     = exp(lm_coef[[1]]),
    thalf  = time[peak],
    p      = 1.0001,
    nbdisp = 1
  )
  theta0[par_missing] <- par_vals[par_missing]
  theta0 <- theta0[par_needed]
  if (!all(is.finite(unlist(theta0)))) {
    warning("`theta0` has missing or infinite values.")
  }

  if (debug_lm) {
    plot(log(cases + 0.1) ~ I(time - time[1]),
      data = lm_data,
      xlab = "time (from start of window)",
      ylab = "log(cases + 0.1)"
    )
    abline(lm_coef)
  }

  list(
    first = first,
    peak = peak,
    last = last,
    wlen = wlen,
    theta0 = theta0
  )
}
