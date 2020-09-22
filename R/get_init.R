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
#'   initial parameter estimates:
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
#'   `theta0` can be `NULL` or a list specifying only a subset
#'   of parameters. In this case, the absent parameters are set
#'   internally as follows:
#'   `r <- beta1`, `x0 <- exp(beta0)`, `K <- sum(cases)`,
#'   `thalf <- times[peak]`, `p <- 1.0001`, `nbdisp <- 1`.
#'   Here, `beta0` and `beta1` are the intercept and slope
#'   of a linear fit to `log(cases)` within the first half
#'   of the fitting window, with zeros in `cases` replaced
#'   with 0.1.
#' @param r_if_leq0 A numeric scalar. Assigned to `r`
#'   if `!r %in% names(theta0)` and `beta1 <= 0`.
#' @param min_wlen An integer scalar. The minimum number
#'   of observations in the fitting window, which must be
#'   at least the number of model parameters (see argument
#'   `theta0`).
#' @param peak An integer in `seq_along(cases)` giving the
#'   index of the peak in `cases`, or `NULL`. See Details.
#' @param first An integer in `seq_along(cases)` giving the
#'   index of the first observation in the fitting window,
#'   or `NULL`. See Details.
#' @param first_level A numeric scalar or `NULL`. Can be
#'   used to define `first` if `first = NULL`. See Details.
#' @param skip_zero A logical scalar. If `TRUE`, then an
#'   attempt is made when defining `first` to ensure that
#'   `cases[first] > 0`. See Details.
#' @param debug_window A logical scalar. If `TRUE`, then
#'   `get_init()` returns early with a list specifying
#'   the chosen fitting window.
#' @param debug_lm A logical scalar. If `TRUE`, then
#'   a diagnostic plot of the linear model (see Value)
#'   is generated in addition to the usual output.
#'
#' @return
#' A list with integer elements `first`, `last`, and `wlen`,
#' such that `first:last` specifies the chosen fitting window,
#' and a list element `theta0` specifying the chosen initial
#' parameter estimates. `names(theta0)` is the subset of
#' `c(r, x0, K, thalf, p, nbdisp)` relevant for the indicated
#' `model` and `error` (see `theta0` in Arguments). Values
#' for parameters already specified in argument `theta0` are
#' retained. Those absent in argument `theta0` are defined as
#' follows:
#'
#' \describe{
#'   \item{`r`,`x0`}{The slope and intercept of a linear
#'     least squares fit to `log(cases_nz)` within the
#'     first half of the chosen fitting window, where
#'     `cases_nz = ifelse(cases == 0, 0.1, cases)`. If
#'     the slope is negative or zero, then `r` is assigned
#'     the value of argument `r_if_leq0`.
#'   }
#'   \item{`K`}{`sum(cases)`.}
#'   \item{`thalf`}{`times[peak]`.}
#'   \item{`p`}{1.0001.}
#'   \item{`nbdisp`}{1.}
#' }
#'
#' `peak` will match its value in the function call, but
#' `first` will differ according to the values of other
#' arguments (see Details). `last` will be equal to `peak`
#' if `model = "exponential"` and equal to
#' `max(length(cases), peak + 1)` otherwise.
#' Finally, `wlen = last-first+1` is simply the number of
#' observations in the fitting window.
#'
#' @details
#' Some details about the selection of fitting windows.
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' get_init(
#'   time = ontario$time,
#'   cases = ontario$new_confirmations,
#'   model = "richards",
#'   error = "nbinom"
#' )
#'
#' @references
#' \insertref{Ma+14}{epigrowthfitTMB}
#'
#' \insertref{Earn+20}{epigrowthfitTMB}
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
                     debug_window = FALSE,
                     debug_lm = FALSE) {
  ### CHECKS ON INPUT -----------------------------------------------

  if (!(is.character(model) && length(model) == 1 &&
          model %in% c("exponential", "logistic", "richards"))) {
    stop("`model` must be an element of ",
         "`c(\"exponential\", \"logistic\", \"richards\")`.")
  }
  if (!(is.character(error) && length(error) == 1 &&
          error %in% c("nbinom", "poisson"))) {
    stop("`error` must be an element of `c(\"nbinom\", \"poisson\")`.")
  }
  npar <- switch(model, exponential = 2, logistic = 3, richards = 4) +
    (error == "nbinom")
  pe <- 1 * (model %in% c("logistic", "richards"))
  if (!(is.numeric(min_wlen) && length(min_wlen) == 1 &&
          is.finite(min_wlen) && min_wlen >= npar)) {
    warning("`min_wlen` invalid, setting `min_wlen <- ", npar,
            "` (number of parameters).",
            call. = FALSE)
    min_wlen <- npar
  } else if (min_wlen %% 1 != 0) {
    warning("`min_wlen` is numeric but not integer, ",
            "setting `min_wlen <- floor(min_wlen)`.",
            call. = FALSE)
    min_wlen <- floor(min_wlen)
  }
  if (missing(time)) {
    stop("Missing argument `time`.")
  } else if (!(is.numeric(time) && length(time) >= npar))  {
    stop("`time` must be numeric and have length no less than ", npar, ".")
  } else if (!all(is.finite(time))) {
    stop("`time` must not contain missing or infinite values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }
  if (missing(cases)) {
    stop("Missing argument `cases`.")
  } else if (!(is.numeric(cases) && length(cases) == length(time))) {
    stop("`cases` must be numeric and have length equal to `length(time)`.")
  } else if (!(all(is.finite(cases)) && all(cases >= 0))) {
    stop("`cases` must not contain missing, infinite, or negative values.")
  }
  if (min_wlen > length(cases)) {
    warning("`min_wlen` exceeds `length(cases)`, ",
            "setting `min_wlen <- length(cases)`.",
            call. = FALSE)
    min_wlen <- length(cases)
  }
  min_peak <- min_wlen - pe
  if (!(is.numeric(peak) && length(peak) == 1 &&
          peak %in% min_peak:length(cases))) {
    stop("`peak` must be an element of `", min_peak, ":length(cases)`.")
  }
  last <- min(length(cases), peak + pe)
  max_first <- last - min_wlen + 1
  if (!is.null(first)) {
    if (!(is.numeric(first) && length(first) == 1 &&
            first %in% seq_len(max_first))) {
      stop("`first` must be an element of `1:", max_first, "`.")
    }
    if (!is.null(first_level)) {
      warning("`first` kept, `first_level` set to `NULL` and ignored.",
              call. = FALSE)
      first_level <- NULL
    }
  }
  if (!is.null(first_level)) {
    if (!(is.numeric(first_level) && length(first_level) == 1 &&
            is.finite(first_level))) {
      warning("`first_level` invalid, set to `NULL` and ignored.",
              call. = FALSE)
      first_level <- NULL
    }
  }
  if (!(is.logical(skip_zero) && length(skip_zero) == 1 && !is.na(skip_zero))) {
    warning("`skip_zero` invalid, set to `TRUE`.", call. = FALSE)
    skip_zero <- TRUE
  }
  if (!is.null(theta0) && !is.list(theta0)) {
    warning("`theta0` invalid, set to `NULL`.", call. = FALSE)
    theta0 <- NULL
  }
  if (!(is.numeric(r_if_leq0) && length(r_if_leq0) == 1 &&
          isTRUE(r_if_leq0 > 0))) {
    warning("`r_if_leq0` invalid, set to 0.1.", call. = FALSE)
    r_if_leq0 <- 0.1
  }
  if (!(is.logical(debug_window) && length(debug_window) == 1 &&
          !is.na(debug_window))) {
    warning("`debug_window` invalid, set to `FALSE`.", call. = FALSE)
    debug_window <- FALSE
  }
  if (!(is.logical(debug_lm) && length(debug_lm) == 1 && !is.na(debug_lm))) {
    warning("`debug_lm` invalid, set to `FALSE`.", call. = FALSE)
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
      ## Set `first` equal to one plus the index of the observation
      ## in `cases[1:(max_first-1)]` less than or equal to
      ## `first_level * max(cases)`. If there is more than one such
      ## index, then choose the greatest. If there is no such index,
      ## then set `first` equal to 1.
      is_below_level <- cases[1:(max_first-1)] < first_level * max(cases)
      first <- if (any(is_below_level)) max(which(is_below_level)) + 1 else 1
    }
  }

  ## Enforce `cases[first] > 0` if desired and possible
  if (skip_zero && first < max_first && cases[first] == 0) {
    is_nz <- cases[(first+1):max_first] > 0
    if (any(is_nz)) {
      new_first <- first + min(which(is_nz))
      warning("`cases[first] = 0` with `first = ", first,
              "`, advancing `first <- ", new_first, "`.",
              call. = FALSE)
      first <- new_first
    } else {
      warning("`cases[first] = 0` with `first = ", first,
              "`, but `first` cannot be advanced.",
              call. = FALSE)
    }
  }

  ## Fitting window length
  wlen <- last - first + 1

  if (debug_window) {
    out <- list(
      first = first,
      peak = peak,
      last = last,
      wlen = wlen,
      min_wlen = min_wlen
    )
    return(out)
  }

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

  ## Fit a linear model to `log(max(0.1, cases))` in the first
  ## half of the fitting window. This accommodates zeros in `cases`.
  ## NB: `max()` isn't actually vectorized so use `ifelse()`.
  m <- max(2, floor(wlen / 2))
  lm_data <- data.frame(time, cases = ifelse(cases == 0, 0.1, cases))
  lm_data <- lm_data[first:(first+m), ]
  lm_coef <- coef(lm(log(cases) ~ I(time - time[1]), data = lm_data))
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
    plot(log(cases) ~ I(time - time[1]),
      data = lm_data,
      xlab = "time (from start of window)",
      ylab = "log(max(0.1, cases))"
    )
    abline(lm_coef)
  }

  list(
    first = first,
    last = last,
    wlen = wlen,
    theta0 = theta0
  )
}
