get_init <- function(time,
                     cases,
                     model = "exponential",
                     error = "poisson",
                     theta = NULL,
                     r_if_leq0 = 0.1,
                     min_wlen = 3 + (model == "richards") + (error == "nbinom"),
                     peak = min_wlen - 1 + which.max(cases[min_wlen:length(cases)]),
                     first = NULL,
                     first_level = NULL,
                     skip_zero = TRUE,
                     debug_window = FALSE,
                     debug_coef = FALSE) {
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
  npar <- 3 + (model == "richards") + (error == "nbinom")
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
  if (!(is.logical(skip_zero) && length(skip_zero) == 1 && is.finite(skip_zero))) {
    warning("`skip_zero` invalid, set to `TRUE`.",
            call. = FALSE)
    skip_zero <- TRUE
  }
  if (!is.null(theta0) && !is.list(theta0)) {
    warning("`theta0` invalid, set to `NULL`.",
            call. = FALSE)
    theta0 <- NULL
  }

  ### FITTING WINDOW ------------------------------------------------

  if (is.null(first)) {
    if (is.null(first_level)) {
      ## Set `first` equal to the index of the minimum of `cases[1:peak]`.
      ## If there is more than one such index, then choose the greatest.
      first <- peak - which.min(cases[peak:1]) + 1
    } else {
      ## Set `first` equal to the index of the observation in `cases[1:peak]`
      ## less than or equal to `first_level * max(cases)`. If there is more
      ## than one such index, then choose the greatest. If there is no such
      ## index, then set `first` equal to 1.
      is_below_level <- cases[1:peak] <= first_level * max(cases)
      first <- if (any(is_below_level)) max(which(is_below_level)) else 1
    }
  }

  ## Enforce `cases[first] > 0` if possible
  if (start_at_nz && cases[first] == 0) {
    is_nz <- cases[first:peak] > 0
    if (any(is_nz)) {
      first <- min(which(is_nz))
      warning("`cases[first] = 0`, setting `first <- ", first, "`.")
    }
  }

  ## Make sure `last - first + 1` is at least `min_wlen`
  if (first > max_first) {
    first <- min(first, max_first)
    warning("`first` exceeds maximum ", max_first,
            ", setting `first <- ", max_first, "`.")
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

  ### THETA ---------------------------------------------------------

  if (is.null(theta0)) {
    theta0 <- list()
  }

  ## Parameters needed by `theta0`
  par_needed <- switch(model,
    exponential = c("r", "K", "x0"),
    logistic    = c("r", "K", "thalf"),
    richards    = c("r", "K", "thalf", "p")
  )
  if (error == "nbinom") {
    par_needed <- c(par_needed, "nbdisp")
  }
  ## Parameters missing in `theta0`
  par_missing <- par_needed[!par_needed %in% theta0]

  ## Fit a linear model to `log(cases + 0.1)` in the first half
  ## of the fitting window, adding 0.1 to allow zeros
  if (length(par_missing) > 0) {
    m <- max(2, floor(wlen / 2))
  }
  lm_data <- data.frame(time, cases)[first+(0:m)]
  lm_coef <- coef(lm(log(cases + 0.1) ~ I(time - time[1]), data = lm_data))
  if (isTRUE(lm_coef[[2]] <= 0)) {
    warning("Estimated `r <= 0`, setting `r` equal to `r_if_leq0`.")
  }

  ## List of parameter values to add to `theta0` if necessary
  par_vals <- list(
    r = if (isTRUE(lm_coef[[2]] <= 0)) r_if_leq0 else lm_coef[[2]],
    K = sum(cases),
    x0 = exp(lm_coeff[[1]]),
    thalf = time[peak],
    p = 1.0001,
    nbdisp = 1
  )
  theta0[par_missing] <- par_vals[par_missing]
  theta0 <- theta0[par_needed]

  if (debug_coef) {
    plot(log(cases + 0.1) ~ I(time - time[1]), data = lm_data)
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
