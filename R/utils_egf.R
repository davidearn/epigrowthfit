#' Get top level nonlinear model parameter names
#'
#' Retrieves the names used internally for top level nonlinear model
#' parameter names.
#'
#' @param object
#'   An \code{"\link{egf_model}"}, \code{"\link{egf}"},
#'   or \code{"\link[=make_tmb_data]{tmb_data}"} object
#'   specifying a top level nonlinear model. Otherwise,
#'   \code{\link{NULL}}.
#' @param link
#'   A \link{logical} flag. If \code{TRUE}, then \code{"name"}
#'   is replaced with \code{"log(name)"} or \code{"logit(name)"}
#'   depending on the link function used internally for the
#'   parameter.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{character} vector giving the subset of names relevant to
#' \code{object}, or the complete set if \code{object = \link{NULL}}.
#'
#' @examples
#' get_names_top(NULL, link = FALSE)
#' get_names_top(NULL, link = TRUE)
#' get_names_top(egf_model(), link = FALSE)
#'
#' @export
get_names_top <- function(object, ...) {
  UseMethod("get_names_top", object)
}

#' @rdname get_names_top
#' @export
get_names_top.default <- function(object, link = TRUE, ...) {
  stopifnot(is.null(object))
  names_top <- c("r", "alpha", "c0", "tinfl", "K",
                 "p", "a", "b", "disp", paste0("w", 1:6))
  if (link) {
    return(string_add_link(names_top))
  }
  names_top
}

#' @rdname get_names_top
#' @export
get_names_top.egf_model <- function(object, link = TRUE, ...) {
  names_top <- switch(object$curve,
    exponential    = c("r", "c0"),
    subexponential = c("alpha", "c0", "p"),
    gompertz       = c("alpha", "tinfl", "K"),
    logistic       = c("r", "tinfl", "K"),
    richards       = c("r", "tinfl", "K", "a")
  )
  if (object$excess) {
    names_top <- c(names_top, "b")
  }
  if (object$family == "nbinom") {
    names_top <- c(names_top, "disp")
  }
  if (object$day_of_week > 0L) {
    names_top <- c(names_top, paste0("w", 1:6))
  }
  if (link) {
    return(string_add_link(names_top))
  }
  names_top
}

#' @rdname get_names_top
#' @export
get_names_top.egf <- function(object, link = TRUE, ...) {
  get_names_top(object$model, link = link)
}

#' @rdname get_names_top
#' @export
get_names_top.tmb_data <- function(object, link = TRUE, ...) {
  names_top <- levels(object$X_info$par)
  if (link) {
    return(names_top)
  }
  string_remove_link(names_top)
}

#' Manipulate top level nonlinear model parameter names
#'
#' Utilities for modifying strings \code{s} and
#' \code{fs = \link{sprintf}("\%s(\%s)", f, s)},
#' where \code{s} is the name used internally
#' for a top level nonlinear model parameter and
#' \code{f} is the name (either \code{"log"} or \code{"logit"})
#' of the corresponding link function.
#'
#' @param s,fs \link[=character]{Character} vectors.
#'
#' @return
#' For strings \code{s}, \code{f}, and \code{fs} as described above,
#' \code{string_get_link(s)} returns \code{f}, \code{string_add_link(s)}
#' returns \code{fs}, \code{string_remove_link(fs)} returns \code{s},
#' and \code{string_extract_link(fs)} returns \code{f}.
#'
#' All functions are vectorized. \code{\link{NA_character_}}
#' is returned for invalid input (see Examples).
#'
#' @examples
#' ## string_get_link("r")
#' ## string_add_link("r")
#' ## string_remove_link("log(r)")
#' ## string_extract_link("log(r)")
#' ## string_extract_link("invalid", "input", "r" , "log(r)")
#'
#' @name string_get_link
#' @keywords internal
NULL

#' @rdname string_get_link
string_get_link <- function(s) {
  ok <- s %in% get_names_top(NULL, link = FALSE)
  s[ok] <- replace(rep_len("log", sum(ok)), s[ok] == "p", "logit")
  s[!ok] <- NA
  s
}

#' @rdname string_get_link
string_add_link <- function(s) {
  ok <- s %in% get_names_top(NULL, link = FALSE)
  s[ok] <- sprintf("%s(%s)", string_get_link(s[ok]), s[ok])
  s[!ok] <- NA
  s
}

#' @rdname string_get_link
string_remove_link <- function(fs) {
  ok <- fs %in% get_names_top(NULL, link = TRUE)
  fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\2", fs[ok])
  fs[!ok] <- NA
  fs
}

#' @rdname string_get_link
string_extract_link <- function(fs) {
  ok <- fs %in% get_names_top(NULL, link = TRUE)
  fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\1", fs[ok])
  fs[!ok] <- NA
  fs
}

#' Get link and inverse link functions
#'
#' Retrieve the link function named by a string, or its inverse.
#'
#' @param f
#'   A \link{character} string naming a link function.
#' @param inverse
#'   A \link{logical} flag. If \code{TRUE}, then the inverse is returned.
#'
#' @return
#' A \link{function}.
#'
#' @examples
#' ## match_link("log")
#' ## match_link("log", inverse = TRUE)
#'
#' @name match_link
#' @keywords internal
#' @importFrom stats plogis qlogis
match_link <- function(f, inverse = FALSE) {
  if (inverse) {
    switch(f,
      identity = identity,
      log = exp,
      logit = function(q) plogis(q),
      stop("Link not implemented.")
    )
  } else {
    switch(f,
      identity = identity,
      log = log,
      logit = function(p) qlogis(p),
      stop("Link not implemented.")
    )
  }
}

#' Construct list of model frames
#'
#' Constructs model frames for use by \code{\link{egf}} from corresponding
#' formulae and data frames. A battery of checks is performed on the input.
#'
#' @inheritParams egf
#'
#' @return
#' A list with elements
#' \code{frame}, \code{frame_windows}, \code{frame_parameters},
#' and \code{frame_append}. See descriptions under \code{\link{egf}}.
#'
#' @keywords internal
#' @importFrom stats terms model.frame na.fail na.pass as.formula complete.cases
make_frames <- function(model,
                        formula, formula_windows, formula_parameters,
                        data, data_windows,
                        subset, subset_windows,
                        na_action, na_action_windows,
                        init, origin, append) {
  ## Utility for checking `formula` and `formula_windows`
  check_ok_formula <- function(formula) {
    s <- deparse(substitute(formula))
    a <- attributes(terms(formula))
    stop_if_not(
      a$response > 0L,
      is.call(lhs <- a$variables[[1L + a$response]]),
      lhs[[1L]] == "cbind",
      length(lhs) == 3L,
      m = sprintf("Left hand side of `%s` must be a call to `cbind` with 2 arguments.", s)
    )
    stop_if_not(
      a$intercept == 1L,
      is.null(a$offset),
      length(a$term.labels) < 2L,
      m = sprintf("Right hand side of `%s` must be 1 or have exactly one term.", s)
    )
    invisible(NULL)
  }

  ## Utility for constructing `frame` and `frame_windows`
  make_frame <- function(formula, data, subset, na.action, drop.unused.levels) {
    ## Build model frame
    group <- length(attr(terms(formula), "term.labels")) == 1L
    tt <- c(if (group) list(formula[[3L]]), as.list(formula[[2L]])[-1L])
    .names <- c(if (!group) "", vapply(tt, deparse, ""))
    tt <- Map(call, c(if (group) "as.factor", "identity", "identity"), tt, USE.NAMES = FALSE)
    cl <- match.call()
    cl[[1L]] <- quote(model.frame)
    cl$formula <- as.formula(call("~", unsplit_terms(tt)), env = environment(formula))
    cl$subset <- subset
    mf <- eval(cl, parent.frame(), baseenv())
    if (!group) {
      mf <- data.frame(rep_len(factor(1), nrow(mf)), mf)
    }
    attr(mf, "names_original") <- .names
    attr(mf, "terms") <- NULL
    mf
  }


  ### Construction step ########################################################

  ### Time series stuff

  check_ok_formula(formula)
  frame <- make_frame(
    formula = formula,
    data = data,
    subset = subset,
    na.action = na.pass
  )
  names(frame) <- c("ts", "time", "x")
  frame <- frame[!is.na(frame$ts), , drop = FALSE]
  stop_if_not(
    nrow(frame) > 0L,
    m = wrap(
      "Data frame constructed from `formula`, `data`, and `subset` ",
      "must have at least one row."
    )
  )

  ### Fitting window stuff

  check_ok_formula(formula_windows)
  frame_windows <- make_frame(
    formula = formula_windows,
    data = data_windows,
    subset = subset_windows,
    na.action = switch(na_action_windows, fail = na.fail, na.pass)

  )
  names(frame_windows) <- c("ts", "start", "end")
  stop_if_not(
    (N <- nrow(frame_windows)) > 0L,
    m = wrap(
      "Data frame constructed from ",
      "`formula_windows`, `data_windows`, and `subset_windows` ",
      "must have at least one row."
    )
  )

  ### Mixed effects stuff

  ## Names of top level nonlinear model parameters
  names_top <- get_names_top(model, link = TRUE)
  p <- length(names_top)

  if (inherits(formula_parameters, "formula")) {
    formula_parameters <- simplify_terms(formula_parameters)
    formula_parameters <- rep_len(list(formula_parameters), p)
    names(formula_parameters) <- names_top
  } else {
    names(formula_parameters) <- vapply(formula_parameters, function(x) deparse(x[[2L]]), "")
    stop_if_not(
      names(formula_parameters) %in% names_top,
      m = wrap(
        "`deparse(formula_parameters[[i]][[2L]])` must be an element of ",
        sprintf("`c(%s)`.", paste(dQuote(names_top, FALSE), collapse = ", "))
      )
    )
    formula_parameters <- lapply(formula_parameters, function(x) simplify_terms(x[-2L]))
    formula_parameters[setdiff(names_top, names(formula_parameters))] <- list(~1)
    formula_parameters <- formula_parameters[names_top]
  }

  ## Warn about fixed effects formulae without an intercept
  ## when relying on default initial values of linear coefficients
  if (is.null(init)) {
    check_ok_intercept <- function(x) {
      a <- attributes(terms(split_effects(x)$fixed))
      a$intercept == 1L || length(a$term.labels) == 0L
    }
    warn_if_not(
      vapply(formula_parameters, check_ok_intercept, FALSE),
      m = wrap(
        "Default initial values for linear fixed effects coefficients ",
        "are not reliable for fixed effects models without an intercept. ",
        "Consider setting `init` explicitly or including an intercept."
      )
    )
  }

  ## Replace `|` operators with `+` operators in formulae
  ## to enable use of `model.frame` machinery
  formula_parameters_no_bars <- lapply(formula_parameters, gsub_bar_plus)

  ## Build model frames
  cl <- call("model.frame",
    formula = NULL,
    data = quote(data_windows),
    subset = subset_windows,
    na.action = quote(switch(na_action_windows, fail = na.fail, na.pass)),
    drop.unused.levels = TRUE
  )
  frame_parameters <- rep_len(list(), p)
  names(frame_parameters) <- names_top
  for (i in seq_len(p)) {
    cl$formula <- formula_parameters_no_bars[[i]]
    frame_parameters[[i]] <- eval(cl, parent.frame(), baseenv())
  }

  ## Replace, e.g., `model.frame(~1, data)`, which is an empty data frame
  ## with 0 length and 0 rows regardless of `data`, with a data frame with
  ## 0 length and the _correct_ number of rows
  len <- lengths(frame_parameters)
  frame_parameters[len == 0L] <- list(data.frame(row.names = seq_len(N)))

  ## Test model frames for rowwise correspondence with `frame_windows`
  stop_if_not(
    vapply(frame_parameters, nrow, 0L) == N,
    m = wrap(
      "Data frames constructed from `formula_windows` and `formula_parameters` ",
      "must have a common number of rows."
    )
  )

  ### Appended stuff

  if (!is.null(append) && is.data.frame(data_windows)) {
    i <- eval(subset_windows, data_windows, environment(formula_windows))
    if (append == ".") {
      j <- setdiff(names(data_windows), unlist(lapply(frame_parameters, names), FALSE, FALSE))
    } else {
      l <- as.list(seq_along(data_windows))
      names(l) <- names(data_windows)
      j <- eval(append, l, baseenv())
    }
    frame_append <- data_windows[i, j]
  } else {
    frame_append <- data.frame(row.names = seq_len(N))
  }


  ### Validation step ##########################################################

  ### Time series stuff

  nf <- as.list(attr(frame, "names_original"))
  names(nf) <- names(frame)

  if (inherits(frame$time, c("Date", "POSIXt"))) {
    frame$time <- julian(frame$time)
  } else {
    stop_if_not(
      is.numeric(frame$time),
      m = sprintf("`%s` must be a numeric, Date, or POSIXt vector.", nf$time)
    )
  }
  stop_if_not(
    is.numeric(frame$x),
    m = sprintf("`%s` must be a numeric vector.", nf$x)
  )

  ### Fitting window stuff

  nfw <- as.list(attr(frame_windows, "names_original"))
  names(nfw) <- names(frame_windows)

  if (inherits(frame_windows$start, c("Date", "POSIXt"))) {
    frame_windows$start <- julian(frame_windows$start)
  } else {
    stop_if_not(
      is.numeric(frame_windows$start),
      m = sprintf("`%s` must be a numeric, Date, or POSIXt vector.", nfw$start)
    )
  }
  if (inherits(frame_windows$end, c("Date", "POSIXt"))) {
    frame_windows$end <- julian(frame_windows$end)
  } else {
    stop_if_not(
      is.numeric(frame_windows$start),
      m = sprintf("`%s` must be a numeric, Date, or POSIXt vector.", nfw$end)
    )
  }


  ### Mixed effects stuff

  get_names_variables <- function(x) {
    vapply(attr(terms(as.formula(call("~", x))), "variables"), deparse, "")[-1L]
  }
  check_ok_bar <- function(bar, data) {
    v <- lapply(bar[-1L], get_names_variables)
    all(vapply(data[v[[1L]]], is.numeric, FALSE),
        vapply(data[v[[2L]]], is.factor,  FALSE))
  }
  check_ok_bars <- function(formula, data) {
    bars <- split_effects(formula)$random
    all(vapply(bars, check_ok_bar, FALSE, data = data))
  }
  stop_if_not(
    mapply(check_ok_bars, formula = formula_parameters, data = frame_parameters),
    m = wrap(
      "`formula_parameters` variables on left and right hand sides ",
      "of `|` must be numeric vectors and factors, respectively."
    )
  )


  ### Subsetting step ##########################################################

  ## Discard fitting windows without observations on all
  ## `formula_windows` and `formula_parameters` variables
  cc <- do.call(complete.cases, c(list(frame_windows), frame_parameters[len > 0L]))
  if (!all(cc)) {
    frame_windows <- frame_windows[cc, , drop = FALSE]
    frame_parameters <- lapply(frame_parameters, `[`, cc, , drop = FALSE)
    frame_append <- frame_append[cc, , drop = FALSE]
    frame_windows$ts <- droplevels(frame_windows$ts)
  }

  ## Discard time series without corresponding fitting windows
  ## and fitting windows without corresponding time series
  isl <- intersect(levels(frame$ts), levels(frame_windows$ts))
  stop_if_not(
    length(isl) > 0L,
    m = "There must be at least one fitting window with corresponding time series data."
  )
  frame$ts <- factor(frame$ts, levels = isl, exclude = NULL)
  i1 <- !is.na(frame$ts)
  frame <- frame[i1, , drop = FALSE]
  frame_windows$ts <- factor(frame_windows$ts, levels = isl, exclude = NULL)
  i2 <- !is.na(frame_windows$ts)
  frame_windows <- frame_windows[i2, , drop = FALSE]
  frame_parameters <- lapply(frame_parameters, `[`, i2, , drop = FALSE)
  frame_append <- frame_append[i2, , drop = FALSE]


  ### Another validation step ##################################################

  ### Time series stuff

  stop_if_not(
    is.finite(frame$time),
    m = sprintf("`%s` must be finite (after coercion to numeric).", nf$time)
  )
  if (do_day_of_week <- (model$day_of_week > 0L)) {
    if (is.double(frame$time)) {
      stop_if_not(
        all.equal(frame$time, z <- round(frame$time)),
        m = sprintf("`%s` must be integer-valued (after coercion to numeric).", nf$time)
      )
      frame$time <- z
    }
    check_ok_diff_time <- function(x) all(diff(x) == 1)
  } else {
    check_ok_diff_time <- function(x) all(diff(x) > 0)
  }
  stop_if_not(
    tapply(frame$time, frame$ts, check_ok_diff_time),
    m = sprintf("`%s` must be increasing%s%s.",
      nf$time,
      if (do_day_of_week) " with one day spacing" else "",
      if (nzchar(nf$ts)) sprintf(" in each level of `%s`", nf$ts) else ""
    )
  )
  stop_if_not(
    frame$x[!is.na(frame$x)] >= 0,
    m = sprintf("`%s` must be non-negative.", nf$x)
  )
  if (is.double(frame$x)) {
    if (any(is_NaN_or_Inf <- is.nan(frame$x) | is.infinite(frame$x))) {
      warning(sprintf("NaN and Inf in `%s` replaced with NA.", nf$x))
      frame$x[is_NaN_or_Inf] <- NA
    }
    if (!isTRUE(all.equal(frame$x, z <- round(frame$x)))) {
      warning(sprintf("Nonintegral elements of `%s` rounded to nearest integer.", nf$x))
    }
    frame$x <- z
  }

  ### Mixed effects stuff

  any_infinite <- function(d) {
    i <- vapply(d, is.double, FALSE)
    any(vapply(d[i], function(x) any(is.infinite(x)), FALSE))
  }
  stop_if_not(
    !vapply(frame_parameters, any_infinite, FALSE),
    m = "Numeric `formula_parameters` variables must not contain Inf or -Inf."
  )


  ### Labeling step ############################################################

  ## Order everything by time series and chronologically within time series
  o1 <- order(frame$ts)
  frame <- frame[o1, , drop = FALSE]
  o2 <- do.call(order, unname(frame_windows))
  frame_windows <- frame_windows[o2, , drop = FALSE]
  frame_parameters <- lapply(frame_parameters, `[`, o2, , drop = FALSE)
  frame_append <- frame_append[o2, , drop = FALSE]

  ## Create enumerated labels for fitting windows as ordered
  N <- nrow(frame_windows)
  frame_windows$window <- gl(N, 1L, labels = sprintf("window_%0*d", 1L + as.integer(log10(N)), seq_len(N)))

  ## Create a factor grouping observations by fitting window
  make_window_segment <- function(d, dw) {
    m0 <- sprintf("Fitting windows (%s, %s]", nfw$start, nfw$end)
    if (nzchar(nf$ts)) {
      m0 <- paste(m0, "in time series", dQuote(dw$ts[1L], FALSE))
    }
    stop_if_not(
      dw$start < dw$end,
      m = paste(m0, sprintf("do not satisfy `%s < %s`.", nfw$start, nfw$end))
    )
    stop_if_not(
      dw$start[-1L] >= dw$end[-nrow(dw)],
      m = paste(m0, "are not disjoint.")
    )
    f <- function(a, b) which(d$time >= a & d$time <= b)[-1L]
    index <- Map(f, a = dw$start, b = dw$end)
    ulindex <- unlist(index, FALSE, FALSE)
    if (na_action == "fail") {
      stop_if_not(
        !anyNA(d$x[ulindex]),
        m = paste(m0, sprintf("contain missing values (instances of NA in `%s`).", nf$x))
      )
    }
    stop_if_not(
      vapply(index, function(i) sum(!is.na(d$x[i])), 0L) > 0L,
      m = paste(m0, sprintf("contain zero observations of `%s`.", nf$x))
    )
    window <- rep_len(factor(NA, levels = levels(dw$window)), nrow(d))
    window[ulindex] <- rep.int(dw$window, lengths(index))
    window
  }
  window_split <- Map(make_window_segment,
    d = split(frame, frame$ts),
    dw = split(frame_windows, frame_windows$ts)
  )
  frame$window <- unsplit(window_split, frame$ts)


  ### Cleaning up ##############################################################

  ## Ignore edge cases
  frame$x[!duplicated(frame$ts)] <- NA

  ## Contract fitting windows to narrowest intervals
  ## containing the same set of observations
  first <- which(!(is.na(frame$window) | duplicated(frame$window))) - 1L
  last  <- which(!(is.na(frame$window) | duplicated(frame$window, fromLast = TRUE)))
  frame_windows$start <- frame$time[first]
  frame_windows$end   <- frame$time[last]

  ## Order variables hierarchically
  frame <- frame[c("ts", "window", "time", "x")]
  frame_windows <- frame_windows[c("ts", "window", "start", "end")]

  ## Drop unused factor levels
  frame_parameters <- lapply(frame_parameters, droplevels)
  frame_append <- droplevels(frame_append)

  ## Discard row names
  row.names(frame) <- NULL
  row.names(frame_windows) <- NULL
  frame_parameters <- lapply(frame_parameters, `row.names<-`, NULL)
  row.names(frame_append) <- NULL

  ## Set attributes
  attr(frame, "first") <- first
  attr(frame, "last")  <- last
  attr(frame, "formula") <- formula
  attr(frame_windows, "formula") <- formula_windows
  frame_parameters <- Map(`attr<-`, frame_parameters, "terms", lapply(formula_parameters, terms))

  list(
    frame = frame,
    frame_windows = frame_windows,
    frame_parameters = frame_parameters,
    frame_append = frame_append
  )
}

#' Construct list of prior objects
#'
#' Validates and processes \code{\link{egf}} argument \code{priors}.
#'
#' @inheritParams egf
#'
#' @return
#' A named \link{list} of \code{"\link{egf_prior}"} objects obtained
#' by evaluating the right hand side of each formula listed in \code{priors}.
#'
#' @keywords internal
make_priors <- function(model, priors) {
  names_top <- get_names_top(model, link = TRUE)
  names(priors) <- vapply(priors, function(x) deparse(x[[2L]]), "")
  stop_if_not(
    names(priors) %in% names_top,
    m = wrap(
      "`deparse(priors[[i]][[2L]])` must be an element of ",
      sprintf("`c(%s)`.", paste(dQuote(names_top, FALSE), collapse = ", "))
    )
  )
  priors <- lapply(priors, function(x) eval(x[[3L]], environment(x), baseenv()))
  stop_if_not(
    vapply(priors, inherits, FALSE, "egf_prior"),
    m = "`priors[[i]][[3L]]` must evaluate to an \"egf_prior\" object."
  )
  priors[setdiff(names_top, names(priors))] <- list(NULL)
  priors[names_top]
}

#' Construct design matrices
#'
#' Utilities for constructing the design matrices \code{X} and \code{Z}
#' associated with a fixed effects model formula and a random effects term,
#' respectively.
#'
#' @param x
#'   For \code{make_X}, a \link{formula} of the form \code{~tt}.
#'   For \code{make_Z}, a \link{call} to binary operator \code{`|`}
#'   of the form \code{(tt | g)}.
#'   Here, \code{tt} is an expression composed of potentially many terms,
#'   while \code{g} is an unevaluated factor, possibly an interaction.
#' @param frame
#'   A \link[=model.frame]{model frame} listing the variables used in \code{x}.
#' @param sparse
#'   A \link{logical} flag. If \code{TRUE}, then the design matrix
#'   is returned in \link[Matrix:sparseMatrix]{sparse} format.
#'
#' @details
#' \code{make_X(x, frame, sparse)} constructs an \code{X} matrix
#' by evaluating \code{\link{model.matrix}(x, frame)}
#' or \code{\link[Matrix]{sparse.model.matrix}(x, frame)}
#' (depending on \code{sparse}) and deleting from the result
#' columns containing only zeros.
#'
#' \code{make_Z(x, frame)} constructs a \code{Z} matrix
#' following steps outlined in \code{\link{vignette}("lmer", "lme4")}.
#' It uses \code{make_X} to construct the so-called raw model matrix
#' from \code{x[[2L]]}. The result is always sparse.
#'
#' @return
#' A matrix with \link{attributes}:
#' \item{contrasts}{
#'   See \code{\link{model.matrix}}.
#'   Absent if the expression \code{tt} in \code{x = ~tt}
#'   or \code{x = (tt | g)} does not involve \link{factor}s.
#' }
#' \item{assign}{
#'   See \code{\link{model.matrix}}.
#'   Indexes \code{\link{labels}(\link{terms}(~tt))}
#'   for the expression \code{tt} in \code{x = ~tt} or \code{x = (tt | g)}.
#' }
#' \item{group}{
#'   (\code{make_Z} only.) For \code{x = (tt | g)}, a \link{factor} of
#'   length \code{\link{ncol}(Z)} with levels \code{\link{levels}(g)},
#'   useful for splitting \code{Z} into group-specific submatrices.
#' }
#'
#' @name make_X
#' @keywords internal
NULL

#' @rdname make_X
#' @importFrom stats model.matrix
#' @importFrom Matrix sparse.model.matrix
make_X <- function(x, frame, sparse) {
  f <- if (sparse) sparse.model.matrix else model.matrix
  X <- f(x, data = frame)
  j <- colSums(abs(X)) > 0
  structure(X[, j, drop = FALSE],
    contrasts = attr(X, "contrasts"),
    assign = attr(X, "assign")[j]
  )
}

#' @rdname make_X
#' @importFrom methods as
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix KhatriRao
#' @importMethodsFrom Matrix t
make_Z <- function(x, frame) {
  X <- make_X(as.formula(call("~", x[[2L]])), frame = frame, sparse = FALSE)
  ng <- vapply(split_interaction(x[[3L]]), deparse, "")
  g <- interaction(frame[ng], drop = TRUE, sep = ":", lex.order = FALSE)
  J <- t(as(g, Class = "sparseMatrix"))
  Z <- t(KhatriRao(t(J), t(X)))
  ## For consistency, we want `model.matrix()`-style names
  ## for group levels
  G <- model.matrix(as.formula(call("~", call("+", 0, x[[3L]]))), data = frame)
  j <- colSums(abs(G)) > 0
  G <- G[, j, drop = FALSE]
  colnames(Z) <- sprintf("(%s | %s)",
    rep(colnames(X), times = ncol(G)),
    rep(colnames(G), each = ncol(X))
  )
  structure(Z,
    contrasts = attr(X, "contrasts"),
    assign = rep(attr(X, "assign"), times = ncol(G)),
    group = gl(ncol(G), ncol(X), labels = levels(g))
  )
}

#' Describe design matrix columns
#'
#' A utility linking each column of a combined design matrix
#' to the corresponding nonlinear or dispersion model parameter,
#' mixed effects formula term, grouping variable, and group.
#'
#' @param x
#'   A named \link{list} of \link{formula}e of the form \code{~tt}
#'   or a named list of \link{call}s to binary operator \code{`|`}
#'   of the form \code{(tt | g)}. \code{\link{names}(x)} must indicate
#'   corresponding top level nonlinear model parameters.
#' @param m
#'   A \link{list} of \code{X} or \code{Z} matrices obtained by
#'   applying \code{\link{make_X}} or \code{\link{make_Z}} to the
#'   elements of \code{x}.
#'
#' @return
#' A \link[=data.frame]{data frame} with rows corresponding to columns
#' of the combined design matrix \code{\link{do.call}(\link{cbind}, m)},
#' and variables:
#' \item{par}{
#'   Top level nonlinear model parameter.
#' }
#' \item{term}{
#'   \link[=deparse]{Deparse}d term from right hand side of \code{`~`}
#'   or left hand side of \code{`|`}. Otherwise, \code{"(Intercept)"}.}
#' \item{group}{
#'   \link[=deparse]{Deparse}d term from right hand side of \code{`|`}.
#'   \code{\link{NA}} if not applicable.
#' }
#' \item{level}{
#'   Level of factor or interaction \code{group}.
#'   \code{\link{NA}} if not applicable.
#' }
#' \item{colname}{
#'   Column name.
#' }
#'
#' @keywords internal
#' @importFrom stats terms as.formula
make_XZ_info <- function(x, m) {
  if (length(x) == 0L) {
    s <- c("par", "term", "group", "level", "colname")
    l <- rep_len(list(factor()), length(s))
    names(l) <- s
    return(as.data.frame(l))
  }
  ncol_m <- vapply(m, ncol, 0L)
  if (inherits(x[[1L]], "formula")) { # X case
    group <- NA_character_
    level <- NA_character_
  } else { # Z case
    bar_lhs <- lapply(x, `[[`, 2L)
    bar_rhs <- lapply(x, `[[`, 3L)
    x <- lapply(bar_lhs, function(x) as.formula(call("~", x)))
    group <- rep.int(vapply(bar_rhs, deparse, ""), ncol_m)
    level <- unlist(lapply(m, attr, "group"), FALSE, FALSE)
  }
  get_term_labels <- function(formula) {
    labels(terms(formula))
  }
  rep_term_labels <- function(term_labels, assign) {
    c("(Intercept)", term_labels)[assign + 1L]
  }
  term <- Map(rep_term_labels,
    term_labels = lapply(x, get_term_labels),
    assign = lapply(m, attr, "assign")
  )
  term <- unlist(term, FALSE, FALSE)
  par <- rep.int(names(x), ncol_m)
  colname <- unlist(lapply(m, colnames), FALSE, FALSE)
  data.frame(par, term, group, level, colname, row.names = NULL, stringsAsFactors = TRUE)
}

#' Create TMB infrastructure
#'
#' Gathers necessary components of a call to \code{\link[TMB]{MakeADFun}}.
#'
#' @param frame,frame_parameters
#'   Elements of the \link{list} output of \code{\link{make_frames}}.
#' @inheritParams egf
#'
#' @return
#' A \link{list} with elements \code{data}, \code{parameters}, \code{map},
#' \code{random}, \code{profile}, \code{DLL}, and \code{silent}.
#'
#' @seealso \code{\link{make_tmb_data}}, \code{\link{make_tmb_parameters}}
#' @keywords internal
make_tmb_args <- function(model, frame, frame_parameters, priors, control,
                          do_fit, init, map) {
  tmb_data <- make_tmb_data(
    model = model,
    frame = frame,
    frame_parameters = frame_parameters,
    priors = priors,
    control = control
  )
  tmb_parameters <- make_tmb_parameters(
    tmb_data = tmb_data,
    model = model,
    frame = frame,
    frame_parameters = frame_parameters,
    do_fit = do_fit,
    init = init
  )
  if (is.null(map)) {
    tmb_map <- list()
  } else {
    len <- lengths(tmb_parameters)
    if (is.logical(map)) {
      map <- replace(gl(length(map), 1L), map, NA)
    }
    stop_if_not(
      is.factor(map),
      length(map) == sum(len),
      m = sprintf("`map` must be a factor or logical vector of length %d.", sum(len))
    )
    tmb_map <- split(map, rep.int(gl(3L, 1L, labels = names(len)), len))
    tmb_map <- lapply(tmb_map, factor, exclude = NA)
  }
  if (has_random(tmb_data)) {
    ## Declare that `b` contains random effects
    tmb_random <- "b"
  } else {
    ## Declare that there are no random effects
    tmb_random <- NULL
    ## Fix `b` and `theta` to NA_real_ since only `beta` is used
    tmb_parameters$b <- tmb_parameters$theta <- NA_real_
    tmb_map$b <- tmb_map$theta <- factor(NA)
  }
  list(
    data = tmb_data,
    parameters = tmb_parameters,
    map = tmb_map,
    random = tmb_random,
    profile = if (control$profile) "beta",
    DLL = "epigrowthfit",
    silent = (control$trace == 0L)
  )
}

#' Construct data objects for C++ template
#'
#' Gathers in a \link{list} data objects to be passed to
#' the package's C++ template via \pkg{TMB}'s \code{DATA_*} macros.
#'
#' @inheritParams make_tmb_args
#'
#' @return
#' [Below,
#' \code{N = \link{nlevels}(frame$window)}
#' is the number of fitting windows,
#' \code{n = N + \link{sum}(!\link{is.na}(frame$window))}
#' is the total number of time points associated with a fitting window, and
#' \code{p = \link{length}(frame_parameters)}
#' is the number of top level nonlinear model parameters.]
#'
#' A \link{list} inheriting from \link{class} \code{"tmb_data"},
#' with elements:
#' \item{time}{
#'   A \link[=double]{numeric} vector of length \code{n} giving
#'   times since the left endpoint of the current fitting window.
#' }
#' \item{time_seg_len}{
#'   An \link{integer} vector of length \code{N} specifying the
#'   length of each fitting window as a number of time points.
#' }
#' \item{x}{
#'   A \link[=double]{numeric} vector of length \code{n-N} giving incidence
#'   in each fitting window. \code{x[i]} in window \code{k} is the number
#'   of cases observed from \code{time[k+i-1]} to \code{time[k+i]}.
#' }
#' \item{day1}{
#'   If \code{model$day_of_week > 0}, then an \link{integer} vector
#'   of length \code{N} indicating the first day of week in each
#'   fitting window, with value \code{i} in \code{0:6} mapping to
#'   the day of week \code{i} days after the reference day specified
#'   by \code{model$day_of_week}. Otherwise, an integer vector of
#'   the same length filled with \code{-1}.
#' }
#' \item{flags}{
#'   A \link{list} with \link{integer} elements, used as flags to
#'   specify the model being estimated and to indicate what blocks
#'   of template code should be run.
#' }
#' \item{Y}{
#'   The \link[=model.offset]{offset} component of the response matrix
#'   in dense format, with \code{N} rows and \code{p} columns.
#' }
#' \item{indices}{
#'   A \link{list} with \link{integer} elements and \link{names} of the
#'   form \code{"index_link_parameter"} (e.g., \code{"index_log_r"}),
#'   giving the column 0-index of top level nonlinear model parameters
#'   (e.g., \code{log(r)}) in the response matrix. Value \code{-1}
#'   is used for parameters not belonging to the model being estimated.
#' }
#' \item{hyperparameters}{
#'   A \link{list} of \code{p} \link{numeric} vectors listing parameters
#'   of the prior distribution (if any) assigned to top level nonlinear
#'   model parameters.
#' }
#' \item{X}{
#'   The fixed effects design matrix in \link[Matrix:sparseMatrix]{sparse}
#'   or dense format, depending on \code{control$sparse_X}, with \code{N} rows.
#' }
#' \item{Xs, Xd}{
#'   If \code{control$sparse_X = TRUE}, then \code{Xs = X} and \code{Xd}
#'   is an empty dense matrix. Otherwise, \code{Xd = X} and \code{Xs} is
#'   an empty sparse matrix.
#' }
#' \item{Z}{
#'   The random effects design matrix in \link[Matrix:sparseMatrix]{sparse}
#'   format, with \code{N} rows. If there are no random effects, then \code{Z}
#'   is an empty sparse matrix.
#' }
#' \item{X_info, Z_info}{
#'   \link[=data.frame]{Data frame}s with \code{\link{ncol}(X)}
#'   and \code{\link{ncol}(Z)} rows, respectively, constructed
#'   by \code{\link{make_XZ_info}}. Row \code{j} describes the
#'   coefficient associated with column \code{j} of \code{X}
#'   or \code{Z}. \code{Z_info} has additional \link{factor}s
#'   \code{cor} and \code{vec} \link{split}ting coefficients
#'   by relation to a common covariance matrix and random vector,
#'   respectively.
#' }
#' \item{beta_index, b_index}{
#'   \link[=integer]{Integer} vectors of length \code{\link{ncol}(X)}
#'   and \code{\link{ncol}(Z)}, respectively, with values in \code{0:(p-1)}.
#'   These split the columns of \code{X} and \code{Z} by relation to
#'   a common top level nonlinear model parameter.
#' }
#' \item{beta_index_tab, b_index_tab}{
#'   \link[=integer]{Integer} vectors of length \code{p} counting the
#'   columns of \code{X} and \code{Z}, respectively, that relate to
#'   a common top level nonlinear model parameter.
#' }
#' \item{block_rows, block_cols}{
#'   \link[=integer]{Integer} vectors together giving the dimensions
#'   of each block of random effects.
#' }
#'
#' @seealso \code{\link[TMB]{MakeADFun}}
#' @keywords internal
#' @importFrom stats formula terms model.matrix model.offset
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix sparse.model.matrix KhatriRao
make_tmb_data <- function(model, frame, frame_parameters, priors, control,
                          do_fit, init) {
  ## Indices of time points associated with fitting windows
  first <- attr(frame, "first")
  last <- attr(frame, "last")
  index <- Map(seq.int, first, last)

  ## Fitting window lengths as numbers of time points
  time_seg_len <- lengths(index)
  N <- length(index)

  ## Time since earliest time point
  ulindex <- unlist(index, FALSE, FALSE)
  time <- frame$time[ulindex] - rep.int(frame$time[first], time_seg_len)

  ## Incidence
  x <- frame$x[!is.na(frame$window)]

  if (model$day_of_week > 0L) {
    ## Date during 24 hours starting at earliest time point
    Date1 <- .Date(frame$time[first])

    ## Day of week on that Date coded as an integer `i` in `0:6`.
    ## Integer `i` maps to the day of week `i` days after a reference day,
    ## which is the day `model$day_of_week` days after Saturday
    ## NB: weekdays(.Date(2L)) == "Saturday"
    day1 <- as.integer(julian(Date1, origin = .Date(2L + model$day_of_week)) %% 7)
  } else {
    day1 <- rep_len(-1L, N)
  }

  ## Top level nonlinear model parameter names
  names_top <- names(frame_parameters)
  p <- length(frame_parameters)

  ## Fixed effects formula and list of random effects terms
  ## extracted from each mixed effects formula
  formula_parameters <- lapply(frame_parameters, function(x) formula(terms(x)))
  effects <- lapply(formula_parameters, split_effects)
  fixed <- lapply(effects, `[[`, "fixed")
  random <- lapply(effects, `[[`, "random")
  random <- Map(`names<-`, random, Map(rep_len, names(random), lengths(random)))

  ## Response matrix, offset component only
  offsets <- lapply(frame_parameters, model.offset)
  offsets[vapply(offsets, is.null, FALSE)] <- list(rep_len(0, N))
  Y <- do.call(cbind, offsets)

  ## Fixed effects design matrices
  Xl <- Map(make_X, x = fixed, frame = frame_parameters, sparse = control$sparse_X)
  X <- do.call(cbind, Xl)

  ## Random effects design matrices
  Zl <- Map(make_Z, x = do.call(c, unname(random)), frame = rep.int(frame_parameters, lengths(random)))
  Z <- do.call(cbind, Zl)

  ## Nonlinear or dispersion model parameter, formula term,
  ## group (g in (terms | g)), and group level associated
  ## with each column of each matrix
  ## FIXME: Preserve contrasts?
  X_info <- make_XZ_info(x = fixed, m = Xl)
  Z_info <- make_XZ_info(x = do.call(c, unname(random)), m = Zl)
  X_info$par <- factor(X_info$par, levels = names_top)
  Z_info$par <- factor(Z_info$par, levels = names_top)

  ## Random effects coefficients factored by relation
  ## to a correlation matrix and to a random vector
  Z_info$cor <- interaction(Z_info[c("term", "group")], drop = TRUE, sep = " | ", lex.order = TRUE)
  Z_info$vec <- interaction(Z_info[c("term", "group", "level")], drop = TRUE, lex.order = TRUE)
  levels(Z_info$vec) <- seq_along(levels(Z_info$vec))

  ## Permutation of random effects coefficients producing
  ## convenient arrangement in random effects blocks,
  ## given column-major order
  o <- do.call(order, unname(Z_info[c("cor", "vec", "par")]))
  Z <- Z[, o, drop = FALSE]
  Z_info <- Z_info[o, , drop = FALSE]
  row.names(Z_info) <- NULL

  ## Coefficients factored by relation to a nonlinear or
  ## dispersion model parameter
  beta_index <- as.integer(X_info$par) - 1L
  b_index <- as.integer(Z_info$par) - 1L

  ## Number of coefficients related to each nonlinear and
  ## dispersion model parameter
  beta_index_tab <- c(table(X_info$par))
  b_index_tab <- c(table(Z_info$par))

  ## Dimensions of random effects blocks whose column vectors
  ## are related by a common covariance matrix
  block_rows <- unname(colSums(table(Z_info$par, Z_info$cor) > 0L))
  block_cols <- unname(rowSums(table(Z_info$cor, Z_info$vec) > 0L))

  ## Empty design matrices
  empty_dense_matrix <- `dim<-`(integer(0L), c(N, 0L))
  empty_sparse_matrix <- sparseMatrix(
    i = integer(0L),
    j = integer(0L),
    x = integer(0L),
    dims = c(N, 0L)
  )

  ## Parameters of priors on nonlinear and dispersion model parameters
  has_prior <- !vapply(priors, is.null, FALSE)
  hyperparameters <- rep_len(list(numeric(0L)), p)
  hyperparameters[has_prior] <- lapply(priors[has_prior], `[[`, "parameters")

  ## Flags
  flag_regularize <- rep_len(-1L, p)
  flag_regularize[has_prior] <- get_flag("prior", vapply(priors[has_prior], `[[`, "", "family"))
  flags <- list(
    flag_curve = get_flag("curve", model$curve),
    flag_excess = as.integer(model$excess),
    flag_family = get_flag("family", model$family),
    flag_day_of_week = as.integer(model$day_of_week > 0L),
    flag_regularize = flag_regularize,
    flag_trace = control$trace,
    flag_sparse_X = as.integer(control$sparse_X),
    flag_predict = 0L
  )

  ## Column indices of top level nonlinear model parameters
  ## in response matrix
  names_top_all <- get_names_top(NULL, link = TRUE)
  indices <- as.list(match(names_top_all, names_top, 0L) - 1L)
  names(indices) <- sub("^(log|logit)\\((.*)\\)$", "index_\\1_\\2", names_top_all)

  out <- list(
    time = time,
    time_seg_len = time_seg_len,
    x = x,
    day1 = day1,
    flags = flags,
    Y = Y,
    indices = indices,
    hyperparameters = hyperparameters,
    X = X,
    Xd = if (control$sparse_X) empty_dense_matrix else X,
    Xs = if (control$sparse_X) X else empty_sparse_matrix,
    X_info = X_info,
    beta_index = beta_index,
    beta_index_tab = beta_index_tab,
    Z = if (is.null(Z)) empty_sparse_matrix else Z,
    Z_info = Z_info,
    b_index = b_index,
    b_index_tab = b_index_tab,
    block_rows = block_rows,
    block_cols = block_cols
  )
  class(out) <- c("tmb_data", "list")
  out
}

#' Construct parameter objects for C++ template
#'
#' Gathers in a \link{list} parameter objects to be passed to
#' the package's C++ template via \pkg{TMB}'s \code{PARAMETER_*}
#' macros during the first likelihood evaluation.
#' See also \code{\link[TMB]{MakeADFun}}.
#'
#' @param tmb_data
#'   A \code{"\link[=make_tmb_data]{tmb_data}"} object.
#' @inheritParams make_tmb_args
#'
#' @details
#' When \code{init = NULL}, naive estimates of top level nonlinear
#' model parameters are obtained for each fitting window as follows:
#' \describe{
#' \item{\code{r}}{
#'   The slope of a linear model fit to \code{\link{log1p}(\link{cumsum}(x)))}.
#' }
#' \item{\code{alpha}}{
#'   \code{r*c0^(1-p)} if \code{curve = "subexponential"},
#'   \code{r/log(K/c0)} if \code{curve = "gompertz"}.
#'   These are the values obtained by setting the per capita growth
#'   rate at time 0 in the subexponential and Gompertz models equal
#'   to \code{r}, substituting the naive estimates of \code{r},
#'   \code{c0}, \code{K}, and \code{p}, and solving for \code{alpha}.
#' }
#' \item{\code{c0}}{
#'   \code{\link{exp}(log_c0)}, where \code{log_c0} is the intercept
#'   of a linear model fit to \code{\link{log1p}(\link{cumsum}(x))}.
#' }
#' \item{\code{tinfl}}{
#'   \code{\link{max}(time)}. This assumes that the fitting window
#'   ends near the time of a peak in incidence.
#' }
#' \item{\code{K}}{
#'   \code{2*\link{sum}(x)}. This assumes that the fitting window
#'   ends near the time of a peak in incidence _and_ that incidence
#'   is roughly symmetric about the peak.
#' }
#' \item{\code{p}}{0.95}
#' \item{\code{a, b, disp, w[123456]}}{1}
#' }
#' The naive estimates are log- or logit-transformed (all but \code{p}
#' use a log link), and the initial value of the \code{"(Intercept)"}
#' coefficient in a top level parameter's fixed effects model
#' (if that coefficient exists) is taken to be the mean of the link scale
#' estimates across fitting windows. The remaining elements of \code{beta}
#' take initial value 0. All elements of \code{b} and \code{theta}
#' take initial value 0, so that, initially, all random vectors follow
#' a standard normal distribution.
#'
#' @return
#' A \link{list} with elements:
#' \item{beta}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{length}(tmb_data$beta_index)} listing
#'   initial values for fixed effects coefficients.
#' }
#' \item{b}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{length}(tmb_data$b_index)} listing initial values
#'   for random effects coefficients (unit variance scale).
#' }
#' \item{theta}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{sum}(f(tmb_data$block_rows))}, where \code{f(n) = n*(n+1)/2},
#'   specifying initial covariance matrices for all random vectors.
#'   Let \code{x} be a segment of \code{theta} of length \code{f(n)},
#'   corresponding to a random vector \code{u} of length \code{n}.
#'   Then \code{x[1:n]} lists log standard deviations of the elements
#'   of \code{u}, and \code{x[-(1:n)]} lists (in row-major order)
#'   the subdiagonal elements of the unit lower triangular matrix
#'   \code{L} in the Cholesky factorization of the correlation matrix
#'   \code{S} of \code{u}:\cr
#'   \code{S = A \link[=matmult]{\%*\%} \link{t}(A)}\cr
#'   \code{A = D^-0.5 \link[=matmult]{\%*\%} L}\cr
#'   \code{D = \link{diag}(diag(L \link[=matmult]{\%*\%} \link{t}(L)))}
#' }
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis terms
make_tmb_parameters <- function(tmb_data, model, frame, frame_parameters,
                                do_fit, init) {
  ## Lengths of parameter objects
  f <- function(n) as.integer(n * (n + 1) / 2)
  len <- c(
    beta  = length(tmb_data$beta_index),
    b     = length(tmb_data$b_index),
    theta = sum(f(tmb_data$block_rows))
  )

  ## If user does not specify a full parameter vector
  if (is.null(init)) {
    ## Initialize each parameter object to a vector of zeros
    init_split <- lapply(len, numeric)

    ## Get names of top level nonlinear model parameters
    ## whose mixed effects formulae have an intercept
    has_intercept <- function(frame) {
      attr(terms(frame), "intercept") == 1L
    }
    nfp <- names(frame_parameters)[vapply(frame_parameters, has_intercept, FALSE)]

    ## For each of these top level nonlinear model parameters,
    ## compute the mean over all fitting windows of the naive estimate,
    ## and assign the result to the coefficient of `beta` corresponding
    ## to "(Intercept)"
    if (length(nfp) > 0L) {
      ## Time series segments
      tx_split <- split(frame[c("time", "x")], frame$window)

      ## Functions computing naive estimates given time series segments
      get_r_c0 <- function(d) {
        n <- max(2, trunc(nrow(d) / 2))
        ab <- try(coef(lm(log1p(cumsum(x)) ~ time, data = d, subset = seq_len(n), na.action = na.omit)), silent = TRUE)
        if (inherits(ab, "try-error") || any(!is.finite(ab))) {
          return(c(0.04, 1))
        }
        c(ab[[2L]], exp(ab[[1L]]))
      }
      get_tinfl <- function(d) {
        max(d$t)
      }
      get_K <- function(d) {
        2 * sum(d$x, na.rm = TRUE)
      }

      ## Naive estimates for each fitting window
      r_c0 <- vapply(tx_split, get_r_c0, c(0, 0))
      r  <- r_c0[1L, ]
      c0 <- r_c0[2L, ]
      tinfl <- vapply(tx_split, get_tinfl, 0)
      K <- vapply(tx_split, get_K, 0)
      p <- 0.95
      alpha <- switch(model$curve,
        subexponential = r * c0^(1 - p),
        gompertz = r / (log(K) - log(c0)),
        NA_real_
      )
      Y_init <- data.frame(r, alpha, c0, tinfl, K, p)
      Y_init[c("a", "b", "disp", paste0("w", 1:6))] <- 1

      ## Link transform
      Y_init[] <- Map(function(x, s) match_link(s)(x),
        x = Y_init,
        s = string_get_link(names(Y_init))
      )
      names(Y_init) <- string_add_link(names(Y_init))

      ## Index of elements of `beta` corresponding to "(Intercept)"
      i1 <- match(nfp, tmb_data$X_info$par, 0L)

      ## Assign means over windows
      init_split$beta[i1] <- colMeans(Y_init[nfp])
    }
  } else {
    ## Validate and split full parameter vector,
    ## producing `list(beta, b, theta)`
    stop_if_not(
      is.numeric(init),
      length(init) == sum(len),
      is.finite(init),
      m = sprintf("`init` must be a finite numeric vector of length %d.", sum(len))
    )
    names(init) <- NULL
    init_split <- split(init, rep.int(gl(3L, 1L, labels = names(len)), len))
  }

  ## Retain all naive estimates when debugging
  if (is.null(init) && !do_fit) {
    attr(init_split, "Y_init") <- as.matrix(Y_init[names(frame_parameters)])
  }
  init_split
}

#' Check for random effects
#'
#' Determines whether an object specifies a random effects model.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=make_tmb_data]{tmb_data}"} object.
#'
#' @return
#' \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
has_random <- function(object) {
  UseMethod("has_random", object)
}

#' @export
has_random.tmb_data <- function(object) {
  ncol(object$Z) > 0L
}

#' @export
has_random.egf <- function(object) {
  has_random(object$tmb_args$data)
}


#' Patch TMB-generated functions
#'
#' Define wrapper functions on top of \code{\link[TMB]{MakeADFun}}-generated
#' functions \code{fn} and \code{gr}, so that function and gradient evaluations
#' can retry inner optimization using fallback methods in the event that the
#' default method (usually \code{\link[TMB]{newton}}) fails.
#'
#' @param fn,gr
#'   \link[=function]{Function}s from a \code{\link[TMB]{MakeADFun}}-generated
#'   \link{list} object to be patched.
#' @param inner_optimizer
#'   A \link{list} of \code{"\link{egf_inner_optimizer}"} objects
#'   specifying inner optimization methods to be tried in turn.
#'
#' @return
#' A patched version of \link{function} \code{fn} or \code{gr}.
#'
#' @name patch_fn
#' @keywords internal
NULL

#' @rdname patch_fn
patch_fn <- function(fn, inner_optimizer) {
  e <- environment(fn)
  if (!exists(".egf_env", where = e, mode = "environment", inherits = FALSE)) {
    e$.egf_env <- new.env(parent = emptyenv())
  }
  e$.egf_env$fn <- fn
  e$.egf_env$inner_optimizer <- inner_optimizer

  last.par <- random <- inner.method <- inner.control <- .egf_env <- NULL # for `check`
  pfn <- function(x = last.par[-random], ...) {
    oim <- inner.method
    oic <- inner.control
    on.exit({
      inner.method <<- oim
      inner.control <<- oic
    })
    for (io in .egf_env$inner_optimizer) {
      inner.method <<- io$method
      inner.control <<- io$control
      v <- .egf_env$fn(x, ...)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
        return(v)
      }
    }
    warning("Unable to evaluate `fn(x)`, returning NaN.")
    NaN
  }
  environment(pfn) <- e
  pfn
}

#' @rdname patch_fn
patch_gr <- function(gr, inner_optimizer) {
  e <- environment(gr)
  if (!exists(".egf_env", where = e, mode = "environment", inherits = FALSE)) {
    e$.egf_env <- new.env(parent = emptyenv())
  }
  e$.egf_env$gr <- gr
  e$.egf_env$inner_optimizer <- inner_optimizer

  last.par <- random <- inner.method <- inner.control <- .egf_env <- NULL # for `check`
  pgr <- function(x = last.par[-random], ...) {
    oim <- inner.method
    oic <- inner.control
    on.exit({
      inner.method <<- oim
      inner.control <<- oic
    })
    n <- length(x)
    for (io in .egf_env$inner_optimizer) {
      inner.method <<- io$method
      inner.control <<- io$control
      v <- .egf_env$gr(x, ...)
      if (is.numeric(v) && length(v) == n && all(is.finite(v))) {
        return(v)
      }
    }
    warning("Unable to evaluate `gr(x)`, returning NaN.")
    NaN
  }
  environment(pgr) <- e
  pgr
}

#' Construct the combined model frame
#'
#' Joins in a single \link[=data.frame]{data frame} all mixed effects
#' \link[=model.frame]{model frames} and further variables originally
#' indicated by \code{\link{egf}} argument \code{append}.
#'
#' @param object An \code{"\link{egf}"} object.
#'
#' @details
#' If a variable name occurs in multiple mixed effects model frames,
#' then only one instance is retained. Except in unusual cases
#' (possible only if model formulae have different formula environments),
#' all instances of a variable name are identical, and no information is lost.
#'
#' Since the data frames being combined each correspond rowwise
#' to \code{object$frame_windows}, so does the result.
#'
#' @return
#' A data frame combining (in the sense of \code{\link{cbind}})
#' all of data frames in the \link{list} \code{object$frame_parameters}
#' and the data frame \code{object$frame_append}.
#'
#' @keywords internal
make_combined <- function(object) {
  stopifnot(inherits(object, "egf"))
  l <- c(unname(object$frame_parameters), list(object$frame_append))
  combined <- do.call(cbind, l)
  combined[!duplicated(names(combined))]
}
