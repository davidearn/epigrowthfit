#' Get top level nonlinear model parameter names
#'
#' Retrieves the names used internally for top level nonlinear model
#' parameters.
#'
#' @param object
#'   An \code{"\link{egf_model}"}, \code{"\link{egf}"},
#'   or \code{"\link[=egf_make_tmb_data]{tmb_data}"} object
#'   specifying a top level nonlinear model. Otherwise,
#'   \code{\link{NULL}}.
#' @param link
#'   A \link{logical} flag. If \code{TRUE},
#'   then \code{"link(name)"} is returned instead of \code{"name"}.
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
#'
#' model <- egf_model()
#' get_names_top(model, link = FALSE)
#' get_names_top(model, link = TRUE)
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
#' for a top level nonlinear model parameter
#' and \code{f} is the name
#' (either \code{"log"} or \code{"logit"})
#' of the corresponding link function.
#'
#' @param s,fs \link[=character]{Character} vectors.
#'
#' @return
#' For strings \code{s}, \code{f}, and \code{fs} as described above,
#' \code{string_get_link(s)} returns \code{f},
#' \code{string_add_link(s)} returns \code{fs},
#' \code{string_remove_link(fs)} returns \code{s}, and
#' \code{string_extract_link(fs)} returns \code{f}.
#'
#' \code{\link{NA_character_}} is returned for invalid values
#' of the argument.
#'
#' @examples
#' ## string_get_link("r")
#' ## string_add_link("r")
#' ## string_remove_link("log(r)")
#' ## string_extract_link("log(r)")
#'
#' ## string_extract_link("invalid string", "r" , "log(r)")
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

#' @importFrom stats terms formula
egf_sanitize_formula <- function(formula) {
  s <- deparse(substitute(formula))
  tt <- terms(formula, simplify = TRUE)
  a <- attributes(tt)
  stop_if_not(
    a$response > 0L,
    is.call(lhs <- a$variables[[1L + a$response]]),
    lhs[[1L]] == "cbind",
    length(lhs) == 3L,
    m = wrap("Left hand side of ", sQuote(s), " must be a call to `cbind` with 2 arguments.")
  )
  stop_if_not(
    a$intercept == 1L,
    is.null(a$offset),
    length(a$term.labels) < 2L,
    m = wrap("Right hand side of ", sQuote(s), " must be 1 or have exactly one term.")
  )
  formula(tt)
}

#' @importFrom stats terms
egf_sanitize_formula_parameters <- function(formula_parameters, model, ignore_intercept) {
  names_top <- get_names_top(model, link = TRUE)
  if (repeated <- inherits(formula_parameters, "formula")) {
    formula_parameters <- simplify_terms(formula_parameters)
    formula_parameters <- rep_len(list(formula_parameters), length(names_top))
    names(formula_parameters) <- names_top
  } else {
    names(formula_parameters) <- vapply(formula_parameters, function(x) deparse(x[[2L]]), "")
    if (!all(names(formula_parameters) %in% names_top)) {
      stop(wrap(
        "`deparse(formula_parameters[[i]][[2L]])` must be an element of ",
        sprintf("c(%s).", paste(dQuote(names_top), collapse = ", "))
      ))
    }
    formula_parameters <- lapply(formula_parameters, function(x) simplify_terms(x[-2L]))
    formula_parameters[setdiff(names_top, names(formula_parameters))] <- list(~1)
    formula_parameters <- formula_parameters[names_top]
  }
  if (!ignore_intercept) {
    check_ok_intercept <- function(x) {
      a <- attributes(terms(split_effects(x)$fixed))
      a$intercept == 1L || length(a$term.labels) == 0L
    }
    if (repeated) {
      ok <- check_ok_intercept(formula_parameters[[1L]])
    } else {
      ok <- all(vapply(formula_parameters, check_ok_intercept, FALSE))
    }
    if (!ok) {
      warning(wrap(
        "Default initial values for linear fixed effects coefficients ",
        "are not reliable for fixed effects models without an intercept. ",
        "Consider setting `init` explicitly or including an intercept."
      ))
    }
  }
  formula_parameters
}

#' @importFrom stats terms as.formula model.frame na.fail na.pass complete.cases
egf_make_frames <- function(model,
                            formula, formula_windows, formula_parameters,
                            data, data_windows,
                            subset, subset_windows,
                            na_action, na_action_windows,
                            append) {
  ## Reused for `frame` and `frame_windows`
  make_frame <- function(formula, data, subset, na.action, drop.unused.levels) {
    group <- length(attr(terms(formula), "term.labels")) == 1L
    tt <- c(if (group) list(formula[[3L]]), as.list(formula[[2L]])[-1L])
    names_original <- c(if (!group) "", vapply(tt, deparse, ""))
    tt <- Map(call, c(if (group) "as.factor", "identity", "identity"), tt, USE.NAMES = FALSE)
    cl <- match.call()
    cl[[1L]] <- quote(model.frame)
    cl$formula <- as.formula(call("~", unsplit_terms(tt)), env = environment(formula))
    cl$subset <- subset
    mf <- eval(cl, parent.frame())
    if (!group) {
      mf <- data.frame(rep_len(factor(1), nrow(mf)), mf)
    }
    attr(mf, "names_original") <- names_original
    attr(mf, "terms") <- NULL
    mf
  }


  ### Construction step ########################################################

  ### Time series stuff

  frame <- make_frame(
    formula = formula,
    data = data,
    subset = subset,
    na.action = na.pass,
    drop.unused.levels = TRUE
  )
  names(frame) <- c("ts", "time", "x")
  frame <- frame[!is.na(frame$ts), , drop = FALSE]
  if (nrow(frame) == 0L) {
    stop(wrap(
      "Data frame constructed from `formula`, `data`, and `subset` ",
      "must have at least one row."
    ))
  }

  ### Fitting window stuff

  frame_windows <- make_frame(
    formula = formula_windows,
    data = data_windows,
    subset = subset_windows,
    na.action = switch(na_action_windows, fail = na.fail, na.pass),
    drop.unused.levels = TRUE
  )
  names(frame_windows) <- c("ts", "start", "end")
  if((N <- nrow(frame_windows)) == 0L) {
    stop(wrap(
      "Data frame constructed from ",
      "`formula_windows`, `data_windows`, and `subset_windows` ",
      "must have at least one row."
    ))
  }

  ### Mixed effects stuff

  ## Build model frames
  cl <- call("model.frame",
    formula = NULL,
    data = quote(data_windows),
    subset = subset_windows,
    na.action = switch(na_action_windows, fail = quote(na.fail), quote(na.pass)),
    drop.unused.levels = TRUE
  )
  frame_parameters <- rep_len(list(), length(formula_parameters))
  names(frame_parameters) <- names(formula_parameters)
  for (i in seq_along(frame_parameters)) {
    cl$formula <- gsub_bar_plus(formula_parameters[[i]])
    frame_parameters[[i]] <- eval(cl, parent.frame())
  }

  ## Replace, e.g., `model.frame(~1, data)`, which is an empty data frame
  ## with 0 length and 0 rows regardless of `data`, with a data frame with
  ## 0 length and the _correct_ number of rows
  len <- lengths(frame_parameters)
  frame_parameters[len == 0L] <- list(data.frame(row.names = seq_len(N)))

  ## Test model frames for rowwise correspondence with `frame_windows`
  if (any(vapply(frame_parameters, nrow, 0L) != N)) {
    stop(wrap(
      "Data frames constructed from `formula_windows` and `formula_parameters` ",
      "must have a common number of rows."
    ))
  }

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
  } else if (!is.numeric(frame$time)) {
    stop(sQuote(nf$time), " must be a numeric, Date, or POSIXt vector.")
  }
  if (!is.numeric(frame$x)) {
    stop(sQuote(nf$x), " must be a numeric vector.")
  }

  ### Fitting window stuff

  nfw <- as.list(attr(frame_windows, "names_original"))
  names(nfw) <- names(frame_windows)

  if (inherits(frame_windows$start, c("Date", "POSIXt"))) {
    frame_windows$start <- julian(frame_windows$start)
  } else if (!is.numeric(frame_windows$start)) {
    stop(sQuote(nfw$start), " must be a numeric, Date, or POSIXt vector.")
  }
  if (inherits(frame_windows$end, c("Date", "POSIXt"))) {
    frame_windows$end <- julian(frame_windows$end)
  } else if (!is.numeric(frame_windows$start)) {
    stop(sQuote(nfw$end), " must be a numeric, Date, or POSIXt vector.")
  }

  ### Mixed effects stuff

  get_names_bar_lhs <- function(formula) {
    bars <- split_effects(formula)$random
    if (length(bars) == 0L) {
      return(character(0L))
    }
    f <- function(x) {
      vapply(attr(terms(as.formula(call("~", x[[2L]]))), "variables"), deparse, "")[-1L]
    }
    unique(unlist(lapply(bars, f), FALSE, FALSE))
  }
  names_bar_lhs <- lapply(formula_parameters, get_names_bar_lhs)
  check_ok_bar_lhs <- function(data, names) {
    f <- function(x) is.double(x) || is.integer(x) || is.logical(x)
    all(vapply(data[names], f, FALSE))
  }
  if (!all(mapply(check_ok_bar_lhs, data = frame_parameters, names = names_bar_lhs))) {
    stop(wrap(
      "`formula_parameters` variables on left hand side of `|` ",
      "must be of double, integer, or logical type."
    ))
  }


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
  lts <- intersect(levels(frame$ts), levels(frame_windows$ts))
  if (length(lts) == 0L) {
    stop("There must be at least one fitting window with corresponding time series data.")
  }
  frame$ts <- factor(frame$ts, levels = lts, exclude = NULL)
  i1 <- !is.na(frame$ts)
  frame <- frame[i1, , drop = FALSE]
  frame_windows$ts <- factor(frame_windows$ts, levels = lts, exclude = NULL)
  i2 <- !is.na(frame_windows$ts)
  frame_windows <- frame_windows[i2, , drop = FALSE]
  frame_parameters <- lapply(frame_parameters, `[`, i2, , drop = FALSE)
  frame_append <- frame_append[i2, , drop = FALSE]


  ### Another validation step ##################################################

  ### Time series stuff

  if (!all(is.finite(frame$time))) {
    stop(sQuote(nf$time), " must be finite (after coercion to numeric).")
  }
  if (do_day_of_week <- (model$day_of_week > 0L)) {
    if (is.double(frame$time)) {
      if (!isTRUE(all.equal(frame$time, z <- round(frame$time)))) {
        stop(sQuote(nf$time), " must be integer-valued (after coercion to numeric).")
      }
      frame$time <- z
    }
    check_ok_diff_time <- function(x) all(diff(x) == 1)
  } else {
    check_ok_diff_time <- function(x) all(diff(x) > 0)
  }
  if (!all(tapply(frame$time, frame$ts, check_ok_diff_time))) {
    stop(wrap(
      sQuote(nf$time), " must be increasing",
      if (do_day_of_week) " with one day spacing" else "",
      if (nzchar(nf$ts)) paste0(" in each level of ", sQuote(nf$ts)) else "",
      "."
    ))
  }
  if (!all(frame$x[!is.na(frame$x)] >= 0)) {
    stop(sQuote(nf$x), " must be non-negative.")
  }
  if (is.double(frame$x)) {
    if (any(is_NaN_or_Inf <- is.nan(frame$x) | is.infinite(frame$x))) {
      warning("NaN and Inf in ", sQuote(nf$x), " replaced with NA.")
      frame$x[is_NaN_or_Inf] <- NA
    }
    if (!isTRUE(all.equal(frame$x, z <- round(frame$x)))) {
      warning("Nonintegral elements of ", sQuote(nf$x), " rounded to nearest integer.")
    }
    frame$x <- z
  }

  ### Mixed effects stuff

  check_ok_bar_lhs_double <- function(data, names) {
    data <- data[names]
    i <- vapply(data, is.double, FALSE)
    f <- function(x) any(is.infinite(x))
    !any(vapply(data[i], f, FALSE))
  }
  if (!all(mapply(check_ok_bar_lhs_double, data = frame_parameters, names = names_bar_lhs))) {
    stop(wrap(
      "Numeric `formula_parameters` variables on left hand side of `|` ",
      "must not contain Inf or -Inf."
    ))
  }


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
  lw <- sprintf("window_%0*d", 1L + as.integer(log10(N)), seq_len(N))
  frame_windows$window <- gl(N, 1L, labels = lw)

  ## Create a factor grouping observations by fitting window
  make_window_segment <- function(d, dw) {
    m0 <- sprintf("Fitting windows (%s, %s]", nfw$start, nfw$end)
    if (nzchar(nf$ts)) {
      m0 <- paste(m0, "in time series", dQuote(dw$ts[1L], FALSE))
    }
    if (!all(dw$start < dw$end)) {
      stop(wrap(m0, " do not satisfy ", nfw$start, "<", nfw$end, "."))
    }
    if (!all(dw$start[-1L] >= dw$end[-nrow(dw)])) {
      stop(wrap(m0, " are not disjoint."))
    }
    f <- function(a, b) which(d$time >= a & d$time <= b)[-1L]
    index <- Map(f, a = dw$start, b = dw$end)
    ulindex <- unlist(index, FALSE, FALSE)
    if (na_action == "fail" && anyNA(d$x[ulindex])) {
      stop(wrap(m0, " contain missing values (instances of NA in ", sQuote(nf$x), ")."))
    }
    if (any(vapply(index, function(i) sum(!is.na(d$x[i])) == 0L, FALSE))) {
      stop(wrap(m0, " contain zero observations of ", sQuote(nf$x), "."))
    }
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

  frame <- frame[c("ts", "window", "time", "x")]
  frame_windows <- frame_windows[c("ts", "window", "start", "end")]

  frame_parameters <- lapply(frame_parameters, droplevels)
  frame_append <- droplevels(frame_append)

  row.names(frame) <- NULL
  row.names(frame_windows) <- NULL
  frame_parameters <- lapply(frame_parameters, `row.names<-`, NULL)
  row.names(frame_append) <- NULL

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

egf_make_priors_top <- function(formula_priors_top, model) {
  names_top <- get_names_top(model, link = TRUE)
  res <- vector(length(names_top), mode = "list")
  names(res) <- names_top
  if (length(formula_priors_top) == 0L) {
    return(res)
  }
  names(formula_priors_top) <- vapply(formula_priors_top, function(x) deparse(x[[2L]]), "")
  if (!all(names(formula_priors_top) %in% names_top)) {
    stop(wrap(
      "`deparse(formula_priors_top[[i]][[2L]])` must be an element of ",
      sprintf("c(%s).", paste(dQuote(names_top, FALSE), collapse = ", "))
    ))
  }
  if (anyDuplicated(names(formula_priors_top)) > 0L) {
    stop(wrap(
      "At least one top level parameter has a multiply defined prior. ",
      "Priors should be specified exactly once (or not at all)."
    ))
  }
  priors <- lapply(formula_priors_top, function(x) eval(x[[3L]], environment(x)))
  if (!all(vapply(priors, inherits, FALSE, "egf_prior"))) {
    stop("`formula_priors_top[[i]][[3L]])` must evaluate to an \"egf_prior\" object.")
  }
  f <- function(x) {
    x$parameters <- lapply(x$parameters, `[[`, 1L)
    x
  }
  res[names(priors)] <- lapply(priors, f)
  res
}

egf_make_priors_bottom <- function(formula_priors_bottom, beta, theta) {
  len <- c(beta = length(beta), theta = length(theta))
  initialize <- function(s) {
    x <- vector(len[[s]], mode = "list")
    names(x) <- enum_dupl_string(rep_len(s, len[[s]]))
    x
  }
  res <- sapply(names(len)[len > 0L], initialize, simplify = FALSE)
  finalize <- function(x) {
    unlist(unname(x), FALSE, TRUE)
  }
  if (length(formula_priors_bottom) == 0L || sum(len) == 0L) {
    return(finalize(res))
  }

  lhs <- lapply(formula_priors_bottom, `[[`, 2L)
  f <- function(s) {
    function(x) {
      (is.name(x) && x == s) ||
        (is.call(x) && (x[[1L]] == "[" || x[[1L]] == "[[") &&
           is.name(x[[2L]]) && x[[2L]] == s)
    }
  }
  classify <- sapply(names(res), f, simplify = FALSE)
  lhs_parameter <- rep_len(factor(NA, levels = names(res)), length(lhs))
  for (i in seq_along(lhs)) {
    for (parameter in names(res)) {
      if (classify[[parameter]](lhs[[i]])) {
        lhs_parameter[i] <- parameter
        break
      }
    }
    if (is.na(lhs_parameter[i])) {
      s <- Reduce(function(x, y) paste(x, "or", y), sQuote(names(res)))
      stop(wrap(
        "`formula_priors_bottom[[i]][[2L]]` must be ", s, " or a call ",
        "to `[` or `[[` subsetting ", s, "."
      ))
    }
  }

  e <- lapply(res, seq_along)
  eval_lhs <- function(x, formula) eval(x, e, environment(formula))
  indices <- Map(eval_lhs, x = lhs, formula = formula_priors_bottom)
  if (sum(lengths(indices)) == 0L) {
    return(finalize(res))
  }

  rhs <- lapply(formula_priors_bottom, `[[`, 3L)
  eval_rhs <- function(x, formula) eval(x, environment(formula))
  priors <- Map(eval_rhs, x = rhs, formula = formula_priors_bottom)
  if (!all(vapply(priors, inherits, FALSE, "egf_prior"))) {
    stop("`formula_priors_bottom[[i]][[3L]]` must evaluate to an \"egf_prior\" object.")
  }
  for (l in split(indices, lhs_parameter)) {
    i <- unlist(l, FALSE, FALSE)
    if (anyNA(i)) {
      stop("Invalid index vector in `formula_priors_bottom`.")
    }
    if (anyDuplicated(i) > 0L) {
      stop(wrap(
        "At least one bottom level parameter has a multiply defined prior. ",
        "Priors should be specified exactly once (or not at all)."
      ))
    }
  }
  decompose <- function(prior, length.out) {
    if (length.out == 0L) {
      return(list())
    }
    prior$parameters <- lapply(prior$parameters, rep_len, length.out = length.out)
    f <- function(i) {
      prior$parameters <- lapply(prior$parameters, `[[`, i)
      prior
    }
    lapply(seq_len(length.out), f)
  }
  priors <- Map(decompose, prior = priors, length.out = lengths(indices))

  for (parameter in names(res)) {
    if (any(m <- lhs_parameter == parameter)) {
      i <- unlist(indices[m], FALSE, FALSE)
      value <- unlist(priors[m], FALSE, FALSE)
      res[[parameter]][i] <- value
    }
  }
  finalize(res)
}

#' Construct design matrices
#'
#' Utilities for constructing the design matrices \code{X} and \code{Z}
#' associated with a fixed effects model formula and a random effects term,
#' respectively.
#'
#' @param x
#'   For \code{egf_make_X}, a \link{formula} of the form \code{~tt}.
#'   For \code{egf_make_Z}, a \link{call} to binary operator \code{`|`}
#'   of the form \code{(tt | g)}.
#'   Here, \code{tt} is an expression composed of potentially many terms,
#'   while \code{g} is an expression, possibly an interaction, indicating
#'   a grouping.
#' @param data
#'   A \link[=model.frame]{model frame} listing the variables used in \code{x}.
#' @param sparse
#'   A \link{logical} flag. If \code{TRUE}, then the design matrix
#'   is returned in \link[Matrix:sparseMatrix]{sparse} format.
#'
#' @details
#' \code{egf_make_X(x, data, sparse)} constructs an \code{X} matrix
#' by evaluating \code{\link{model.matrix}(x, data)}
#' or \code{\link[Matrix]{sparse.model.matrix}(x, data)}
#' (depending on \code{sparse}) and deleting from the result
#' columns containing only zeros.
#'
#' \code{egf_make_Z(x, data)} constructs a \code{Z} matrix
#' following steps outlined in \code{\link{vignette}("lmer", "lme4")}.
#' It uses \code{egf_make_X} to construct the so-called raw
#' model matrix from \code{x[[2L]]}. The result is always sparse.
#'
#' @return
#' A matrix with \link{attributes}:
#' \item{assign}{
#'   See \code{\link{model.matrix}}.
#'   Indexes \code{\link{labels}(\link{terms}(~tt))}
#'   for the expression \code{tt} in \code{x = ~tt} or \code{x = (tt | g)}.
#' }
#' \item{contrasts}{
#'   See \code{\link{model.matrix}}.
#'   Absent if the expression \code{tt} in \code{x = ~tt}
#'   or \code{x = (tt | g)} does not contain terms that are
#'   \link{character} vectors or \link{factor}s.
#' }
#' \item{group}{
#'   (\code{egf_make_Z} only.) For \code{x = (tt | g)}, a \link{factor}
#'   of length \code{\link{ncol}(Z)} with levels \code{\link{levels}(g)},
#'   useful for splitting \code{Z} into group-specific submatrices.
#' }
#'
#' @name egf_make_X
#' @keywords internal
NULL

#' @rdname egf_make_X
#' @importFrom stats model.matrix
egf_make_X <- function(x, data, sparse) {
  if (sparse) {
    X <- Matrix::sparse.model.matrix(x, data = data)
    j <- Matrix::colSums(abs(X)) > 0
  } else {
    X <- model.matrix(x, data = data)
    j <- colSums(abs(X)) > 0
  }
  ## FIXME: setting attributes on an S4 object is probably a bad idea...
  structure(X[, j, drop = FALSE],
    assign = attr(X, "assign")[j],
    contrasts = attr(X, "contrasts")
  )
}

#' @rdname egf_make_X
#' @importFrom methods as
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix KhatriRao
#' @importMethodsFrom Matrix t
egf_make_Z <- function(x, data) {
  X <- egf_make_X(as.formula(call("~", x[[2L]])), data = data, sparse = FALSE)
  ng <- vapply(split_interaction(x[[3L]]), deparse, "")
  g <- interaction(data[ng], drop = TRUE, sep = ":", lex.order = FALSE)
  J <- t(as(g, Class = "sparseMatrix"))
  Z <- t(KhatriRao(t(J), t(X)))
  ## For consistency, we want `model.matrix`-style names for group levels
  G <- model.matrix(as.formula(call("~", call("+", 0, x[[3L]]))), data = data)
  j <- colSums(abs(G)) > 0
  G <- G[, j, drop = FALSE]
  colnames(Z) <- sprintf("(%s | %s)",
    rep(colnames(X), times = ncol(G)),
    rep(colnames(G), each = ncol(X))
  )
  ## FIXME: setting attributes on an S4 object is probably a bad idea...
  structure(Z,
    assign = rep(attr(X, "assign"), times = ncol(G)),
    contrasts = attr(X, "contrasts"),
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
#'   applying \code{\link{egf_make_X}} or \code{\link{egf_make_Z}}
#'   to the elements of \code{x}.
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
egf_make_XZ_info <- function(x, m) {
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

#' Construct data objects for C++ template
#'
#' Gathers in a \link{list} data objects to be passed to
#' the package's C++ template via \pkg{TMB}'s \code{DATA_*} macros.
#'
#' @inheritParams egf
#' @param frame,frame_parameters
#'   Model frames obtained from the list output of \code{egf_make_frames}.
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
#'   by \code{\link{egf_make_XZ_info}}. Row \code{j} describes the
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
#' @importFrom Matrix sparseMatrix sparse.model.matrix
egf_make_tmb_data <- function(model, frame, frame_parameters, control, fit, init) {
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
  Xl <- Map(egf_make_X, x = fixed, data = frame_parameters, sparse = control$sparse_X)
  X <- do.call(cbind, Xl)

  ## Random effects design matrices
  Zl <- Map(egf_make_Z, x = do.call(c, unname(random)), data = rep.int(frame_parameters, lengths(random)))
  Z <- do.call(cbind, Zl)

  ## Nonlinear or dispersion model parameter, formula term,
  ## group (g in (terms | g)), and group level associated
  ## with each column of each matrix
  ## FIXME: Preserve contrasts?
  X_info <- egf_make_XZ_info(x = fixed, m = Xl)
  Z_info <- egf_make_XZ_info(x = do.call(c, unname(random)), m = Zl)
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

  ## Flags
  flags <- list(
    flag_curve = get_flag("curve", model$curve),
    flag_excess = as.integer(model$excess),
    flag_family = get_flag("family", model$family),
    flag_day_of_week = as.integer(model$day_of_week > 0L),
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
#' @inheritParams egf
#' @param tmb_data
#'   A \code{"\link[=egf_make_tmb_data]{tmb_data}"} object.
#' @param
#'   Model frames obtained from the list output of \code{egf_make_frames}.
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
#' take initial value 0. All elements of \code{theta} and \code{b}
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
#' \item{b}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{length}(tmb_data$b_index)} listing initial values
#'   for random effects coefficients (unit variance scale).
#' }
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis terms
egf_make_tmb_parameters <- function(tmb_data, model, frame, frame_parameters,
                                    fit, init) {
  ## Lengths of parameter objects
  f <- function(n) as.integer(n * (n + 1) / 2)
  len <- c(
    beta  = length(tmb_data$beta_index),
    theta = sum(f(tmb_data$block_rows)),
    b     = length(tmb_data$b_index)
  )

  ## If user does not specify a full parameter vector
  if (is.null(init)) {
    ## Initialize each parameter object to a vector of zeros
    init_split <- lapply(len, numeric)

    ## Identify top level nonlinear model parameters
    ## whose mixed effects formulae have an intercept
    has_intercept <- function(frame) {
      attr(terms(frame), "intercept") == 1L
    }
    i1 <- vapply(frame_parameters, has_intercept, FALSE)

    ## For each of these top level nonlinear model parameters,
    ## compute the mean over all fitting windows of the naive estimate,
    ## and assign the result to the coefficient of `beta` corresponding
    ## to "(Intercept)"
    if (any(i1)) {
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

      ## Identify elements of `beta` corresponding to "(Intercept)"
      ## and assign means over fitting windows
      nfp <- names(frame_parameters)[i1]
      index <- match(nfp, tmb_data$X_info$par, 0L)
      init_split$beta[index] <- colMeans(Y_init[nfp])
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

  ## Retain all naive estimates if debugging
  if (is.null(init) && !fit) {
    attr(init_split, "Y_init") <- as.matrix(Y_init[names(frame_parameters)])
  }
  init_split
}

egf_make_tmb_args <- function(model, frame, frame_parameters, control, fit, init, map) {
  tmb_data <- egf_make_tmb_data(
    model = model,
    frame = frame,
    frame_parameters = frame_parameters,
    control = control
  )
  tmb_parameters <- egf_make_tmb_parameters(
    tmb_data = tmb_data,
    model = model,
    frame = frame,
    frame_parameters = frame_parameters,
    fit = fit,
    init = init
  )
  if (is.null(map)) {
    tmb_map <- list()
  } else {
    len <- lengths(tmb_parameters)
    stop_if_not(
      is.factor(map),
      length(map) == sum(len),
      m = sprintf("`map` must be a factor of length %d.", sum(len))
    )
    tmb_map <- split(map, rep.int(gl(3L, 1L, labels = names(len)), len))
    tmb_map <- lapply(tmb_map, factor)
  }
  if (has_random(tmb_data)) {
    ## Declare that `b` contains random effects
    tmb_random <- "b"
  } else {
    ## Declare that there are no random effects
    tmb_random <- NULL
    ## Fix `theta` and `b` to NA_real_ since only `beta` is used
    tmb_parameters$theta <- tmb_parameters$b <- NA_real_
    tmb_map$theta <- tmb_map$b <- factor(NA)
  }
  list(
    data = tmb_data,
    parameters = tmb_parameters,
    map = tmb_map,
    random = tmb_random,
    profile = if (control$profile) "beta" else NULL,
    DLL = "epigrowthfit",
    silent = (control$trace == 0L)
  )
}

egf_update_tmb_args <- function(tmb_args, priors_top, priors_bottom) {
  f <- function(priors) {
    n <- length(priors)
    has_prior <- !vapply(priors, is.null, FALSE)
    flag_regularize <- rep_len(-1L, n)
    flag_regularize[has_prior] <- get_flag("prior", vapply(priors[has_prior], `[[`, "", "family"))
    hyperparameters <- rep_len(list(numeric(0L)), n)
    hyperparameters[has_prior] <- lapply(priors[has_prior], function(x) unlist(x$parameters, FALSE, FALSE))
    list(flag_regularize = flag_regularize, hyperparameters = hyperparameters)
  }
  l <- f(priors_top)
  tmb_args$data$flags$flag_regularize_top <- l$flag_regularize
  tmb_args$data$hyperparameters_top <- l$hyperparameters
  l <- f(priors_bottom)
  tmb_args$data$flags$flag_regularize_bottom <- l$flag_regularize
  tmb_args$data$hyperparameters_bottom <- l$hyperparameters
  tmb_args
}

#' Check for random effects
#'
#' Determines whether an object specifies a random effects model.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf_make_tmb_data]{tmb_data}"}
#'   object.
#'
#' @return
#' \code{TRUE} or \code{FALSE}.
#'
#' @name has_random
#' @keywords internal
has_random <- function(object) {
  UseMethod("has_random", object)
}

#' @rdname has_random
#' @export
has_random.tmb_data <- function(object) {
  ncol(object$Z) > 0L
}

#' @rdname has_random
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
#' @name egf_patch_fn
#' @keywords internal
NULL

#' @rdname egf_patch_fn
egf_patch_fn <- function(fn, inner_optimizer) {
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
    NaN # no warning to avoid duplication of `optim` and `nlminb` warnings
  }
  environment(pfn) <- e
  pfn
}

#' @rdname egf_patch_fn
egf_patch_gr <- function(gr, inner_optimizer) {
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
    NaN # warning because scalar result is unexpected
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
egf_make_combined <- function(object) {
  stopifnot(inherits(object, "egf"))
  l <- c(unname(object$frame_parameters), list(object$frame_append))
  combined <- do.call(cbind, l)
  combined[!duplicated(names(combined))]
}
