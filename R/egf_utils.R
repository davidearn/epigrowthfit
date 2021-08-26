#' Utilities for data validation and processing
#'
#' Functions used inside of \code{\link{egf}} to validate and process
#' user-supplied data. Most edge cases should be caught and handled here,
#' so that functions called downstream can operate without excessive
#' conditional logic.
#'
#' @inheritParams egf
#'
#' @return
#' \describe{
#' \item{egf_sanitize_formula}{
#'   A \link{formula}.
#' }
#' \item{egf_sanitize_formula_parameters}{
#'   A \link{list} of \link{formula}e.
#' }
#' \item{egf_make_frames}{
#'   A (recursive) \link{list} of \link[=data.frame]{data frames}.
#' }
#' \item{egf_make_priors_top, egf_make_priors_bottom}{
#'   A \link{list} of \code{"\link{egf_prior"}} objects.
#' }
#' }
#'
#' @noRd
NULL

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
  names_top <- egf_get_names_top(model, link = TRUE)
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
    i <- eval_subset(subset_windows, data_windows, environment(formula_windows))
    if (append == ".") {
      j <- setdiff(names(data_windows), unlist(lapply(frame_parameters, names), FALSE, FALSE))
    } else {
      j <- eval_append(append, data_windows, baseenv())
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
  names_top <- egf_get_names_top(model, link = TRUE)
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
#' Utilities for constructing the design matrix \code{X} and \code{Z}
#' associated with a fixed effects model formula \code{~tt} or random
#' effects term \code{(tt | g)}, respectively.
#'
#' @param fixed
#'   A fixed effects \link{formula} \code{~tt}.
#' @param random
#'   A random effects term \code{(tt | g)}
#'   (\link{call} to binary operator \code{`|`}).
#' @param data
#'   A \link[=model.frame]{model frame} listing the variables
#'   used in \code{fixed} or \code{random}.
#' @param sparse
#'   A \link{logical} flag. If \code{TRUE}, then the design matrix
#'   is returned in \link[Matrix:sparseMatrix]{sparse} format.
#'
#' @details
#' \code{egf_make_X(fixed, data, sparse)}
#' constructs an \code{X} matrix by evaluating
#' \code{\link{model.matrix}(fixed, data)} or
#' \code{\link[Matrix]{sparse.model.matrix}(fixed, data)}
#' (depending on \code{sparse})
#' and deleting from the result columns containing only zeros.
#'
#' \code{egf_make_Z(x, data)} constructs a \code{Z} matrix
#' following the steps outlined in \code{\link{vignette}("lmer", "lme4")}.
#' It uses \code{egf_make_X} to construct the so-called raw model matrix
#' from \code{random[[2L]]}.
#' The result is always sparse.
#'
#' @return
#' A \link[=matrix]{dense} or \link[Matrix:sparseMatrix]{sparse} matrix
#' with \link{attributes}:
#' \item{assign}{
#'   See \code{\link{model.matrix}}.
#'   Indexes \code{\link{labels}(\link{terms}(~tt))}
#'   for the expression \code{tt} in \code{fixed ~ tt}
#'   or \code{random = (tt | g)}.
#' }
#' \item{contrasts}{
#'   See \code{\link{model.matrix}}.
#'   Absent if the expression \code{tt} in \code{fixed = ~tt}
#'   or \code{random = (tt | g)} does not contain terms
#'   evaluating to \link{character} vectors or \link{factor}s.
#' }
#' \item{index}{
#'   (\code{egf_make_Z} only.)
#'   A \link{factor} of length \code{\link{ncol}(Z)} with levels
#'   \code{\link{levels}(g)}, useful for splitting the columns of
#'   \code{Z} by group level.
#' }
#'
#' @noRd
NULL

#' @importFrom stats model.matrix
egf_make_X <- function(fixed, data, sparse) {
  if (sparse) {
    X <- Matrix::sparse.model.matrix(fixed, data = data)
    j <- Matrix::colSums(abs(X)) > 0
  } else {
    X <- model.matrix(fixed, data = data)
    j <- colSums(abs(X)) > 0
  }
  ## FIXME: setting attributes on an S4 object is potentially dangerous...
  structure(X[, j, drop = FALSE],
    assign = attr(X, "assign")[j],
    contrasts = attr(X, "contrasts")
  )
}

#' @importFrom methods as
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix KhatriRao
#' @importMethodsFrom Matrix t
egf_make_Z <- function(random, data) {
  X <- egf_make_X(as.formula(call("~", random[[2L]])), data = data, sparse = FALSE)
  ng <- vapply(split_interaction(random[[3L]]), deparse, "")
  g <- interaction(data[ng], drop = TRUE, sep = ":", lex.order = FALSE)
  J <- t(as(g, Class = "sparseMatrix"))
  Z <- t(KhatriRao(t(J), t(X)))
  ## For consistency, we want `model.matrix`-style names for group levels
  G <- model.matrix(as.formula(call("~", call("+", 0, random[[3L]]))), data = data)
  j <- colSums(abs(G)) > 0
  G <- G[, j, drop = FALSE]
  colnames(Z) <- sprintf("(%s | %s)",
    rep(colnames(X), times = ncol(G)),
    rep(colnames(G), each = ncol(X))
  )
  ## FIXME: setting attributes on an S4 object is potentially dangerous...
  structure(Z,
    assign = rep(attr(X, "assign"), times = ncol(G)),
    contrasts = attr(X, "contrasts"),
    index = gl(ncol(G), ncol(X), labels = levels(g))
  )
}

#' Combine design matrices
#'
#' Utilities for combining parameter-specific fixed effects design matrices
#' \code{X} and term-specific random effects design matrices \code{Z}.
#'
#' @param fixed
#'   A named \link{list} of fixed effects \link{formula}e \code{~tt}.
#'   \code{\link{names}(fixed)} must indicate corresponding top level
#'   nonlinear model parameters.
#' @param X
#'   A \link{list} of \code{X} matrices obtained by applying
#'   \code{egf_make_X} to the elements of \code{fixed}.
#' @param random
#'   A named list of random effects terms \code{(tt | g)}
#'   (\link{call}s to binary operator \code{`|`}).
#'   \code{\link{names}(random)} must indicate corresponding top level
#'   nonlinear model parameters.
#' @param Z
#'   A \link{list} of \code{Z} matrices obtained by applying
#'   \code{egf_make_Z} to the elements of \code{random}.
#'
#' @return
#' \code{egf_combine_X} returns the result of
#' \code{\link{do.call}(\link{cbind}, X)}.
#' \code{egf_combine_Z} returns the result of
#' \code{\link{do.call}(\link{cbind}, Z)},
#' with columns ordered
#' by relation to a common random effects term \code{(tt | g)}
#' (order of terms depends on \code{random}), then
#' by relation to a common level of grouping variable \code{g}
#' (order of levels of interactions is reverse lexicographic), then
#' by top level nonlinear model parameter
#' (order of parameters is \code{\link{unique}(\link{names}(random))}).
#'
#' In both cases, \link[=attributes]{attribute} \code{info} is
#' a \link[=data.frame]{data frame} with one row per column of
#' the returned matrix, storing details about the corresponding
#' linear coefficients:
#' \item{bottom}{
#'   Bottom level mixed effects model parameter.
#' }
#' \item{top}{
#'   Top level nonlinear model parameter whose fitted value is a function
#'   of \code{bottom}.
#' }
#' \item{term}{
#'   \link[=deparse]{Deparse}d term from \code{tt}, or \code{"(Intercept)"}.
#' }
#' \item{group}{
#'   (\code{egf_combine_Z} only.)
#'   \link[=deparse]{Deparse}d expression \code{g}.
#' }
#' \item{level}{
#'   (\code{egf_combine_Z} only.)
#'   \link[=levels]{Level} of \link{factor} or interaction indicated by
#'   \code{group}.
#' }
#' \item{colname}{
#'   Column name in design matrix.
#' }
#' \item{cor}{
#'   (\code{egf_combine_Z} only.)
#'   Correlation matrix.
#'   This is the interaction of \code{term} and \code{group}.
#' }
#' \item{vec}{
#'   (\code{egf_combine_Z} only.)
#'   Random vector.
#'   This is the interaction of \code{term}, \code{group}, and \code{level}.
#' }
#'
#' @noRd
NULL

#' @importFrom stats terms
egf_combine_X <- function(fixed, X) {
  stopifnot(length(fixed) == length(X))

  nc <- vapply(X, ncol, 0L)
  get_term_labels <- function(formula) {
    labels(terms(formula))
  }
  rep_term_labels <- function(term_labels, assign) {
    c("(Intercept)", term_labels)[1L + assign]
  }

  bottom <- enum_dupl_string(rep_len("beta", sum(nc)))
  top <- rep.int(gl(length(fixed), 1L, labels = names(fixed)), nc)
  term <- Map(rep_term_labels,
    term_labels = lapply(fixed, get_term_labels),
    assign = lapply(X, attr, "assign")
  )
  term <- factor(unlist(term, FALSE, FALSE))

  res <- do.call(cbind, X)
  attr(res, "info") <- data.frame(bottom, top, term, colname = colnames(res), row.names = NULL, stringsAsFactors = FALSE)
  res
}

#' @importFrom stats as.formula terms
egf_combine_Z <- function(random, Z) {
  stopifnot(length(random) == length(Z))
  if (length(random) == 0L) {
    return(NULL)
  }

  random_formula <- lapply(random, function(x) as.formula(call("~", x[[2L]])))
  random_group <- lapply(random, `[[`, 3L)

  nc <- vapply(Z, ncol, 0L)
  get_term_labels <- function(formula) {
    labels(terms(formula))
  }
  rep_term_labels <- function(term_labels, assign) {
    c("(Intercept)", term_labels)[1L + assign]
  }

  bottom <- enum_dupl_string(rep_len("b", sum(nc)))
  top <- factor(rep.int(names(random), nc), levels = unique(names(random)))
  term <- Map(rep_term_labels,
    term_labels = lapply(random_formula, get_term_labels),
    assign = lapply(Z, attr, "assign")
  )
  term <- factor(unlist(term, FALSE, FALSE))
  group <- factor(rep.int(vapply(random_group, deparse, ""), nc))
  level <- unlist(lapply(Z, attr, "index"), FALSE, FALSE)

  res <- do.call(cbind, Z)
  info <- data.frame(top, term, group, level, colname = colnames(res), row.names = NULL, stringsAsFactors = FALSE)

  o <- do.call(order, unname(info[c("term", "group", "level", "top")]))
  res <- res[, o, drop = FALSE]
  info <- info[o, , drop = FALSE]
  attr(res, "info") <- data.frame(bottom, info, row.names = NULL, stringsAsFactors = FALSE)
  res
}
