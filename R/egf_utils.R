#' Get nonlinear model parameter names
#'
#' Returns the names used internally for nonlinear model parameters.
#'
#' @param curve,excess,distr,weekday
#'   Character or logical flags specifying an incidence model;
#'   see [egf()]. Alternatively, `curve` can be an `"egf"` object
#'   returned by [egf()] or a `"tmb_data"` object returned by
#'   [make_tmb_data()], in which case `excess`, `distr`, and
#'   `weekday` are ignored.
#' @param link
#'   A logical scalar. If `TRUE`, then a prefix indicating the
#'   link function used internally (either `"log_"` or `"logit_"`)
#'   is prepended to each parameter name.
#'
#' @return
#' A subset of `"r"`, `"alpha"`, `"c0"`, `"tinfl"`, `"K"`,
#' `"p"`, `"a"`, `"b"`, `"nbdisp"`, and `paste0("w", 1:6)`,
#' with or without prefixes depending on `link`.
#'
#' @examples
#' a <- list(
#'   curve = "logistic",
#'   excess = FALSE,
#'   distr = "nbinom",
#'   weekday = TRUE
#' )
#' pn_l0  <- do.call(get_par_names, c(a, list(link = FALSE)))
#' pn_l1  <- do.call(get_par_names, c(a, list(link = TRUE)))
#' identical(pn_l0, sub("^(log|logit)_", "", pn_l1))
#'
#' @export
get_par_names <- function(curve = NULL, excess = NULL,
                          distr = NULL, weekday = NULL,
                          link = TRUE) {
  UseMethod("get_par_names", curve)
}

#' @export
get_par_names.default <- function(curve = NULL, excess = NULL,
                                  distr = NULL, weekday = NULL,
                                  link = TRUE) {
  if (all(vapply(list(curve, distr, excess, weekday), is.null, FALSE))) {
    pn <- c("r", "alpha", "c0", "tinfl", "K", "p", "a", "b", "nbdisp",
            paste0("w", 1:6))
  } else {
    pn <- character(0L)
    if (!is.null(curve)) {
      a <- switch(curve,
        exponential    = c("r", "c0"),
        logistic       = c("r", "tinfl", "K"),
        richards       = c("r", "tinfl", "K", "a"),
        subexponential = c("alpha", "c0", "p"),
        gompertz       = c("alpha", "c0", "K")
      )
      pn <- c(pn, a)
    }
    if (!is.null(excess) && excess) {
      pn <- c(pn, "b")
    }
    if (!is.null(distr) && distr == "nbinom") {
      pn <- c(pn, "nbdisp")
    }
    if (!is.null(weekday) && weekday) {
      pn <- c(pn, paste0("w", 1:6))
    }
  }
  if (link) {
    return(add_link_string(pn))
  }
  pn
}

#' @export
get_par_names.tmb_data <- function(curve, excess, distr, weekday, link = TRUE) {
  pn <- levels(curve$Zinfo$par)
  if (link) {
    return(pn)
  }
  remove_link_string(pn)
}

#' @export
get_par_names.egf <- function(curve, excess, distr, weekday, link = TRUE) {
  get_par_names(curve$tmb_args$data, link = link)
}

#' Manipulate nonlinear model parameter names
#'
#' Prepend, strip, and retrieve link prefixes
#' from nonlinear model parameter names.
#'
#' @param s A character vector.
#'
#' @return
#' `add_link_string(s)` replaces valid `s[i]` with
#' `paste(get_link_string(s[i]), s[i], sep = "_")`.
#'
#' `remove_link_string(s)` replaces valid `s[i]` with
#' `sub("^(log|logit)_.*$", "\\1", s[i])`.
#'
#' `get_link_string(s)` replaces valid `s[i]` with
#' `"log"` or `"logit"` according to the link function
#' used internally for `s[i]`.
#'
#' In all cases, invalid `s[i]` are replaced with
#' `NA_character_`.
#'
#' @examples
#' get_link_string("r")
#' add_link_string("r")
#' remove_link_string("log_r")
#'
#' @noRd
NULL

add_link_string <- function(s) {
  ok <- s %in% get_par_names(link = FALSE)
  s[ok] <- paste(get_link_string(s[ok]), s[ok], sep = "_")
  s[!ok] <- NA_character_
  s
}

remove_link_string <- function(s) {
  ok <- s %in% get_par_names(link = TRUE)
  s[ok] <- sub("^(log|logit)_", "", s[ok])
  s[!ok] <- NA_character_
  s
}

get_link_string <- function(s) {
  ok <- s %in% get_par_names(link = FALSE)
  s[ok] <- ifelse(s[ok] == "p", "logit", "log")
  if (any(!ok)) { # avoids infinite recursion with get_par_names()
    ok_also <- !ok & s %in% get_par_names(link = TRUE)
    s[ok_also] <- sub("^(log|logit)_.*$", "\\1", s[ok_also])
    s[!(ok | ok_also)] <- NA_character_
  }
  s
}

#' Get link and inverse link
#'
#' Retrieve the link or inverse link function associated with a string.
#'
#' @param s A character string, either `"log"` or `"logit"`.
#'
#' @return
#' A function.
#'
#' @examples
#' get_link("log")
#' get_inverse_link("logit")
#' get_link("logit")
#' get_inverse_link("log")
#'
#' @noRd
NULL

#' @importFrom stats qlogis
get_link <- function(s) {
  switch(s,
    log = log,
    logit = function(p) qlogis(p),
    stop("Link not implemented.")
  )
}

#' @importFrom stats plogis
get_inverse_link <- function(s) {
  switch(s,
    log = exp,
    logit = function(q) plogis(q),
    stop("Link not implemented.")
  )
}

#' Validate model formulae and construct model frames
#'
#' Constructs model frames to be used by [egf()],
#' while performing myriad checks on the input.
#'
#' @inheritParams egf
#'
#' @details
#' `frame_ts` is not constructed by [stats::model.frame()]. Without
#' loss of generality, if `formula_ts = y ~ x | ts`, then `frame_ts`
#' has variables `x` (Date), `y` (numeric), `ts` (factor),
#' and `.window` (factor). It is obtained (roughly) by evaluating the
#' expressions for `x`, `y`, and `ts` in `data` (with `data` enclosed
#' by `environment(formula_ts)`), joining the resulting vectors and
#' `window` in a data frame, discarding rows belonging to time series
#' without fitting windows, and ordering the remaining rows according
#' to `order(ts)`.
#'
#' `frame_glmm` model frames are constructed by [stats::model.frame()].
#' They each have `nlevels(window) < nrow(frame_ts)` rows, because all
#' nonlinear model parameters (and by consequence all `formula_glmm`
#' variables) are required to be constant within fitting windows.
#'
#' @return
#' A list with elements:
#' \item{`frame_ts`}{
#'   The model frame for time series formula `formula_ts`
#'   (see Details), storing `terms(formula_ts)` as an attribute.
#' }
#' \item{`frame_glmm`}{
#'   A list containing the model frame for each mixed effects formula
#'   in `formula_glmm` (after completion). `frame_glmm[[i]]` stores
#'   `terms(formula_glmm[[i]])` as an attribute.
#' }
#'
#' @export
#' @importFrom stats terms model.frame na.pass as.formula
make_frames <- function(formula_ts, formula_glmm, data, window,
                        curve, excess, distr, weekday, na_action,
                        par_init) {
  ## Nonlinear model parameter names
  pn <- get_par_names(curve = curve, excess = excess,
                      distr = distr, weekday = weekday, link = TRUE)
  p <- length(pn)
  min_window_length <- 1L + p

  ## Check data
  stop_if_not(
    inherits(data, c("data.frame", "list", "environment")),
    m = "`data` must be a data frame, list, or environment."
  )

  ## Check time series formula
  stop_if_not(
    inherits(formula_ts, "formula"),
    m = "`formula_ts` must be a formula."
  )
  formula_ts <- expand_terms(formula_ts)
  a <- attributes(terms(formula_ts))
  stop_if_not(
    a$response > 0L,
    a$intercept == 1L,
    is.null(a$offset),
    length(a$term.labels) == 1L,
    m = paste0(
      "To be parsed correctly, `formula_ts` must have\n",
      "a response, intercept, no offsets, and exactly\n",
      "one term."
    )
  )

  ## Check time series variables
  try_eval1 <- function(expr) {
    try(eval(expr, envir = data, enclos = environment(formula_ts)))
  }
  cv <- deparse(formula_ts[[2L]])
  cases <- try_eval1(formula_ts[[2L]])
  stop_if_not(
    is.vector(cases, "numeric"),
    m = "Left hand side of `formula_ts`\nmust evaluate to a numeric vector."
  )
  stop_if_not(
    cases[!is.na(cases)] >= 0,
    m = sprintf("`%s` must be non-negative.", cv)
  )
  if (is.double(cases)) {
    ## Round to nearest integer, with warning if distance
    ## exceeds tolerance
    is_integer_within_tol <- function(x, tol = sqrt(.Machine$double.eps)) {
      abs(x - round(x)) < tol
    }
    cases[is.infinite(cases)] <- NA
    warn_if_not(
      is_integer_within_tol(cases),
      m = sprintf("Non-integer numeric elements of `%s` rounded.", cv)
    )
    cases <- round(cases)
  }
  if (is_bar(formula_ts[[3L]])) {
    tsv <- deparse(formula_ts[[3L]][[3L]])
    ts <- try_eval1(formula_ts[[3L]][[3L]])
    stop_if_not(
      is.factor(ts),
      m = "Right hand side of `|` in `formula_ts`\nmust evaluate to a factor."
    )
  } else {
    tsv <- "ts"
    ts <- rep.int(factor(1), length(cases))
    formula_ts[[3L]] <- call("|", formula_ts[[3L]], as.name(tsv))
  }
  dv <- deparse(formula_ts[[3L]][[2L]])
  date <- try_eval1(formula_ts[[3L]][[2L]])
  stop_if_not(
    inherits(date, "Date"),
    m = "Left hand side of `|` in `formula_ts`\nmust evaluate to a Date vector."
  )
  stop_if_not(
    diff(range(lengths(list(date, cases, ts)))) == 0L,
    m = sprintf("`%s`, `%s`, and `%s`\nmust have a common length.", dv, cv, tsv)
  )
  ts <- droplevels(ts, exclude = NA)
  stop_if_not(
    !anyNA(date),
    !anyNA(ts),
    m = sprintf("`%s` and `%s` must not have missing values.", dv, tsv)
  )
  stop_if_not(
    table(ts) >= min_window_length,
    m = sprintf("Time series must have length %d or greater.", min_window_length)
  )
  stop_if_not(
    tapply(date, ts, function(x) all(diff(x) > 0)),
    m = sprintf("`%s` must be increasing in each level of `%s`.", dv, tsv)
  )

  ## Construct time series model frame
  frame_ts <- data.frame(date, cases, ts)
  names(frame_ts) <- c(dv, cv, tsv)

  ## Check fitting windows
  stop_if_not(
    is.factor(window),
    length(window) == nrow(frame_ts),
    m = sprintf("`window` must be a factor of length `length(%s)`.", dv)
  )
  window <- droplevels(window, exclude = NA)
  stop_if_not(
    nlevels(window) > 0L,
    m = "`window` must have at least one used level."
  )
  stop_if_not(
    table(window) >= min_window_length,
    m = sprintf("Used levels of `window` must have\nfrequency %d or greater.", min_window_length)
  )
  stop_if_not(
    nlevels(interaction(ts, window, drop = TRUE)) == nlevels(window),
    m = sprintf("`%s` must be constant in each level of `window`.", tsv)
  )
  if (weekday) {
    stop_if_not(
      tapply(date, window, function(x) all(diff(x) == 1)),
      n = sprintf("weekday = TRUE: `%s` must have 1-day spacing\nin all fitting windows.", dv)
    )
  }
  if (na_action == "fail") {
    stop_if_not(
      !tapply(cases, window, function(x) anyNA(x[-1L])),
      m = sprintf("na_action = \"fail\": `%s` has missing values\nin at least one fitting window.", cv)
    )
  } else {
    stop_if_not(
      tapply(cases, window, function(x) sum(!is.na(x[-1L])) >= min_window_length - 1L),
      m = sprintf("`%s` has insufficient data\n(fewer than %d observations)\nin at least one fitting window.", cv, min_window_length - 1L)
    )
  }

  ## Check mixed effects formulae
  if (inherits(formula_glmm, "formula") && length(formula_glmm) == 2L) {
    formula_glmm <- rep.int(list(expand_terms(formula_glmm)), p)
    names(formula_glmm) <- pn
  } else {
    stop_if_not(
      inherits(formula_glmm, "list"),
      length(formula_glmm) > 0L,
      vapply(formula_glmm, inherits, FALSE, "formula"),
      lengths(formula_glmm) == 2L,
      !is.null(names(formula_glmm)),
      names(formula_glmm) %in% pn,
      m = paste0(
        "`formula_glmm` must be a formula of the form\n",
        "`~terms` or a named list of such formulae with\n",
        "`names(formula_glmm)` a subset of\n",
        "`get_par_names(curve, distr, excess, weekday, link = TRUE)`."
      )
    )
    formula_glmm[setdiff(pn, names(formula_glmm))] <- list(~1)
    formula_glmm <- lapply(formula_glmm[pn], expand_terms)
  }
  if (is.null(par_init)) {
    is_intercept_ok <- function(formula) {
      a <- attributes(terms(split_effects(formula)$fixed))
      a$intercept == 1L || length(a$term.labels) == 0L
    }
    warn_if_not(
      vapply(formula_glmm, is_intercept_ok, FALSE),
      m = paste0(
        "Default initial values for coefficients of\n",
        "fixed effects models without an intercept\n",
        "are not reliable. Consider including an\n",
        "intercept or setting `par_init` explicitly."
      )
    )
  }

  ## Construct mixed effects model frames
  formula_glmm_no_bars <- lapply(formula_glmm, gsub_bar_plus)
  frame_glmm <- lapply(formula_glmm_no_bars, model.frame,
    data = data,
    na.action = na.pass,
    drop.unused.levels = TRUE
  )
  fix_empty <- function(d, nrow) {
    if (!all(dim(d) == c(0L, 0L))) {
      return(d)
    }
    structure(as.data.frame(seq_len(nrow))[-1L], terms = attr(d, "terms"))
  }
  frame_glmm <- lapply(frame_glmm, fix_empty, nrow = nrow(frame_ts))
  stop_if_not(
    vapply(frame_glmm, nrow, 0L) == nrow(frame_ts),
    m = sprintf("`formula_glmm` variables must have length `length(%s)`.", dv)
  )

  ## Discard time series without fitting windows
  droplevels_subset <- function(x, subset) {
    if (is.null(dim(x))) {
      return(droplevels(x[subset]))
    }
    droplevels(x[subset, , drop = FALSE])
  }
  frame_ts[[tsv]] <- factor(frame_ts[[tsv]],
    levels = levels(droplevels_subset(frame_ts[[tsv]], !is.na(window)))
  )
  keep <- !is.na(frame_ts[[tsv]])
  frame_ts <- frame_ts[keep, , drop = FALSE]
  frame_glmm <- lapply(frame_glmm, droplevels_subset, keep)
  window <- window[keep]

  ## Order rows by time series
  ord <- order(frame_ts[[tsv]])
  frame_ts <- frame_ts[ord, , drop = FALSE]
  frame_glmm <- lapply(frame_glmm, `[`, ord, , drop = FALSE)
  window <- window[ord]

  ## Order levels of `window` by occurrence
  ## NB: Methods depend on this order, so edit with care
  window <- factor(window, levels = unique(window), exclude = NA)

  ## Check that fitting windows are contiguous
  stop_if_not(
    tapply(seq_along(window), window, function(i) all(diff(i) == 1L)),
    m = "Fitting windows must be contiguous."
  )

  ## Discard rows of mixed effects model frames
  ## not belonging to a fitting window
  frame_glmm <- lapply(frame_glmm, droplevels_subset, !is.na(window))

  ## Check mixed effects variables
  get_var_names <- function(x) {
    vapply(attr(terms(as.formula(call("~", x))), "variables")[-1L], deparse, "")
  }
  is_bar_ok <- function(bar, frame) {
    vn <- lapply(bar[-1L], get_var_names)
    all(vapply(frame[vn[[1L]]], is.vector, FALSE, "numeric"),
        vapply(frame[vn[[2L]]], is.factor, FALSE))
  }
  all_bars_ok <- function(formula, frame) {
    bars <- split_effects(formula)$random
    all(vapply(bars, is_bar_ok, FALSE, frame))
  }
  stop_if_not(
    mapply(all_bars_ok, formula = formula_glmm, frame = frame_glmm),
    m = paste0(
      "`formula_glmm` variables on left and right hand\n",
      "sides of `|` must evaluate to numeric vectors and\n",
      "factors, respectively."
    )
  )
  stop_if_not(
    !vapply(frame_glmm, anyNA, FALSE),
    m = "`formula_glmm` variables must not have missing values."
  )
  any_infinite <- function(d) {
    i <- vapply(d, is.double, FALSE)
    any(vapply(d[i], function(x) any(is.infinite(x)), FALSE))
  }
  stop_if_not(
    !vapply(frame_glmm, any_infinite, FALSE),
    m = "Numeric `formula_glmm` variables must be finite."
  )
  ## FIXME: Require exact equality for double vectors?
  stop_if_not(
    vapply(split(do.call(cbind, frame_glmm), window[!is.na(window)]), is_constant, FALSE),
    m = "`formula_glmm` variables must be constant\nin each level of `window`."
  )
  stop_if_not(
    !vapply(frame_glmm[lengths(frame_glmm) > 0L], is_constant, FALSE),
    m = "`formula_glmm` variables must not be constant\nacross all levels of `window`."
  )

  ## Keep just one row of mixed effects model frames
  ## from each fitting window, as we have required
  ## mixed effects variables to be constant in each window
  frame_glmm <- lapply(frame_glmm, droplevels_subset, !duplicated(window[!is.na(window)]))

  ## Append `window` to time series model frame
  frame_ts <- cbind(frame_ts, .window = window)

  ## Discard row names
  row.names(frame_ts) <- NULL
  frame_glmm <- lapply(frame_glmm, `row.names<-`, NULL)

  list(
    frame_ts = structure(frame_ts, terms = terms(formula_ts)),
    frame_glmm = Map(`attr<-`, frame_glmm, "terms", lapply(formula_glmm, terms))
  )
}

#' Construct design matrices
#'
#' Construct the `X` matrix associated with a fixed effects
#' model formula or the `Z` matrix associated with a random
#' effects term.
#'
#' @param x
#'   For `make_X()`, a fixed effects formula of the form `~terms`.
#'   For `make_Z()`, a random effects term of the form `(terms | group)`.
#' @param frame
#'   A model frame (a data frame with an appropriate `terms` attribute)
#'   defining the variables used in `x`.
#' @param sparse
#'   A logical scalar. If `TRUE`, then the design matrix is returned
#'   in sparse format.
#'
#' @details
#' `make_X(x, frame, sparse)` is equivalent
#' to `stats::model.matrix(x, frame)` or
#' `Matrix::sparse.model.matrix(x, frame)`
#' (depending on `sparse`), except that
#' columns containing only zeros are dropped
#' from the result.
#'
#' `make_Z()` follows the steps outlined in
#' `vignette("lmer", "lme4")` to construct a
#' `Z` matrix for `x`, using `make_X()` to
#' construct the so-called raw model matrix
#' from `x[[2L]]`. The result is always sparse.
#'
#' @return
#' A matrix with attributes:
#' \item{`contrasts`, `assign`}{
#'   See [stats::model.matrix()].
#'   `contrasts` will be absent if `foo` in `x = ~foo` or
#'   `x = (foo | group)` does not use factors. `assign` indexes
#'   `labels(terms(~foo))` for `x = ~foo` and `x = (foo | group)`.
#' }
#' \item{`group`}{
#'   `make_Z()` only. For `x = (foo | group)`, a factor of length
#'   `ncol(Z)` with levels `levels(group)`. Used to split `Z` back
#'   into group-specific submatrices.
#' }
#'
#' @noRd
NULL

#' @importFrom stats model.matrix
#' @importFrom Matrix sparse.model.matrix
make_X <- function(x, frame, sparse) {
  f <- if (sparse) sparse.model.matrix else model.matrix
  X <- f(x, data = frame)
  if (ncol(X) == 0L) {
    return(X)
  }
  j <- (colSums(abs(X)) > 0L)
  structure(X[, j, drop = FALSE],
    contrasts = attr(X, "contrasts"),
    assign = attr(X, "assign")[j]
  )
}

#' @importFrom methods as
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix KhatriRao
make_Z <- function(x, frame) {
  X <- make_X(as.formula(call("~", x[[2L]])), frame = frame, sparse = FALSE)
  gv <- vapply(split_interaction(x[[3L]]), deparse, "")
  g <- interaction(frame[gv], drop = TRUE, sep = ":")
  J <- t(as(g, Class = "sparseMatrix"))
  Z <- t(KhatriRao(t(J), t(X)))
  ## For consistency, we want model.matrix()-style names
  ## for the group levels
  G <- model.matrix(as.formula(call("~", call("+", 0, x[[3L]]))), data = frame)
  colnames(Z) <- sprintf("(%s | %s)",
    rep(colnames(X), times = nlevels(g)),
    rep(colnames(G), each = ncol(X))
  )
  structure(Z,
    contrasts = attr(X, "contrasts"),
    assign = rep(attr(X, "assign"), times = nlevels(g)),
    group = gl(nlevels(g), ncol(X), labels = levels(g))
  )
}

#' Describe design matrix columns
#'
#' A utility allowing combined design matrices and parameter vectors to
#' be split by nonlinear model parameter, term, group, or interactions
#' thereof.
#'
#' @param xl
#'   A named list of fixed effects formulae of the form `~terms`, or
#'   a named list of random effects terms of the form `(terms | group)`.
#'   `names(xl)` must specify the respective nonlinear model parameters.
#' @param ml
#'   A list of `X` or `Z` matrices obtained by applying [make_X()] or
#'   [make_Z()] to the to the elements of `xl`.
#'
#' @return
#' A data frame with `sum(vapply(ml, ncol, 0L))` giving information
#' about the columns of `do.call(cbind, ml)`. It lists factors:
#' \item{`par`}{
#'   Nonlinear model parameter.
#' }
#' \item{`term`}{
#'   Deparsed term from right hand side of ```~```
#'   or left hand side of ```|```, or `"(Intercept)"`.}
#' \item{`group`}{
#'   Deparsed term from right hand side of ```|```.
#'   `NA` if not applicable.
#' }
#' \item{`colname`}{
#'   Column name.
#' }
#'
#' @noRd
#' @importFrom stats terms as.formula
make_XZ_info <- function(xl, ml) {
  if (length(xl) == 0L) {
    cn <- c("colname", "par", "term", "group")
    m <- matrix(character(0L), ncol = 4L, dimnames = list(NULL, cn))
    return(data.frame(m, stringsAsFactors = TRUE))
  }
  ml_ncol <- vapply(ml, ncol, 0L)
  if (inherits(xl[[1L]], "formula")) { # X case
    group <- factor(NA)
  } else { # Z case
    bar_lhs <- lapply(xl, `[[`, 2L)
    bar_rhs <- lapply(xl, `[[`, 3L)
    xl <- lapply(bar_lhs, function(x) as.formula(call("~", x)))
    group <- rep(factor(vapply(bar_rhs, deparse, "")), ml_ncol)
  }
  get_term_labels <- function(formula) {
    labels(terms(formula))
  }
  rep_factor_term_labels <- function(term_labels, assign) {
    factor(c("(Intercept)", term_labels))[assign + 1L]
  }
  term <- unlist(Map(rep_factor_term_labels,
    term_labels = lapply(xl, get_term_labels),
    assign = lapply(ml, attr, "assign")
  ))
  par <- rep(factor(names(xl)), ml_ncol)
  colname <- unlist(lapply(ml, colnames))
  data.frame(par, term, group, colname, row.names = NULL, stringsAsFactors = FALSE)
}

#' Construct data objects for C++ template
#'
#' Gathers in a list data objects to be passed to the C++ template
#' via TMB's `DATA_` macros. See [TMB::MakeADFun()].
#'
#' @inheritParams make_tmb_args
#'
#' @return
#' A `"tmb_data"` object, which is a list with elements:
#' \item{`t`}{
#'   An integer vector of length `sum(!is.na(frame_ts$.window))`
#'   giving time as a number of days since the earliest time point
#'   in the current fitting window.
#' }
#' \item{`x`}{
#'   An integer vector of length
#'   `sum(!is.na(frame_ts$.window))-nlevels(frame_ts$.window)`
#'   giving incidence in each fitting window. `x[i]` in window `k` is
#'   the number of cases observed between times `t[k+i-1]` and `t[k+i]`.
#' }
#' \item{`wlen`}{
#'   An integer vector of length `nlevels(frame_ts$.window)`
#'   giving the length of each fitting window as a number of
#'   intervals (one less than the number of time points).
#' }
#' \item{`dow0`}{
#'   An integer vector giving the first weekday in each fitting
#'   window, with values `i` in `0:6` mapping to the weekday `i`
#'   days after the reference weekday specified by `weekday_ref`.
#' }
#' \item{`Yo`}{
#'   The offsets matrix in dense format.
#' }
#' \item{`X`}{
#'   The fixed effects design matrix in sparse or dense format,
#'   depending on `sparse_X`.
#' }
#' \item{`Xs`, `Xd`}{
#'   If `sparse_X = TRUE`, then `Xs = X` and `Xd` is an empty dense
#'   matrix. Otherwise, `Xd = X` and `Xs` is an empty sparse matrix.
#' }
#' \item{`Z`}{
#'   The random effects design matrix in sparse format.
#'   If there are no random effects, then `Z` is an empty matrix.
#' }
#' \item{`X_info`, `Z_info`}{
#'   Data frames with `ncol(X)` and `ncol(Z)` rows, respectively,
#'   containing factors `colname`, `par`, `term`, and `group`.
#'   Row `j` describes the coefficient associated with column `j`
#'   of `X` or `Z`. `Zinfo` has additional factors `cor` and `vec`
#'   splitting coefficients by relation to a common correlation
#'   matrix and random vector, respectively.
#' }
#' \item{`fncoef`, `rncoef`}{
#'   Integer vectors of length `p` counting the columns of `X` and `Z`,
#'   respectively, pertaining to a common nonlinear model parameter.
#' }
#' \item{`fpar`, `rpar`}{
#'   Integer vectors of length `ncol(X)` and `ncol(Z)`, respectively,
#'   with values in `0:(p-1)`. These split the columns of `X` and `Z`
#'   by relation to a common nonlinear model parameter.
#' }
#' \item{`rnrow`, `rncol`}{
#'   Integer vectors together giving the dimensions of each block of
#'   the random effects matrix.
#' }
#' \item{`curve_flag`, `excess_flag`, `distr_flag`, `weekday_flag`, `sparse_X_flag`}{
#'   Integer flags referencing `curve`, `excess`, `distr`, `weekday`,
#'   and `sparse_X`.
#' }
#' \item{`predict_flag`}{
#'   An integer flag set equal to 0 so that prediction code is not run.
#' }
#' Additional integer elements of the form `j_link_parameter`
#' (e.g., `j_log_r`) give the 0-index of nonlinear model parameter
#' names (e.g., `"log_r"`) in `levels(Xinfo$par)`. The value -1 is
#' used for parameters not belonging to the nonlinear model being fit.
#'
#' @keywords internal
#' @importFrom stats terms model.matrix model.offset
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix sparse.model.matrix KhatriRao
make_tmb_data <- function(frame_ts, frame_glmm,
                          curve, excess, distr, weekday, weekday_ref,
                          sparse_X) {
  ## Nonlinear model parameter names
  pn <- get_par_names(curve = curve, excess = excess,
                      distr = distr, weekday = weekday, link = TRUE)
  p <- length(pn)

  ## Discard dates not belonging to a fitting window
  frame_ts <- frame_ts[!is.na(frame_ts$.window), , drop = FALSE]
  window <- frame_ts$.window

  ## Time as number of days since end of earliest date
  date_to_integer <- function(x) as.integer(julian(x, origin = x[[1L]]))
  t <- do.call(c, unname(tapply(frame_ts[[1L]], window, date_to_integer, simplify = FALSE)))

  ## Incidence without unused first elements
  firsts <- which(!duplicated(window))
  x <- frame_ts[[2L]][-firsts]

  ## Fitting window length as a number of intervals
  ## (one minus number of time points)
  wlen <- as.integer(table(window)) - 1L

  ## First weekday: i in {0,...,6} maps to reference day + i
  d0 <- frame_ts[[1L]][firsts + 1L]
  dow0 <- as.integer(julian(d0, origin = as.Date("2021-02-28") + weekday_ref) %% 7) # "2021-02-28" was a Saturday

  ## Fixed effects formula and list of random effects terms
  ## extracted from each mixed effects formula
  formula_glmm <- lapply(frame_glmm, terms)
  effects <- lapply(formula_glmm, split_effects)
  fixed <- lapply(effects, `[[`, "fixed")
  random <- lapply(effects, `[[`, "random")
  random <- Map(`names<-`, random, Map(rep.int, names(random), lengths(random)))

  ## Offsets matrix
  offsets <- lapply(frame_glmm, model.offset)
  offsets[vapply(offsets, is.null, FALSE)] <- list(rep.int(0, nlevels(window)))
  Yo <- do.call(cbind, offsets)

  ## Fixed effects design matrices
  Xl <- Map(make_X, x = fixed, frame = frame_glmm, sparse = sparse_X)
  X <- do.call(cbind, Xl)

  ## Random effects design matrices
  Zl <- Map(make_Z, x = do.call(c, unname(random)), frame = rep.int(frame_glmm, lengths(random)))
  Z <- do.call(cbind, Zl)

  ## Nonlinear model parameter, formula term, and grouping term
  ## (g in (r | g)) associated with each column of each matrix
  ## FIXME: Preserve contrasts?
  X_info <- make_XZ_info(xl = fixed, ml = Xl)
  Z_info <- make_XZ_info(xl = do.call(c, unname(random)), ml = Zl)
  X_info$par <- factor(X_info$par, levels = pn)
  Z_info$par <- factor(Z_info$par, levels = pn)

  ## Random effects coefficients factored by relation
  ## to a correlation matrix and to a random vector
  Z_info$cor <- interaction(Z_info[c("term", "group")], drop = TRUE)
  Z_info$vec <- interaction(Z_info[c("term", "group", "colname")], drop = TRUE)
  levels(Z_info$cor) <- seq_len(nlevels(Z_info$cor))
  levels(Z_info$vec) <- seq_len(nlevels(Z_info$vec))

  ## Permutation of random effects coefficients producing
  ## correct arrangement in random effects blocks, given
  ## column-major order
  ord <- do.call(order, Z_info[c("cor", "vec", "par")])
  Z <- Z[, ord, drop = FALSE]
  Z_info <- Z_info[ord, , drop = FALSE]

  ## Number of coefficients for each nonlinear model parameter
  fncoef <- as.integer(table(X_info$par))
  rncoef <- as.integer(table(Z_info$par))

  ## Coefficients factored by nonlinear model parameter
  fpar <- as.integer(X_info$par) - 1L
  rpar <- as.integer(Z_info$par) - 1L

  ## Dimensions of random effects blocks whose column vectors
  ## are related by a covariance matrix
  rnrow <- as.integer(colSums(table(Z_info$par, Z_info$cor) > 0L))
  rncol <- as.integer(rowSums(table(Z_info$cor, Z_info$vec) > 0L))

  ## Empty design matrices
  empty_dense_matrix <- matrix(integer(0L), nrow(X), 0L)
  empty_sparse_matrix <- sparseMatrix(
    i = integer(0L),
    j = integer(0L),
    x = integer(0L),
    dims = c(nrow(X), 0L)
  )

  l1 <- list(
    t = t,
    x = x,
    wlen = wlen,
    dow0 = dow0,
    Yo = Yo,
    X = X,
    Xs = if (sparse_X) X else empty_sparse_matrix,
    Xd = if (sparse_X) empty_dense_matrix else X,
    Z = if (is.null(Z)) empty_sparse_matrix else Z,
    X_info = X_info,
    Z_info = Z_info,
    fpar = fpar,
    rpar = rpar,
    fncoef = fncoef,
    rncoef = rncoef,
    rnrow = rnrow,
    rncol = rncol,
    curve_flag = get_flag("curve", curve),
    excess_flag = 1L * excess,
    distr_flag = get_flag("distr", distr),
    weekday_flag = 1L * weekday,
    sparse_X_flag = 1L * sparse_X,
    predict_flag = 0L
  )
  pn0 <- get_par_names(link = TRUE)
  l2 <- as.list(match(pn0, pn, 0L) - 1L)
  names(l2) <- sprintf("j_%s", pn0)
  structure(c(l1, l2), class = c("tmb_data", "list"))
}

#' Construct parameter objects for C++ template
#'
#' Gathers in a list parameter objects to be passed to the C++ template
#' via TMB's `PARAMETER_` macros during the first likelihood evaluation.
#' See [TMB::MakeADFun()].
#'
#' @param tmb_data
#'   A `"tmb_data"` object returned by [make_tmb_data()].
#' @inheritParams make_tmb_args
#'
#' @details
#' When `par_init = NULL`, naive estimates of nonlinear model
#' parameters are obtained for each fitting window as follows:
#' \describe{
#' \item{`r`}{
#'   The slope of a linear model fit to `log1p(cumsum(x)))`.
#' }
#' \item{`alpha`}{
#'   `r*c0^(1-p)` if `curve = "subexponential"`,
#'   `r/log(K/c0)` if `curve = "gompertz"`.
#'   These are the values obtained by setting the instantaneous
#'   exponential growth rate at time 0 in the subexponential and
#'   Gompertz models equal to `r`, substituting the naive estimates
#'   of `r`, `c0`, `K`, and `p`, and solving for `alpha`.
#' }
#' \item{`c0`}{
#'   `exp(log_c0)`, where `log_c0` is the intercept of a linear model
#'   fit to `log1p(cumsum(x))`.
#' }
#' \item{`tinfl`}{
#'   `max(t)`, assuming that the fitting window ends near the time
#'   of a peak in interval incidence.
#' }
#' \item{`K`}{
#'   `2*sum(x)`, assuming that the fitting window ends near the time
#'   of a peak in interval incidence _and_ that interval incidence is
#'   roughly symmetric about the peak.
#' }
#' \item{`p`}{0.8}
#' \item{`a`, `b`, `nbdisp`, `w[123456]`}{1}
#' }
#' The naive estimates are log- or logit-transformed (all but `p`
#' use a log link), and, for each nonlinear model parameter, the
#' initial value of the `"(Intercept)"` coefficient in its fixed
#' effects model (if any) is taken to be the mean of the link scale
#' estimates across fitting windows. All other mixed effects model
#' parameters take initial value 0.
#'
#' @return
#' A list with elements:
#' \item{`beta`}{
#'   A numeric vector of length `sum(tmb_data$fncoef)`
#'   listing initial values for the fixed effects coefficients.
#' }
#' \item{`b`}{
#'   If there are random effects, then a numeric vector of length
#'   `sum(tmb_data$rncoef)` listing initial values for the (unit
#'   variance) random effects coefficients. Otherwise, `NA_real_`.
#' }
#' \item{`log_sd_b`}{
#'   If there are random effects, then a numeric vector of
#'   length `sum(tmb_data$rnrow)` listing initial values
#'   for the log standard deviations of the random effects
#'   coefficients. Otherwise, `NA_real_`.
#' }
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis terms
make_tmb_parameters <- function(tmb_data, frame_ts, frame_glmm,
                                curve, par_init) {
  ## Lengths of parameter objects
  len <- c(
    beta = sum(tmb_data$fncoef),
    b = sum(tmb_data$rncoef),
    log_sd_b = sum(tmb_data$rnrow)
  )

  ## If user does not specify a full parameter vector
  if (is.null(par_init)) {
    ## Initialize each parameter object to a vector of zeros
    ## of appropriate length
    par_init_split <- Map(rep.int, times = len, x = 0)

    ## Get names of nonlinear model parameters whose mixed
    ## effects formula has an intercept
    has_intercept <- function(frame) {
      attr(terms(frame), "intercept") == 1L
    }
    pn1 <- names(frame_glmm)[vapply(frame_glmm, has_intercept, FALSE)]

    ## For each of these nonlinear model parameters, compute
    ## and assign a default value to the coefficient of `beta`
    ## corresponding to "(Intercept)". The default value is
    ## the mean over fitting windows of the naive estimates.
    if (length(pn1) > 0L) {
      ## Time series split by fitting window
      window <- frame_ts$.window[!is.na(frame_ts$.window)]
      firsts <- which(!duplicated(window))
      tx_split <- split(data.frame(t = tmb_data$t[-firsts], x = tmb_data$x), window[-firsts])

      ## Functions for computing naive estimates given a window
      get_r_c0 <- function(d) {
        h <- max(2, trunc(nrow(d) / 2))
        ab <- tryCatch(
          coef(lm(log1p(cumsum(x)) ~ t, data = d, subset = seq_len(h), na.action = na.omit)),
          error = function(e) c(0, 0.1)
        )
        c(ab[[2L]], exp(ab[[1L]]))
      }
      get_tinfl <- function(d) {
        max(d$t)
      }
      get_K <- function(d) {
        2 * sum(d$x, na.rm = TRUE)
      }

      ## Naive estimates for all windows
      r_c0 <- vapply(tx_split, get_r_c0, c(0, 0))
      r  <- r_c0[1L, ]
      c0 <- r_c0[2L, ]
      tinfl <- vapply(tx_split, get_tinfl, 0)
      K <- vapply(tx_split, get_K, 0)
      p <- 0.8
      alpha <- switch(curve,
        subexponential = r * c0^(1 - p),
        gompertz = r / log(K / c0),
        NA_real_
      )
      Y_init <- data.frame(r, alpha, c0, tinfl, K, p)
      Y_init[c("a", "b", "nbdisp", paste0("w", 1:6))] <- 1

      ## Link transform
      Y_init[] <- Map(function(x, s) get_link(s)(x),
        x = Y_init,
        s = get_link_string(names(Y_init))
      )
      names(Y_init) <- add_link_string(names(Y_init))

      ## Index of elements of `beta` corresponding to "(Intercept)"
      i1 <- match(pn1, tmb_data$X_info$par)

      ## Assign means over windows
      par_init_split$beta[i1] <- colMeans(Y_init[pn1])
    }

  ## If user specifies a full parameter vector
  } else {
    ## Validate and split full parameter vector,
    ## producing list(beta, b, log_sd_b)
    stop_if_not(
      is.vector(par_init, "numeric"),
      length(par_init) == sum(len),
      is.finite(par_init),
      m = sprintf("`par_init` must be a finite numeric vector of length %d.", sum(len))
    )
    names(par_init) <- NULL
    par_init_split <- split(par_init, rep.int(gl(3L, 1L, labels = names(len)), len))
  }

  ## Replace unused parameter objects with NA_real_
  if (!has_random(tmb_data)) {
    par_init_split[c("b", "log_sd_b")] <- list(NA_real_)
  }
  par_init_split
}

#' Create TMB infrastructure
#'
#' Gathers necessary components of a call to [TMB::MakeADFun()].
#'
#' @param frame_ts,frame_glmm
#'   Elements of the list output of [make_frames()].
#' @inheritParams egf
#'
#' @return
#' A list with elements `data`, `parameters`, `map`, `random`,
#' `DLL`, and `silent`.
#'
#' @seealso [make_tmb_data()], [make_tmb_parameters()]
#' @keywords internal
make_tmb_args <- function(frame_ts, frame_glmm,
                          curve, distr, excess, weekday, weekday_ref,
                          sparse_X, par_init) {
  tmb_data <- make_tmb_data(
    frame_ts = frame_ts,
    frame_glmm = frame_glmm,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    weekday_ref = weekday_ref,
    sparse_X = sparse_X
  )
  tmb_parameters <- make_tmb_parameters(
    tmb_data = tmb_data,
    frame_ts = frame_ts,
    frame_glmm = frame_glmm,
    par_init = par_init
  )
  if (has_random(tmb_data)) {
    ## Nothing to fix, since all parameter objects are used
    tmb_map <- list()
    ## Declare that `b` contains random effects
    tmb_random <- "b"
  } else {
    ## Fix `b` and `log_sd_b` (in this case to NA_real_),
    ## since only `beta` is used
    tmb_map <- list(log_sd_b = factor(NA), b = factor(NA))
    ## Declare that there are no random effects
    tmb_random <- NULL
  }
  list(
    data = tmb_data,
    parameters = tmb_parameters,
    map = tmb_map,
    random = tmb_random,
    DLL = "epigrowthfit",
    silent = TRUE
  )
}

#' Check for random effects
#'
#' Determine whether an object specifies a random effects model.
#'
#' @param object
#'   An `"egf"` object returned by [egf()] or a `"tmb_data"` object
#'   returned by [make_tmb_data()].
#'
#' @return
#' `TRUE` or `FALSE`.
#'
#' @noRd
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

#' Optimize a TMB model object
#'
#' Evaluates a call to an optimizer constructed from the components
#' of a TMB model object.
#'
#' @param tmb_out The list output of [TMB::MakeADFun()].
#' @inheritParams egf
#'
#' @return
#' The list output of [stats::nlminb()], [stats::nlm()],
#' or [stats::optim()], depending on `method`.
#'
#' @noRd
#' @importFrom stats nlminb nlm optim
optim_tmb_out <- function(tmb_out, method, ...) {
  if (method == "nlminb") {
    nlminb(
      start = tmb_out$par,
      objective = tmb_out$fn,
      gradient = tmb_out$gr,
      ...
    )
  } else if (method == "nlm") {
    nlm(
      f = structure(tmb_out$fn, gradient = tmb_out$gr),
      p = tmb_out$par,
      ...
    )
  } else {
    optim(
      par = tmb_out$par,
      fn = tmb_out$fn,
      gr = tmb_out$gr,
      method = method,
      ...
    )
  }
}

#' Split up `ADREPORT()`ed variables
#'
#' Extracts reported variables from an `"sdreport"` object and
#' returns them in a list.
#'
#' @param sdreport An `"sdreport"` object returned by [TMB::sdreport()].
#'
#' @details
#' When reconstructing reported matrices, which are returned as vectors,
#' assume column-major order.
#'
#' @return
#' A named list of data frames with variables `estimate` and `se`
#' giving estimates and standard errors.
#'
#' @noRd
split_sdreport <- function(sdreport) {
  ssdr <- summary(sdreport, select = "report")
  colnames(ssdr) <- c("estimate", "se")
  lapply(split(as.data.frame(ssdr), rownames(ssdr)), "row.names<-", NULL)
}
