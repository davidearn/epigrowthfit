#' Get nonlinear model parameter names
#'
#' Returns the names used internally for nonlinear model parameters.
#'
#' @param curve,distr,excess,weekday
#'   Character or logical flags specifying an incidence model;
#'   see [egf()]. Alternatively, `curve` can be an `"egf"` object
#'   returned by [egf()] or a `"tmb_data"` object returned by
#'   [make_tmb_data()], in which case `distr`, `excess`, and
#'   `weekday` are ignored.
#' @param link
#'   A logical scalar. If `TRUE`, then a prefix indicating the
#'   link function used internally (either `"log_"` or `"logit_"`)
#'   is prepended to each parameter name.
#'
#' @return
#' A subset of `"r"`, `"alpha"`, `paste0("w", 1:6)`, `"c0"`, `"tinfl"`,
#' `"K"`, `"p"`, `"a"`, `"b"`, and `"nbdisp"`, with or without prefixes
#' depending on `link`.
#'
#' @examples
#' curve <- "exponential"
#' distr <- "pois"
#' excess <- FALSE
#' weekday <- FALSE
#' pn_l0  <- get_par_names(curve = curve, distr = distr,
#'                         excess = excess, weekday = weekday,
#'                         link = FALSE)
#' pn_l1 <- get_par_names(curve = curve, distr = distr,
#'                        excess = excess, weekday = weekday,
#'                        link = TRUE)
#' identical(pn_l0, sub("^(log|logit)_", "", pn_l1))
#'
#' @export
get_par_names <- function(curve = NULL, distr = NULL,
                          excess = NULL, weekday = NULL,
                          link = TRUE) {
  UseMethod("get_par_names", curve)
}

#' @export
get_par_names.default <- function(curve = NULL, distr = NULL,
                                  excess = NULL, weekday = NULL,
                                  link = TRUE) {
  if (all(vapply(list(curve, distr, excess, weekday), is.null, FALSE))) {
    pn <- c("r", "alpha", paste0("w", 1:6),
            "c0", "tinfl", "K", "p", "a", "b", "nbdisp")
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
    if (!is.null(distr) && distr == "nbinom") {
      pn <- c(pn, "nbdisp")
    }
    if (!is.null(excess) && excess) {
      pn <- c(pn, "b")
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
get_par_names.tmb_data <- function(curve, distr, excess, weekday, link = TRUE) {
  pn <- colnames(curve$fid)
  if (link) {
    return(pn)
  }
  remove_link_string(pn)
}

#' @export
get_par_names.egf <- function(curve, distr, excess, weekday, link = TRUE) {
  get_par_names(curve$tmb_args$data, link = link)
}

#' Manipulate nonlinear model parameter names
#'
#' Add, remove, and extract prefixes from nonlinear model
#' parameter names.
#'
#' @param s A character vector.
#'
#' @return
#' `add_link_string()` returns `s` with prefixes
#' `"log_"` and `"logit_"` prepended.
#'
#' `remove_link_string()` returns `s` with prefixes
#' `"log_"` and `"logit_"` stripped.
#'
#' `extract_link_string()` returns the prefixes of `s`
#' without underscore, either `"log"` or `"logit"`.
#'
#' @examples
#' add_link_string("r")
#' remove_link_string("log_r")
#' extract_link_string("log_r)
#'
#' @noRd
NULL

add_link_string <- function(s) {
  ok <- s %in% get_par_names(link = FALSE)
  if (!any(ok)) {
    return(s)
  }
  fmt <- ifelse(s[ok] == "p", "logit_%s", "log_%s")
  s[ok] <- sprintf(fmt, s[ok])
  s
}

remove_link_string <- function(s) {
  ok <- s %in% get_par_names(link = TRUE)
  s[ok] <- sub("^(log|logit)_", "", s[ok])
  s
}

extract_link_string <- function(s) {
  ok <- s %in% get_par_names(link = TRUE)
  s[ok] <- sub("^(log|logit)_.*$", "\\1", s[ok])
  s[!ok] <- NA_character_
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
#' get_inverse_link("logit)
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
#' @return
#' A list with elements:
#' \item{`frame_ts`}{
#'   If `formula_ts = y ~ x | ts`, then a data frame containing
#'   `x` (Date), `y` (numeric), and `ts` (factor) evaluated in
#'   `data` with `environment(formula_ts)` enclosing.
#'   `split(frame_ts, frame_ts[[3L]])` splits `frame_ts` by time series.
#' }
#' \item{`frame_glmm`}{
#'   A data frame with `nrow(frame_ts)` rows containing all
#'   of the variables present in `formula_glmm`, obtained by
#'   evaluating the elements of
#'   `as.list(attr(terms(formula_glmm[[i]]), "variables"))[-1]`
#'   in `data` with `environment(formula_glmm[[i]])` enclosing,
#'   for all `i`. Only uniquely named variables are preserved,
#'   so care should be taken if any variable not found in `data`
#'   is multiply defined in the formula environments.
#' }
#' \item{`window`}{
#'   A factor of length `nrow(frame_ts)` splitting both
#'   `frame_ts` and `frame_glmm` by fitting window.
#' }
#' `frame_ts` and `frame_glmm` have `terms(formula_ts)`
#' and `lapply(formula_glmm, terms)` as attributes.
#'
#' @export
#' @importFrom stats terms model.frame
make_frames <- function(formula_ts, formula_glmm, data, window,
                        curve, distr, excess, weekday, na_action) {
  ## Nonlinear model parameter names
  pn <- get_par_names(curve = curve, distr = distr, excess = excess, weekday = weekday, link = TRUE)
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
  stop_if_not(
    length(formula_ts) == 3L,
    length(split_terms(formula_ts)) == 1L,
    m = "`formula_ts` must have a response and exactly one term."
  )

  ## Check time series
  eval0 <- function(x) try(eval(x, envir = data, enclos = environment(formula_ts)))
  cv <- deparse(formula_ts[[2L]])
  cases <- eval0(formula_ts[[2L]])
  stop_if_not(
    is.vector(cases, "numeric"),
    m = "Left hand side of `formula_ts`\nmust evaluate to a numeric vector."
  )
  stop_if_not(
    cases[!is.na(cases)] >= 0,
    m = sprintf("`%s` must be non-negative.", cv)
  )
  if (!is.integer(cases)) {
    if (any(is.infinite(cases))) {
      warning(sprintf("Infinite numeric elements of `%s` replaced with NA.", cv))
      cases[is.infinite(cases)] <- NA
    }
    if (!all(is_integer_within_tol(cases))) {
      warning(sprintf("Non-integer numeric elements of `%s` truncated.", cv))
      cases <- trunc(cases)
    }
  }
  if (is_bar(formula_ts[[3L]])) {
    tsv <- deparse(formula_ts[[3L]][[3L]])
    ts <- eval0(formula_ts[[3L]][[3L]])
    stop_if_not(
      is.factor(ts),
      m = "Right hand side of `|` in `formula_ts`\nmust evaluate to a factor."
    )
  } else {
    tsv <- "ts"
    ts <- rep.int(factor(1), length(date))
    formula_ts[[3L]] <- call("|", formula_ts[[3L]], as.name(tsv))
  }
  dv <- deparse(formula_ts[[3L]][[2L]])
  date <- eval0(formula_ts[[3L]][[2L]])
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
  frame_ts <- data.frame(ts, date, cases)
  names(frame_ts) <- c(tsv, dv, cv)

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
  if (is_bar(formula_ts[[3L]])) {
    stop_if_not(
      nlevels(interaction(ts, window, drop = TRUE)) == nlevels(window),
      m = sprintf("`%s` must be constant in each level of `window`.", tsv)
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

  ## Construct mixed effects model frame
  formula_glmm_no_bars <- lapply(formula_glmm, gsub_bar_plus)
  frame_glmm_list <- lapply(formula_glmm_no_bars, model.frame,
    data = data,
    na.action = na.pass,
    drop.unused.levels = TRUE
  )
  frame_glmm_list <- frame_glmm_list[lengths(frame_glmm_list) > 0L]
  stop_if_not(
    vapply(frame_glmm_list, nrow, 0L) == nrow(frame_ts),
    m = sprintf("`formula_glmm` variables must have length `length(%s)`.", dv)
  )
  frame_glmm <- do.call(cbind, unname(frame_glmm_list))
  frame_glmm <- frame_glmm[unique(names(frame_glmm))]
  stop_if_not(
    !anyNA(frame_glmm),
    m = "`formula_glmm` variables must not have missing values."
  )
  stop_if_not(
    !vapply(frame_glmm[vapply(frame_glmm, is.numeric, FALSE)], function(x) any(is.infinite(x)), FALSE),
    m = "Numeric `formula_glmm` variables must be finite."
  )
  stop_if_not(
    vapply(split(frame_glmm, window), function(d) nrow(unique(d)) == 1L, FALSE),
    m = "`formula_glmm` variables must be constant\nin each level of `window`."
  )
  stop_if_not(
    vapply(frame_glmm[!is.na(window), , drop = FALSE], function(x) length(unique(x)) > 1L, FALSE),
    m = "`formula_glmm` variables must not be constant\nacross all levels of `window`."
  )
  eval_to_factor <- function(expr, enclos) {
    is.factor(try(eval(expr, envir = data, enclos = enclos)))
  }
  all_bar_rhs_eval_to_factor <- function(x) {
    tl <- split_terms(x)
    tl_is_bar <- vapply(tl, is_bar, FALSE)
    bar_rhs <- lapply(tl[tl_is_bar], `[[`, 3L)
    all(vapply(bar_rhs, eval_to_factor, FALSE, enclos = environment(x)))
  }
  stop_if_not(
    vapply(formula_glmm, all_bar_rhs_eval_to_factor, FALSE),
    m = "Formula variables on right hand side of `|`\nmust be factors."
  )

  ## Order rows by time series
  ord <- order(frame_ts[[tsv]])
  frame_ts <- frame_ts[ord, ]
  frame_glmm <- frame_glmm[ord, ]
  window <- window[ord]

  ## Check that fitting windows are contiguous
  stop_if_not(
    tapply(seq_len(N), window, function(i) all(diff(i) == 1L)),
    m = "Fitting windows must be contiguous."
  )

  ## Order `levels(window)` by occurrence
  ## NB: Certain methods depend on this order, so edit with care
  window <- factor(window, levels = unique(window), exclude = NA)

  ## Discard time series without fitting windows
  frame_ts[[tsv]] <- factor(frame_ts[[tsv]],
    levels = levels(droplevels(frame_ts[[tsv]][!is.na(window)]))
  )
  keep <- !is.na(frame_ts[[tsv]])
  frame_ts <- droplevels(frame_ts[keep, , drop = FALSE])
  frame_glmm <- droplevels(frame_glmm[keep, , drop = FALSE])
  window <- window[keep]

  ## Discard row names
  row.names(frame_ts) <- NULL
  row.names(frame_glmm) <- NULL

  list(
    frame_ts = structure(frame_ts, terms = terms(formula_ts)),
    frame_glmm = structure(frame_glmm, terms = lapply(formula_glmm, terms)),
    window = window
  )
}

#' @importFrom stats model.matrix
#' @importFrom Matrix sparse.model.matrix
make_X <- function(x, frame_glmm, sparse_X) {
  f <- if (sparse_X) sparse.model.matrix else model.matrix
  X <- f(x, data = frame_glmm)
  a <- attr(X, "assign")
  j <- which(colSums(abs(X)) > 0L)
  X <- X[, j, drop = FALSE] # drop zero columns
  attr(X, "assign") <- a[j]
  X
}

#' @importFrom methods as
#' @importFrom stats model.matrix
#' @importFrom Matrix KhatriRao
make_Z <- function(x, frame_glmm, env) {
  X <- make_X(as.formula(call("~", x[[2L]]), env = env), frame_glmm = frame_glmm, sparse_X = FALSE)
  g <- interaction(frame_glmm[vapply(split_interaction(x[[3L]]), deparse, "")], drop = TRUE, sep = ":")
  J <- t(as(g, Class = "sparseMatrix")) # vignette("lmer", "lme4")
  Z <- t(KhatriRao(t(J), t(X))) # vignette("lmer", "lme4")
  ## For consistency, we want model.matrix()-style names
  ## for the group levels
  G <- model.matrix(as.formula(call("~", call("+", 0, x[[3L]])), env = env), data = frame_glmm)
  colnames(Z) <- sprintf("(%s | %s)",
    rep(colnames(X), times = nlevels(g)),
    rep(colnames(G), each = ncol(X))
  )
  attr(Z, "assign") <- rep.int(attr(X, "assign"), nlevels(g))
  Z
}

#' @importFrom stats terms
make_XZinfo <- function(xl, ml, s) {
  ml_ncol <- vapply(ml, ncol, 0L)
  if (inherits(xl[[1L]], "formula")) { # X case
    group <- factor(NA)
  } else { # Z case
    bar_lhs <- lapply(xl, `[[`, 2L)
    bar_rhs <- lapply(xl, `[[`, 3L)
    xl <- lapply(bar_lhs, function(x) as.formula(call("~", x)))
    names(xl) <- rep.int(s, length(xl))
    group <- factor(rep.int(vapply(bar_rhs, deparse, ""), ml_ncol))
  }
  get_term_labels <- function(formula) {
    labels(terms(formula))
  }
  rep_term_labels <- function(term_labels, assign) {
    c("(Intercept)", term_labels)[assign + 1L]
  }
  colname <- unlist(lapply(ml, colnames))
  par <- factor(rep.int(names(xl), ml_ncol))
  term <- factor(unlist(Map(rep_term_labels,
    term_labels = lapply(xl, get_term_labels),
    assign = lapply(ml, attr, "assign")
  )))
  data.frame(colname, par, term, group, row.names = NULL, stringsAsFactors = FALSE)
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
#' \item{`window`}{
#'   A copy of the so-named argument with missing values removed.
#' }
#' \item{`t`}{
#'   An integer vector of length `nrow(frame_ts)` giving time as
#'   a number of days since the earliest date in the relevant
#'   fitting window.
#' }
#' \item{`x`}{
#'   An integer vector of length `nrow(frame_ts)`. Within each fitting
#'   window, the first element is `NA` and the remaining elements
#'   are such that `x[i]` is the number of cases observed from the
#'   end of day `t[i-1]` to the end of day `t[i]`.
#' }
#' \item{`wl`}{
#'   An integer vector of length `nlevels(window)` giving the length
#'   of each fitting window. Equal to `as.integer(table(window))`.
#' }
#' \item{`dow0`}{
#'   An integer vector giving the day-of-week of the earliest date
#'   in each fitting window, with Sunday mapping to 0, Monday mapping
#'   to 1, and so on.
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
#'   The random effects design matrix in sparse format. If there are
#'   no random effects, then `Z` an empty matrix.
#' }
#' \item{`Xinfo`,`Zinfo`}{
#'   Data frames with `ncol(X)` and `ncol(Z)` rows, respectively,
#'   and variables `par`, `term`, `group` (`Zinfo` only) and `coef`
#'   (equal to `colnames(X)` or `colnames(Z)`, respectively). Each
#'   row describes the coefficient associated with the corresponding
#'   column of `X` or `Z`. `X[, Xinfo$par == u, drop = FALSE]` and
#'   `Z[, Zinfo$par == u, drop = FALSE]` extract the fixed and random
#'   effects design matrices for nonlinear model parameter `u`.
#' }
#' \item{`fnc`, `rnc`}{
#'   Integer vectors listing the number of columns of `X` and `Z`,
#'   respectively, pertaining to each nonlinear model parameter.
#'   Equal to
#'   `as.integer(table(Xinfo$par))` and `as.integer(table(Zinfo$par))`.
#' }
#' \item{`curve_flag`, `distr_flag`, `excess_flag`, `weekday_flag`, `sparse_X_flag`}{
#'   Integer flags referencing `curve`, `distr`, `excess`, `weekday`,
#'   and `sparse_X`.
#' }
#' \item{`predict_flag`}{
#'   An integer flag set equal to 0 so that prediction code is not run.
#' }
#' Additional integer elements of the form `j_link_parameter`
#' (e.g., `j_log_r`) give the index of nonlinear model parameters
#' (e.g., `"log_r"`) in `levels(Xinfo$par)`. Indexing starts at 0,
#' and -1 indicates that the parameter does not belong to the
#' nonlinear model being fit.
#'
#' @keywords internal
#' @importFrom stats terms model.matrix
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix sparse.model.matrix KhatriRao
make_tmb_data <- function(frame_ts, frame_glmm, window,
                          curve, distr, excess, weekday, sparse_X) {
  ## Nonlinear model parameter names (with prefixes)
  pn <- get_par_names(curve = curve, distr = distr, excess = excess, weekday = weekday, link = TRUE)
  p <- length(pn)

  ## Subset rows belonging to a fitting window
  keep <- !is.na(window)
  frame_ts <- droplevels(frame_ts[keep, , drop = FALSE])
  frame_glmm <- droplevels(frame_glmm[keep, , drop = FALSE])
  window <- window[keep]

  ## Compute time as number of days since start of fitting window
  date_to_integer <- function(x) as.integer(julian(x, origin = x[1L]))
  t <- do.call(c, unname(tapply(frame_ts[[2L]], window, date_to_integer, simplify = FALSE)))

  ## Set unused first observations to NA
  x <- do.call(c, unname(tapply(frame_ts[[3L]], window, `[<-`, 1L, NA, simplify = FALSE)))

  ## Get fitting window lengths
  wl <- as.integer(table(window))

  ## Get day-of-week of first date, mapping Sunday to 0, etc.
  d0 <- do.call(c, unname(tapply(frame_ts[[2L]], window, `[`, 1L, simplify = FALSE)))
  dow0 <- as.integer(julian(d0, origin = as.Date("1995-03-05")) %% 7) # "1995-03-05" was a Sunday

  ## Extract fixed effects formula and list of random effects terms
  ## from each mixed effects formula
  effects <- lapply(terms(frame_glmm), split_effects)
  fixed <- lapply(effects, `[[`, "fixed")
  random <- lapply(effects, `[[`, "random")

  ## Construct fixed effects design matrix
  Xl <- lapply(fixed, make_X, frame_glmm = frame_glmm, sparse_X = sparse_X)
  X <- do.call(cbind, Xl)

  ## Contruct random effects design matrix
  make_Zl <- function(xl, frame_glmm, env) {
    lapply(xl, make_Z, frame_glmm = frame_glmm, env = env)
  }
  Zll <- Map(make_Zl, xl = random, frame_glmm = list(frame_glmm), env = lapply(fixed, environment))
  Z <- do.call(cbind, do.call(c, Zll))

  ## Record nonlinear model parameter, formula term, and grouping term
  ## (g in (r | g)) associated with each coefficient of each matrix
  ## FIXME: Preserve contrasts?
  Xinfo <- make_XZinfo(xl = fixed, ml = Xl)
  ok <- (lengths(random) > 0L)
  Zinfo <- do.call(rbind, Map(make_XZinfo,
    xl = random[ok],
    ml = Zll[ok],
    s = names(random)[ok],
    USE.NAMES = FALSE
  ))

  ## Get number of coefficients for each parameter
  fnc <- as.integer(table(Xinfo$par))
  rnc <- if (is.null(Z)) rep.int(0L, p) else as.integer(table(Zinfo$par))

  ## Empty design matrices
  empty_dense_matrix <- matrix(integer(0L), nrow(X), 0L)
  empty_sparse_matrix <- sparseMatrix(
    i = integer(0L),
    j = integer(0L),
    x = integer(0L),
    dims = c(nrow(X), 0L)
  )

  l1 <- list(
    window = window,
    t = t,
    x = x,
    wl = wl,
    dow0 = dow0,
    X = X,
    Xs = if (sparse_X) X else empty_sparse_matrix,
    Xd = if (sparse_X) empty_dense_matrix else X,
    Z = if (is.null(Z)) empty_sparse_matrix else Z,
    Xinfo = Xinfo,
    Zinfo = Zinfo,
    fnc = fnc,
    rnc = rnc,
    curve_flag = get_flag("curve", curve),
    distr_flag = get_flag("distr", distr),
    excess_flag = 1L * excess,
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
#' When `par_init = NULL`, naive estimates of incidence model parameters
#' are obtained for each fitting window as follows:
#' \describe{
#' \item{`r`}{
#'   The slope of a linear model fit to `log1p(c(0, cumsum(x[-1])))`.
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
#'   fit to `log1p(c(0, cumsum(x[-1])))`.
#' }
#' \item{`tinfl`}{
#'   `max(t)`, assuming that the fitting window ends near the time
#'   of a peak in interval incidence.
#' }
#' \item{`K`}{
#'   `2*sum(x[-1])`, assuming that the fitting window ends near the
#'   time of a peak in interval incidence _and_ that interval incidence
#'   is roughly symmetric about the peak.
#' }
#' \item{`p`}{0.8}
#' \item{`a`, `b`, `nbdisp`, `w[123456]`}{1}
#' }
#' The naive estimates are log-transformed (all but `p`) or
#' logit-transformed (`p` only), and coefficients are obtained
#' from the within-group means (assuming zero-sum contrasts)
#' to form segments of the parameter object `beta`.
#'
#' When `par_init = NULL` and there are random effects,
#' `log_sd_b` and `b` are for simplicity initialized as zero vectors.
#'
#' @return
#' A list with elements:
#' \item{`beta`}{
#'   A numeric vector of length
#'   `sum(tmb_data$fnc * rowSums(tmb_data$fid))`
#'   listing initial values for the fixed effects coefficients.
#' }
#' \item{`log_sd_b`}{
#'   If there are random effects, then a numeric vector of length
#'   `sum(tmb_data$rid)` listing initial values for the log standard
#'   deviations of the random effects coefficients.
#'   Otherwise, `NA_real_`.
#' }
#' \item{`b`}{
#'   If there are random effects, then a numeric vector of length
#'   `sum(tmb_data$rnc * rowSums(tmb_data$rid))` listing initial
#'   values for the random effects coefficients.
#'   Otherwise, `NA_real_`.
#' }
#'
#' @examples
#' example("make_tmb_data", package = "epigrowthfit")
#' tmb_parameters <- epigrowthfit:::make_tmb_parameters(
#'   tmb_data = tmb_data,
#'   par_init = NULL
#' )
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis
make_tmb_parameters <- function(tmb_data, par_init) {
  any_re <- has_random(tmb_data)

  ## If user did not specify a parameter vector
  if (is.null(par_init)) {
    ts_split  <- split(as.data.frame(tmb_data[c("t", "x")]), tmb_data$window)
    fgrp_red <- tmb_data$fgrp[!duplicated(tmb_data$window), , drop = FALSE]

    get_log_r_log_c0 <- function(d) {
      h <- max(2, nrow(d) / 2)
      x <- d$t
      y <- log1p(c(0L, cumsum(d$x)))
      ab <- tryCatch(
        coef(lm(y ~ x, data = data.frame(x, y), subset = seq_len(h), na.action = na.omit)),
        error = function(e) c(0, 0.1)
      )
      c(log(ab[[2L]]), ab[[1L]])
    }
    get_log_tinfl <- function(d) log(max(d$t))
    get_log_K <- function(d) log(2) + log(sum(d$x[-1L], na.rm = TRUE))

    log_r_log_c0 <- vapply(ts_split, get_log_r_log_c0, c(0, 0))
    log_tinfl <- vapply(ts_split, get_log_tinfl, 0)
    log_K <- vapply(ts_split, get_log_K, 0)
    log_r  <- log_r_log_c0[1L, ]
    log_c0 <- log_r_log_c0[2L, ]
    p <- 0.8
    log_alpha <- switch(curve,
      subexponential = log_r + (1 - p) * log_c0,
      gompertz = log_r - log(log_K - log_c0),
      NA_real_
    )

    Y_init <- data.frame(
      log_r      = log_r,
      log_alpha  = log_alpha,
      log_c0     = log_c0,
      log_tinfl  = log_tinfl,
      log_K      = log_K,
      logit_p    = qlogis(p)
    )
    Y_init[sprintf("log_%s", c("a", "b", "nbdisp", paste0("w", 1:6)))] <- 0

    pn <- get_par_names(tmb_data, link = TRUE)
    get_term_labels <- function(par_name) {
      rownames(tmb_data$fid)[tmb_data$fid[, par_name] == 1L]
    }
    tl <- lapply(pn, get_term_labels)
    make_beta_segment <- function(par_name, term_label) {
      if (tmb_data$fvar[term_label] == "(Intercept)") {
        ## Within-group means
        m <- unname(tapply(Y_init[[par_name]], fgrp_red[[term_label]], mean))
        ## Mean of means
        mm <- mean(m)
        return(c(mm, m[-length(m)] - mm))
      }
      rep.int(0, tmb_data$fnc[[term_label]])
    }
    beta <- unlist(Map(make_beta_segment,
      par_name = rep.int(pn, lengths(tl)),
      term_label = unlist(tl),
      USE.NAMES = FALSE
    ))
    if (any_re) {
      log_sd_b <- rep.int(0, sum(tmb_data$rid))
      b <- rep.int(0, sum(tmb_data$rnc * rowSums(tmb_data$rid)))
    }

  ## If user specified a full parameter vector
  } else {
    p <- c(
      beta = sum(tmb_data$fnc * rowSums(tmb_data$fid)),
      log_sd_b = if (any_re) sum(tmb_data$rid),
      b = if (any_re) sum(tmb_data$rnc * rowSums(tmb_data$rid))
    )
    stop_if_not(
      is.numeric(par_init),
      length(par_init) == sum(p),
      is.finite(par_init),
      m = sprintf("`par_init` must be a finite numeric vector of length %d.", sum(p))
    )
    names(par_init) <- NULL
    beta <- par_init[seq_len(p[1L])]
    if (any_re) {
      log_sd_b <- par_init[seq.int(to = sum(p[1:2]), length.out = p[2L])]
      b <- par_init[seq.int(to = sum(p), length.out = p[3L])]
    }
  }

  if (!any_re) {
    log_sd_b <- b <- NA_real_
  }
  list(beta = beta, log_sd_b = log_sd_b, b = b)
}

#' Create TMB infrastructure
#'
#' Gathers necessary components of a call to [TMB::MakeADFun()].
#'
#' @param frame_ts,frame_glmm,window
#'   Elements of the list output of [make_frames()].
#' @inheritParams egf
#'
#' @return
#' A list with elements `data`, `parameters`, `map`, `random`,
#' `DLL`, and `silent`.
#'
#' @seealso [make_tmb_data()], [make_tmb_parameters()]
#' @keywords internal
make_tmb_args <- function(frame_ts, frame_glmm, data, window,
                          curve, distr, excess, weekday, sparse_X,
                          par_init) {
  tmb_data <- make_tmb_data(
    frame_ts = frame_ts,
    frame_glmm = frame_glmm,
    data = data,
    window = window,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    sparse_X = sparse_X
  )
  tmb_parameters <- make_tmb_parameters(
    tmb_data = tmb_data,
    par_init = par_init
  )
  if (has_random(tmb_data)) {
    tmb_map <- list()
    tmb_random <- "b"
  } else {
    tmb_map <- list(log_sd_b = factor(NA), b = factor(NA))
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
#' Determine whether an object specifies random effects.
#'
#' @param object
#'   An `"egf"` object returned by [egf()] or a `"tmb_data"` object
#'   returned by [make_tmb_data()].
#'
#' @return
#' `TRUE` if `object` specifies random effects and `FALSE` otherwise.
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

#' Enumerate TMB-generated parameter names
#'
#' Enumerates the names of a parameter vector generated by
#' [TMB::MakeADFun()], making them unique.
#'
#' @param par
#'   A TMB-generated parameter vector or, more generally,
#'   a named numeric vector such that
#'   `all(names(par) %in% c("beta", "log_sd_b", "b"))`.
#'
#' @return
#' `par` with `names(par)` replaced with names of the form
#' `"beta[i]"`, `"log_sd_b[j]"`, or `"b[k]"` (see Examples).
#'
#' @examples
#' n <- 10L
#' par <- numeric(n)
#' names(par) <- sample(c("beta", "log_sd_b", "b"), n, replace = TRUE)
#' rename_par(par)
#'
#' @noRd
rename_par <- function(par) {
  for (s in c("beta", "log_sd_b", "b")) {
    i <- which(names(par) == s)
    names(par)[i] <- sprintf("%s[%d]", s, seq_along(i))
  }
  par
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
