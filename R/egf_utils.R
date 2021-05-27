#' Get nonlinear model parameter names
#'
#' Retrieves the names used internally for nonlinear model parameters.
#'
#' @param curve,excess,distr,weekday
#'   For the default method, atomic scalar flags specifying
#'   an incidence model; see [egf()]. Set to `NULL` or `NA`
#'   to retrieve or suppress, respectively, all names associated
#'   with an argument. For other methods, `curve` is an object
#'   returned by [egf()] or [make_tmb_data()].
#' @param link
#'   A logical scalar. If `TRUE`, then `"name"` is replaced
#'   with `"log(name)"` or `"logit(name)"` depending on the
#'   link function used internally for the parameter.
#' @param ...
#'   Arguments passed to methods by the generic function.
#'
#' @return
#' If `link = FALSE`, then a subset of `"r"`, `"alpha"`,
#' `"c0"`, `"tinfl"`, `"K"`, `"p"`, `"a"`, `"b"`, `"nbdisp"`,
#' and `paste0("w", 1:6)`.
#'
#' @examples
#' ## For a specific model
#' get_par_names(
#'   curve = "logistic",
#'   excess = FALSE,
#'   distr = "nbinom",
#'   weekday = 0L,
#'   link = FALSE
#' )
#'
#' ## Across all possible models
#' get_par_names(NULL, NULL, NULL, NULL, link = FALSE)
#'
#' ## With link
#' get_par_names(NULL, NULL, NULL, NULL, link = TRUE)
#'
#' @export
get_par_names <- function(curve, ...) {
  UseMethod("get_par_names", curve)
}

#' @rdname get_par_names
#' @export
get_par_names.default <- function(curve = NULL, excess = NULL,
                                  distr = NULL, weekday = NULL,
                                  link = TRUE, ...) {
  if (is.null(curve)) {
    pn <- c("r", "alpha", "c0", "tinfl", "K", "p", "a")
  } else {
    pn <- switch(curve,
      exponential    = c("r", "c0"),
      logistic       = c("r", "tinfl", "K"),
      richards       = c("r", "tinfl", "K", "a"),
      subexponential = c("alpha", "c0", "p"),
      gompertz       = c("alpha", "c0", "K"),
      character(0L)
    )
  }
  if (is.null(excess) || isTRUE(excess)) {
    pn <- c(pn, "b")
  }
  if (is.null(distr) || isTRUE(distr == "nbinom")) {
    pn <- c(pn, "nbdisp")
  }
  if (is.null(weekday) || isTRUE(weekday > 0L)) {
    pn <- c(pn, paste0("w", 1:6))
  }
  if (link) {
    return(add_link_string(pn))
  }
  pn
}

#' @rdname get_par_names
#' @export
get_par_names.tmb_data <- function(curve, link = TRUE, ...) {
  pn <- levels(curve$X_info$par)
  if (link) {
    return(pn)
  }
  remove_link_string(pn)
}

#' @rdname get_par_names
#' @export
get_par_names.egf <- function(curve, link = TRUE, ...) {
  get_par_names(curve$tmb_args$data, link = link)
}

#' Manipulate nonlinear model parameter names
#'
#' Retrieve, add, and remove the link component of a string
#' giving the name of a nonlinear model parameter.
#'
#' @param s A character vector.
#'
#' @return
#' `get_link_string(s)` replaces valid `s[i]` with
#' `"log"` or `"logit"` according to the link function
#' used internally for `s[i]`.
#'
#' `add_link_string(s)` replaces valid `s[i]` with
#' `sprintf("%s(%s)", get_link_string(s[i]), s[i])`.
#'
#' `remove_link_string(s)` replaces valid `s[i]` with
#' `sub("^(log|logit)\\((.*)\\)$", "\\2", s[i])`.
#'
#' In all cases, invalid `s[i]` are replaced with
#' `NA_character_`.
#'
#' @examples
#' # add_link_string("r")
#' # remove_link_string("log(r)")
#' # get_link_string("r")
#' # get_link_string("log(r)")
#'
#' @name get_link_string
#' @keywords internal
NULL

#' @rdname get_link_string
get_link_string <- function(s) {
  ok <- s %in% get_par_names(curve = NULL, link = FALSE)
  s[ok] <- replace(rep_len("log", sum(ok)), s[ok] == "p", "logit")
  if (any(!ok)) { # avoids infinite recursion with `get_par_names()`
    ok_also <- !ok & s %in% get_par_names(curve = NULL, link = TRUE)
    s[ok_also] <- sub("^(log|logit)\\((.*)\\)$", "\\1", s[ok_also])
    s[!(ok | ok_also)] <- NA_character_
  }
  s
}

#' @rdname get_link_string
add_link_string <- function(s) {
  ok <- s %in% get_par_names(curve = NULL, link = FALSE)
  s[ok] <- sprintf("%s(%s)", get_link_string(s[ok]), s[ok])
  s[!ok] <- NA_character_
  s
}

#' @rdname get_link_string
remove_link_string <- function(s) {
  ok <- s %in% get_par_names(curve = NULL, link = TRUE)
  s[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\2", s[ok])
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
#' # get_link("log")
#' # get_inverse_link("logit")
#' # get_link("logit")
#' # get_inverse_link("log")
#'
#' @name get_link
#' @keywords internal
NULL

#' @rdname get_link
#' @importFrom stats qlogis
get_link <- function(s) {
  switch(s,
    log = log,
    logit = function(p) qlogis(p),
    stop("Link not implemented.")
  )
}

#' @rdname get_link
#' @importFrom stats plogis
get_inverse_link <- function(s) {
  switch(s,
    log = exp,
    logit = function(q) plogis(q),
    stop("Link not implemented.")
  )
}

#' Construct model frames
#'
#' Constructs time series and mixed effects model frames
#' to be used by [egf()], while completing a battery of
#' checks on the supplied formulae and data.
#'
#' @inheritParams egf
#'
#' @return
#' A list with elements `endpoints`, `frame_ts`, `frame_par`,
#' and `frame_append`. See descriptions in [egf()].
#'
#' @export
#' @importFrom stats terms model.frame na.pass as.formula
make_frames <- function(formula_ts, formula_par,
                        data_ts, data_par,
                        endpoints, origin,
                        curve, excess, distr, weekday,
                        na_action, init, append) {
  ### Check time series formula ###########################

  ## Test for formula of the form `x ~ time` or `x ~ time | ts`
  stop_if_not(
    inherits(formula_ts, "formula"),
    length(formula_ts) == 3L,
    m = "`formula_ts` must be a formula with a response."
  )
  formula_ts <- expand_terms(formula_ts)
  a <- attributes(terms(formula_ts))
  stop_if_not(
    a$intercept == 1L,
    is.null(a$offset),
    length(a$term.labels) == 1L,
    m = paste0(
      "To be parsed correctly, `formula_ts` must have\n",
      "an intercept, no offsets, and exactly one term."
    )
  )


  ### Check mixed effects formulae ########################

  ## Nonlinear model parameter names
  pn <- get_par_names(curve = curve, excess = excess,
                      distr = distr, weekday = weekday, link = TRUE)
  p <- length(pn)

  ## If a formula is supplied, then recycle to a list of length `p`
  if (inherits(formula_par, "formula") && length(formula_par) == 2L) {
    formula_par <- rep_len(list(expand_terms(formula_par)), p)
    names(formula_par) <- pn
  ## Otherwise, test for a valid list of formulae and complete
  ## with the default `~1` to obtain a list of length `p`
  } else {
    stop_if_not(
      is.list(formula_par),
      length(formula_par) > 0L,
      vapply(formula_par, inherits, FALSE, "formula"),
      lengths(formula_par) == 3L,
      (nfp <- vapply(formula_par, function(x) deparse1(x[[2L]]), "")) %in% pn,
      m = paste0(
        "`formula_par` must be a formula of the form `~terms`\n",
        "or a list of formulae of the form `par ~ terms`,\n",
        "with `deparse1(par)` an element of\n",
        "`get_par_names(curve, excess, distr, weekday, link = TRUE)`."
      )
    )
    formula_par <- lapply(formula_par, `[`, -2L)
    names(formula_par) <- nfp
    formula_par[setdiff(pn, names(formula_par))] <- list(~1)
    formula_par <- lapply(formula_par[pn], expand_terms)
  }

  ## Warn about fixed effects formulae without an intercept
  ## when relying on default initial values of linear coefficients
  if (is.null(init)) {
    is_intercept_ok <- function(formula) {
      a <- attributes(terms(split_effects(formula)$fixed))
      a$intercept == 1L || length(a$term.labels) == 0L
    }
    warn_if_not(
      vapply(formula_par, is_intercept_ok, FALSE),
      m = paste0(
        "Default initial values for fixed effects coefficients\n",
        "are unreliable when the fixed effects model does not\n",
        "have an intercept. Consider setting `init` explicitly\n",
        "or including an intercept."
      )
    )
  }


  ### Check data ##########################################

  stop_if_not(
    vapply(list(data_ts, data_par, endpoints), inherits, FALSE, c("data.frame", "list", "environment")),
    m = paste0(
      "`data_ts`, `data_par`, and `endpoints`\n",
      "must be data frames, lists, or environments."
    )
  )


  ### Check time series variables #########################

  try_eval1 <- function(expr) {
    try(eval(expr, envir = data_ts, enclos = environment(formula_ts)))
  }

  ## Validate incidence variable
  cn <- deparse1(formula_ts[[2L]])
  cases <- try_eval1(formula_ts[[2L]])
  stop_if_not(
    is.numeric(cases),
    (n <- length(cases)) > 0L,
    cases[!is.na(cases)] >= 0,
    m = sprintf("`%s` must evaluate to a nonempty,\nnon-negative, numeric vector.", cn)
  )
  if (is.double(cases)) {
    ## Round observations to nearest integer,
    ## with warning if change exceeds tolerance
    cases[is.infinite(cases)] <- NA
    warn_if_not(
      all.equal(cases, z <- round(cases)),
      m = sprintf("Non-integer numeric elements of `%s`\nrounded to nearest integer.", cn)
    )
    cases <- z
  }

  ## If formula uses a `|` operator, then validate grouping variable
  ## on right hand side
  if (is_bar(formula_ts[[3L]])) {
    gn <- deparse1(formula_ts[[3L]][[3L]])
    ts <- try_eval1(formula_ts[[3L]][[3L]])
    stop_if_not(
      is.factor(ts),
      m = sprintf("`%s` must evaluate to a factor.", gn)
    )
    stop_if_not(
      length(ts) == n,
      m = sprintf("`length(%s)` must equal `length(%s)`.", gn, cn)
    )
    stop_if_not(
      !anyNA(ts <- droplevels(ts, exclude = NA)),
      m = sprintf("`%s` must not have missing values.", gn)
    )
  ## Otherwise, assume that there is only one time series
  } else {
    gn <- "ts"
    ts <- rep_len(factor(1), n)
    formula_ts[[3L]] <- call("|", formula_ts[[3L]], as.name(gn))
  }
  lts <- levels(ts)
  M <- length(lts) # tentative number of time series

  ## Validate time variable
  tn <- deparse1(formula_ts[[3L]][[2L]])
  time <- try_eval1(formula_ts[[3L]][[2L]])
  if (inherits(time, "Date")) {
    time <- julian(time, origin = origin)
  }
  stop_if_not(
    is.numeric(time),
    m = sprintf("`%s` must evaluate to a numeric or Date vector.", tn)
  )
  stop_if_not(
    length(time) == n,
    m = sprintf("`length(%s)` must equal `length(%s)`.", tn, cn)
  )
  stop_if_not(
    !anyNA(time),
    m = sprintf("`%s` must not have missing values.", tn)
  )
  if (weekday > 0L && is.double(time)) {
    ## Weekday effect estimation requires integer time points,
    ## corresponding to 23:59:59 local time
    stop_if_not(
      all.equal(time, z <- round(time)),
      m = sprintf("weekday > 0: `%s` must be an integer vector\n(in the sense of `all.equal(%s, round(%s))`).", tn, tn, tn)
    )
    time <- z
  }
  stop_if_not(
    tapply(time, ts, function(x) all(diff(x) > 0)),
    m = sprintf("`%s` must be increasing in each level of `%s`.", tn, gn)
  )


  ### Check fitting windows ###############################

  stop_if_not(
    c("start", "end") %in% names(endpoints),
    m = "`endpoints` must have variables `start` and `end`."
  )
  start <- endpoints$start
  end <- endpoints$end
  if (inherits(start, "Date")) {
    start <- julian(start, origin = origin)
  }
  if (inherits(end, "Date")) {
    end <- julian(end, origin = origin)
  }
  stop_if_not(
    is.numeric(start),
    is.numeric(end),
    m = "In `endpoints`: `start` and `end` must be numeric or Date vectors."
  )
  stop_if_not(
    (N <- length(start)) > 0L, # tentative number of fitting windows
    length(end) == N,
    m = "In `endpoints`: `start` and `end` must have equal, nonzero lengths."
  )
  stop_if_not(
    !anyNA(start),
    !anyNA(end),
    m = "In `endpoints`: `start` and `end` must not have missing values."
  )
  stop_if_not(
    all(endpoints$start < endpoints$end),
    m = "In `endpoints`: `start` must precede `end`."
  )
  if (M > 1L) {
    epts <- try(eval(formula_ts[[3L]][[3L]], endpoints, environment(formula_ts)))
    stop_if_not(
      is.factor(epts),
      m = sprintf("In `endpoints`: `%s` must evaluate to a factor.", gn)
    )
    stop_if_not(
      length(epts) == N,
      m = sprintf("In `endpoints`: `start` and `%s` must have equal lengths.", gn)
    )
    stop_if_not(
      !anyNA(epts),
      m = sprintf("In `endpoints`: `%s` must not have missing values.", gn)
    )
    epts <- factor(epts, levels = lts)
    stop_if_not(
      !anyNA(epts),
      m = sprintf("In `endpoints`: No time series data supplied\nfor at least one level of `%s`.", gn)
    )
  } else {
    epts <- rep_len(factor(lts), N)
  }


  ### Construct mixed effects model frames ################

  formula_par_no_bars <- lapply(formula_par, gsub_bar_plus)
  frame_par <- lapply(formula_par_no_bars, model.frame,
    data = data_par,
    na.action = na.pass,
    drop.unused.levels = TRUE
  )

  ## Test that all model frames have one row per fitting window.
  ## For formulae without variables, such as `~1`, `model.frame()`
  ## returns zero-length model frames with zero rows. These should
  ## be excluded from the test.
  is_empty <- lengths(frame_par) == 0L
  stop_if_not(
    vapply(frame_par[!is_empty], nrow, 0L) == N,
    m = "`formula_par` variables must have length\nequal to `nrow(endpoints)`."
  )

  ## Still prefer the zero-length model frames to have one row
  ## per fitting window
  frame_par[is_empty] <- list(data.frame(row.names = seq_len(N)))


  ### Check mixed effects variables #######################

  ## Validate variable classes
  get_var_names <- function(x) {
    vapply(attr(terms(as.formula(call("~", x))), "variables")[-1L], deparse1, "")
  }
  is_bar_ok <- function(bar, data) {
    v <- lapply(bar[-1L], get_var_names)
    all(vapply(data[v[[1L]]], is.numeric, FALSE),
        vapply(data[v[[2L]]], is.factor,  FALSE))
  }
  all_bars_ok <- function(formula, data) {
    bars <- split_effects(formula)$random
    all(vapply(bars, is_bar_ok, FALSE, data))
  }
  stop_if_not(
    mapply(all_bars_ok, formula = formula_par, data = frame_par),
    m = paste0(
      "`formula_par` variables on left and right hand\n",
      "sides of `|` must evaluate to numeric vectors\n",
      "and factors, respectively."
    )
  )
  frame_par <- lapply(frame_par, droplevels)

  ## Test non-integer numeric variables for `Inf`
  any_infinite <- function(d) {
    i <- vapply(d, is.double, FALSE)
    any(vapply(d[i], function(x) any(is.infinite(x)), FALSE))
  }
  stop_if_not(
    !vapply(frame_par, any_infinite, FALSE),
    m = "Numeric `formula_par` variables must be finite\n(not Inf or -Inf)."
  )


  ### Construct clean table of fitting windows ##########

  endpoints <- data.frame(ts = epts, start, end)

  ## Discard fitting windows without observations
  ## on all mixed effects variables
  cc <- do.call(complete.cases, frame_par)
  if (na_action[2L] == "fail" && !all(cc)) {
    stop("na_action = \"fail\": `formula_par` variables\nmust not have missing values.")
  }
  warn_if_not(
    all(cc),
    m = "Fitting windows with incomplete mixed effects data discarded."
  )
  stop_if_not(
    any(cc),
    m = "At least one fitting window must have complete mixed effects data."
  )
  endpoints <- endpoints[cc, , drop = FALSE]
  frame_par <- lapply(frame_par, `[`, cc, , drop = FALSE)
  Ncc <- nrow(endpoints) # updated number of fitting windows

  ## Order remaining fitting windows by time series
  ## and chronologically within time series
  o <- do.call(order, unname(endpoints))
  endpoints <- endpoints[o, , drop = FALSE]
  frame_par <- lapply(frame_par, `[`, o, , drop = FALSE)

  ## Assign labels
  endpoints$window <- gl(Ncc, 1L, labels = sprintf("window_%0*d", nchar(Ncc), seq_len(Ncc)))
  endpoints <- endpoints[c("ts", "window", "start", "end")]


  ### Resume checking fitting windows #####################

  ## Create a factor `window` such that `split(time, window)`
  ## splits `time` by fitting window
  make_window_segment <- function(time, endpoints) {
    n <- length(time)
    N <- nrow(endpoints)
    lw <- levels(endpoints$window)
    window <- rep_len(factor(NA, levels = lw), n)
    if (n == 0L || N == 0L) {
      return(window)
    }
    stop_if_not(
      endpoints$start[-1L] > endpoints$end[-N],
      m = "In `endpoints`: Intervals [start, end]\nmust be disjoint within time series."
    )
    f <- function(t0, t1) which(time >= t0 & time <= t1)
    il <- Map(f, t0 = endpoints$start, t1 = endpoints$end)
    stop_if_not(
      lengths(il) >= 2L,
      m = "In `endpoints`: Intervals [start, end]\nmust contain at least two time points."
    )
    window[unlist(il)] <- rep.int(endpoints$window, lengths(il))
    window
  }
  window_split <- Map(make_window_segment,
    time  = split(time, ts),
    endpoints = split(endpoints, endpoints$ts)
  )
  window <- unsplit(window_split, ts)

  if (na_action[1L] == "fail") {
    stop_if_not(
      !tapply(cases, window, function(x) anyNA(x[-1L])),
      m = sprintf("na_action = \"fail\": `%s` has missing values\nin at least one fitting window.", cn)
    )
  } else {
    stop_if_not(
      tapply(cases, window, function(x) sum(!is.na(x[-1L])) > 0L),
      m = sprintf("Insufficient data: `%s` has zero observations\nin at least one fitting window.", cn)
    )
  }
  if (weekday > 0L) {
    stop_if_not(
      tapply(time, window, function(x) all(diff(x) == 1)),
      n = sprintf("weekday > 0: `%s` must have 1-day spacing\nin all fitting windows.", tn)
    )
  }

  ## Shrink fitting windows to range of time points
  endpoints$start <- c(tapply(time, window, min))
  endpoints$end <- c(tapply(time, window, max))


  ### Construct time series model frame ###################

  frame_ts <- data.frame(ts, window, time, x = cases)

  if (M > 1L) {
    ## Omit rows from time series without fitting windows
    endpoints$ts <- droplevels(endpoints$ts)
    frame_ts$ts <- factor(frame_ts$ts, levels = levels(endpoints$ts))
    frame_ts <- frame_ts[!is.na(frame_ts$ts), , drop = FALSE]

    ## Order remaining rows by time series
    ## NB: Chronological order within time series already checked
    frame_ts <- frame_ts[order(frame_ts$ts), , drop = FALSE]
  }


  ### Preserve user-specified `data_par` variables ########

  if (is.data.frame(data_par) && nrow(data_par) == N) {
    if (is.null(append)) {
      append <- seq_along(data_par)
    } else {
      append <- append_to_index(append, data_par, parent.frame())
    }
    frame_append <- data_par[which(cc)[o], append, drop = FALSE]
  } else {
    frame_append <- data.frame(row.names = seq_len(Ncc))
  }


  ### Set attributes ######################################

  row.names(endpoints) <- NULL
  row.names(frame_ts) <- NULL
  frame_par <- lapply(frame_par, `row.names<-`, NULL)
  row.names(frame_append) <- NULL

  attr(frame_ts, "terms") <- terms(formula_ts)
  frame_par <- Map(`attr<-`, frame_par, "terms", lapply(formula_par, terms))

  attr(endpoints, "origin") <- origin
  attr(frame_ts, "origin") <- origin

  list(
    endpoints = endpoints,
    frame_ts = frame_ts,
    frame_par = frame_par,
    frame_append = frame_append
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
#'   For `make_Z()`, a random effects term of the form `(terms | g)`.
#' @param frame
#'   A model frame (i.e., a data frame with an appropriate `terms`
#'   attribute) containing the variables used in `x`.
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
#'   `x = (foo | g)` does not use factors. `assign` indexes
#'   `labels(terms(~foo))` for `x = ~foo` and `x = (foo | g)`.
#' }
#' \item{`group`}{
#'   (`make_Z()` only.) For `x = (foo | g)`, a factor of length
#'   `ncol(Z)` with levels `levels(g)`, useful for splitting `Z`
#'   into group-specific submatrices.
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
  if (ncol(X) == 0L) {
    return(X)
  }
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
  gn <- vapply(split_interaction(x[[3L]]), deparse1, "")
  g <- interaction(frame[gn], drop = TRUE, sep = ":", lex.order = FALSE)
  J <- t(as(g, Class = "sparseMatrix"))
  Z <- t(KhatriRao(t(J), t(X)))
  ## For consistency, we want `model.matrix()`-style names
  ## for the group levels
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
#' A utility allowing columns of combined design matrices
#' and elements of combined parameter vectors to be split
#' by nonlinear model parameter, formula term, group, and
#' group level.
#'
#' @param xl
#'   A named list of fixed effects formulae of the form
#'   `~terms`, or a named list of random effects terms of
#'   the form `(terms | g)`. `names(xl)` must specify the
#'   respective nonlinear model parameters.
#' @param ml
#'   A list of `X` or `Z` matrices obtained by applying
#'   [make_X()] or [make_Z()] to the to the elements of `xl`.
#'
#' @return
#' A data frame with rows corresponding to columns of
#' `do.call(cbind, ml)`, and variables:
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
#' \item{`level`}{
#'   Level of the interaction specified by `group`.
#'   `NA` if not applicable.
#' }
#' \item{`colname`}{
#'   Column name.
#' }
#'
#' @keywords internal
#' @importFrom stats terms as.formula
make_XZ_info <- function(xl, ml) {
  if (length(xl) == 0L) {
    cn <- c("par", "term", "group", "level", "colname")
    m <- matrix(character(0L), ncol = 5L, dimnames = list(NULL, cn))
    return(as.data.frame(m, stringsAsFactors = TRUE))
  }
  ml_ncol <- vapply(ml, ncol, 0L)
  if (inherits(xl[[1L]], "formula")) { # X case
    group <- NA_character_
    level <- NA_character_
  } else { # Z case
    bar_lhs <- lapply(xl, `[[`, 2L)
    bar_rhs <- lapply(xl, `[[`, 3L)
    xl <- lapply(bar_lhs, function(x) as.formula(call("~", x)))
    group <- rep(vapply(bar_rhs, deparse1, ""), ml_ncol)
    level <- unlist(lapply(ml, attr, "group"))
  }
  get_term_labels <- function(formula) {
    labels(terms(formula))
  }
  rep_term_labels <- function(term_labels, assign) {
    c("(Intercept)", term_labels)[assign + 1L]
  }
  term <- unlist(Map(rep_term_labels,
    term_labels = lapply(xl, get_term_labels),
    assign = lapply(ml, attr, "assign")
  ))
  par <- rep(names(xl), ml_ncol)
  colname <- unlist(lapply(ml, colnames))
  data.frame(par, term, group, level, colname, row.names = NULL, stringsAsFactors = TRUE)
}

#' Construct data objects for C++ template
#'
#' Gathers in a list data objects to be passed to the C++ template
#' via TMB's `DATA_` macros. See [TMB::MakeADFun()].
#'
#' @inheritParams make_tmb_args
#'
#' @return
#' \[Below,
#' `n = sum(!is.na(frame_ts$window))`
#' denotes the total number of time points belonging a fitting window,
#' `N = nlevels(frame_ts$window)`
#' denotes the number of fitting windows, and
#' `p = length(frame_par)`
#' denotes the number of nonlinear model parameters.\]
#'
#' A `"tmb_data"` object, which is a list with elements:
#' \item{`t`}{
#'   A numeric vector of length `n` giving time as a number of days
#'   since the earliest time point in the current fitting window.
#' }
#' \item{`t_seg_len`}{
#'   An integer vector of length `N` specifying the length of each
#'   fitting window as a number of time points.
#' }
#' \item{`x`}{
#'   A numeric vector of length `n-N` giving incidence in each
#'   fitting window. `x[i]` in window `k` is the number of cases
#'   observed from time `t[k+i-1]` to time `t[k+i]`.
#' }
#' \item{`dow`}{
#'   If `weekday > 0L`, then an integer vector of length `N`
#'   indicating the first weekday in each fitting window, with
#'   value `i` in `0:6` mapping to the weekday `i` days after
#'   the reference weekday specified by `weekday`. Otherwise,
#'   an integer vector of the same length filled with `-1`.
#' }
#' \item{`Yo`}{
#'   The offsets matrix in dense format, with `N` rows and
#'   `p` columns.
#' }
#' \item{`X`}{
#'   The fixed effects design matrix in sparse or dense format,
#'   depending on `sparse_X`, with `N` rows.
#' }
#' \item{`Z`}{
#'   The random effects design matrix in sparse format, with
#'   `N` rows. If there are no random effects, then `Z` is an
#'   empty sparse matrix.
#' }
#' \item{`Xs`, `Xd`}{
#'   If `sparse_X = TRUE`, then `Xs = X` and `Xd` is an empty dense
#'   matrix. Otherwise, `Xd = X` and `Xs` is an empty sparse matrix.
#' }
#' \item{`X_info`, `Z_info`}{
#'   Data frames with `ncol(X)` and `ncol(Z)` rows, respectively.
#'   Row `j` describes the coefficient associated with column `j`
#'   of `X` or `Z`. `Z_info` has additional factors `cor` and `vec`
#'   splitting coefficients by relation to a common covariance
#'   matrix and random vector, respectively.
#' }
#' \item{`beta_seg_len`, `b_seg_len`}{
#'   Integer vectors of length `p` counting the columns of `X` and `Z`,
#'   respectively, pertaining to a common nonlinear model parameter.
#' }
#' \item{`beta_seg_index`, `b_seg_index`}{
#'   Integer vectors of length `ncol(X)` and `ncol(Z)`, respectively,
#'   with values in `0:(p-1)`. These split the columns of `X` and `Z`
#'   by relation to a common nonlinear model parameter.
#' }
#' \item{`block_rows`, `block_cols`}{
#'   Integer vectors together giving the dimensions of each block of
#'   the random effects matrix.
#' }
#' \item{`curve_flag`, `excess_flag`, `distr_flag`, `weekday_flag`, `sparse_X_flag`, `trace_flag`}{
#'   Integer flags referencing `curve`, `excess`, `distr`, `weekday > 0`,
#'   `sparse_X`, and `trace`.
#' }
#' \item{`predict_flag`}{
#'   An integer flag set equal to 0 so that prediction code is not run.
#' }
#' Additional integer elements of the form `j_link_parameter`
#' (e.g., `j_log_r`) give the 0-index of nonlinear model parameter
#' names (e.g., `"log(r)"`) in `names(frame_par)`. The value `-1`
#' is used for parameters not belonging to the model being fit.
#'
#' @keywords internal
#' @importFrom stats formula terms model.matrix model.offset
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix sparse.model.matrix KhatriRao
make_tmb_data <- function(frame_ts, frame_par,
                          curve, excess, distr, weekday,
                          sparse_X, trace) {
  ## Nonlinear model parameter names
  pn <- names(frame_par)
  p <- length(frame_par)

  ## Time series model frame excluding rows not associated
  ## with a fitting window
  origin <- attr(frame_ts, "origin")
  frame_ts <- frame_ts[!is.na(frame_ts$window), , drop = FALSE]

  ## Fitting window length as a number of time points
  t_seg_len <- tabulate(frame_ts$window)
  N <- length(t_seg_len)

  ## Time as number of days since earliest time point
  firsts <- which(!duplicated(frame_ts$window))
  t <- frame_ts$time - rep.int(frame_ts$time[firsts], t_seg_len)

  ## Incidence without unused first elements
  x <- frame_ts$x[-firsts]

  if (weekday > 0L) {
    ## Weekday during 1-day interval starting at earliest time point
    ## NB: When `weekday > 0`, time points are integer numbers
    ##     of days since 23:59:59 on a reference date, so this
    ##     1-day interval is guaranteed to occur on a single weekday

    ## Date during that interval
    ## NB: Reason for adding 1 is that the earliest time points are
    ##     23:59:59 on dates `origin + time[firsts]`, so the 1-day
    ##     intervals that we seek occur on the following dates
    date_day1 <- origin + frame_ts$time[firsts] + 1
    ## Weekday during that interval, coded as an integer `i` in `0:6`
    ## NB: `i` maps to the weekday `i` days after the reference
    ##     weekday, which is the weekday `weekday` days after an
    ##     arbitrary Saturday, in this case `.Date(2)`
    weekday_day1 <- julian(date_day1, origin = .Date(2 + weekday)) %% 7
    dow <- as.integer(weekday_day1)
  } else {
    dow <- rep_len(-1L, N)
  }

  ## Fixed effects formula and list of random effects terms
  ## extracted from each mixed effects formula
  formula_par <- lapply(frame_par, function(x) formula(terms(x)))
  effects <- lapply(formula_par, split_effects)
  fixed <- lapply(effects, `[[`, "fixed")
  random <- lapply(effects, `[[`, "random")
  random <- Map(`names<-`, random, Map(rep_len, names(random), lengths(random)))

  ## Offsets matrix
  offsets <- lapply(frame_par, model.offset)
  offsets[vapply(offsets, is.null, FALSE)] <- list(rep_len(0, N))
  Yo <- do.call(cbind, offsets)

  ## Fixed effects design matrices
  Xl <- Map(make_X, x = fixed, frame = frame_par, sparse = sparse_X)
  X <- do.call(cbind, Xl)

  ## Random effects design matrices
  Zl <- Map(make_Z, x = do.call(c, unname(random)), frame = rep.int(frame_par, lengths(random)))
  Z <- do.call(cbind, Zl)

  ## Nonlinear model parameter, formula term, group (g in (r | g)),
  ## and group level associated with each column of each matrix
  ## FIXME: Preserve contrasts?
  X_info <- make_XZ_info(xl = fixed, ml = Xl)
  Z_info <- make_XZ_info(xl = do.call(c, unname(random)), ml = Zl)
  X_info$par <- factor(X_info$par, levels = pn)
  Z_info$par <- factor(Z_info$par, levels = pn)

  ## Random effects coefficients factored by relation
  ## to a correlation matrix and to a random vector
  Z_info$cor <- interaction(Z_info[c("term", "group")], drop = TRUE, sep = " | ", lex.order = TRUE)
  Z_info$vec <- interaction(Z_info[c("term", "group", "level")], drop = TRUE, lex.order = TRUE)
  levels(Z_info$vec) <- seq_along(levels(Z_info$vec))

  ## Permutation of random effects coefficients producing
  ## correct arrangement in random effects blocks, given
  ## column-major order
  o <- do.call(order, unname(Z_info[c("cor", "vec", "par")]))
  Z <- Z[, o, drop = FALSE]
  Z_info <- Z_info[o, , drop = FALSE]

  ## Number of coefficients for each nonlinear model parameter
  beta_seg_len <- as.integer(table(X_info$par))
  b_seg_len <- as.integer(table(Z_info$par))

  ## Coefficients factored by nonlinear model parameter
  beta_seg_index <- as.integer(X_info$par) - 1L
  b_seg_index <- as.integer(Z_info$par) - 1L

  ## Dimensions of random effects blocks whose column vectors
  ## are related by a covariance matrix
  block_rows <- as.integer(colSums(table(Z_info$par, Z_info$cor) > 0L))
  block_cols <- as.integer(rowSums(table(Z_info$cor, Z_info$vec) > 0L))

  ## Empty design matrices
  empty_dense_matrix <- `dim<-`(integer(0L), c(N, 0L))
  empty_sparse_matrix <- sparseMatrix(
    i = integer(0L),
    j = integer(0L),
    x = integer(0L),
    dims = c(N, 0L)
  )

  l1 <- list(
    t = t,
    t_seg_len = t_seg_len,
    x = x,
    dow = dow,
    Yo = Yo,
    X = X,
    Z = if (is.null(Z)) empty_sparse_matrix else Z,
    Xs = if (sparse_X) X else empty_sparse_matrix,
    Xd = if (sparse_X) empty_dense_matrix else X,
    X_info = X_info,
    Z_info = Z_info,
    beta_seg_len = beta_seg_len,
    b_seg_len = b_seg_len,
    beta_seg_index = beta_seg_index,
    b_seg_index = b_seg_index,
    block_rows = block_rows,
    block_cols = block_cols,
    curve_flag = get_flag("curve", curve),
    excess_flag = as.integer(excess),
    distr_flag = get_flag("distr", distr),
    weekday_flag = as.integer(weekday > 0L),
    sparse_X_flag = as.integer(sparse_X),
    trace_flag = trace,
    predict_flag = 0L
  )
  pn0 <- get_par_names(curve = NULL, link = TRUE)
  l2 <- as.list(match(pn0, pn, 0L) - 1L)
  names(l2) <- sub("^(log|logit)\\((.*)\\)$", "j_\\1_\\2", pn0)
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
#' When `init = NULL`, naive estimates of nonlinear model parameters
#' are obtained for each fitting window as follows:
#' \describe{
#' \item{`r`}{
#'   The slope of a linear model fit to `log1p(cumsum(x)))`.
#' }
#' \item{`alpha`}{
#'   `r*c0^(1-p)` if `curve = "subexponential"`,
#'   `r/log(K/c0)` if `curve = "gompertz"`.
#'   These are the values obtained by setting the per capita growth
#'   rate at time 0 in the subexponential and Gompertz models equal
#'   to `r`, substituting the naive estimates of `r`, `c0`, `K`, and
#'   `p`, and solving for `alpha`.
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
#' effects model (if one exists) is taken to be the mean of the
#' link scale estimates across fitting windows.
#' The remaining elements of `beta` take initial value 0.
#' All elements of `b` and `theta` take initial value 0,
#' so that, initially, all random vectors follow a standard
#' normal distribution.
#'
#' @return
#' A list with elements:
#' \item{`beta`}{
#'   A numeric vector of length `sum(tmb_data$beta_seg_len)`
#'   listing initial values for the fixed effects coefficients.
#' }
#' \item{`b`}{
#'   If there are random effects, then a numeric vector of
#'   length `sum(tmb_data$b_seg_len)` listing initial values
#'   for the (unit variance scale) random effects coefficients.
#'   If there are no random effects, then `b = NA_real_`.
#' }
#' \item{`theta`}{
#'   If there are random effects, then a numeric vector of length
#'   `sum(f(tmb_data$block_rows))`, where `f(n) = n*(n+1)/2`,
#'   specifying initial covariance matrices for all random vectors.
#'   Let `u` be a segment of `theta` of length `f(n)`,
#'   corresponding to a random vector `x` of length `n`.
#'   Then `u[1:n]` lists log standard deviations of the elements
#'   of `x`, and `u[-(1:n)]` lists (in row-major order) the
#'   subdiagonal elements of the unit lower triangular matrix `L`
#'   in the Cholesky factorization
#'   `S = A %*% t(A)`,
#'   `A = D^-0.5 %*% L`,
#'   `D = diag(diag(L %*% t(L)))`
#'   of the correlation matrix `S` of `x`.
#'   If there are no random effects, then `theta = NA_real_`.
#' }
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis terms
make_tmb_parameters <- function(tmb_data, frame_ts, frame_par,
                                curve, init, debug) {
  ## Lengths of parameter objects
  f <- function(n) n * (n + 1L) / 2L
  len <- c(
    beta  = sum(tmb_data$beta_seg_len),
    b     = sum(tmb_data$b_seg_len),
    theta = sum(f(tmb_data$block_rows))
  )

  ## If user does not specify a full parameter vector
  if (is.null(init)) {
    ## Initialize each parameter object to a vector of zeros
    init_split <- lapply(len, numeric)

    ## Get names of nonlinear model parameters whose mixed
    ## effects formula has an intercept
    has_intercept <- function(frame) {
      attr(terms(frame), "intercept") == 1L
    }
    pn1 <- names(frame_par)[vapply(frame_par, has_intercept, FALSE)]

    ## For each of these nonlinear model parameters, compute
    ## the mean over all fitting windows of the naive estimate,
    ## and assign the result to the the coefficient of `beta`
    ## corresponding to "(Intercept)".
    if (length(pn1) > 0L) {
      ## Time series split by fitting window
      window <- frame_ts$window[!is.na(frame_ts$window)]
      firsts <- which(!duplicated(window))
      tx_split <- split(data.frame(t = tmb_data$t[-firsts], x = tmb_data$x), window[-firsts])

      ## Functions for computing naive estimates
      ## given a time series segment
      get_r_c0 <- function(d) {
        n <- max(2, trunc(nrow(d) / 2))
        ab <- tryCatch(
          coef(lm(log1p(cumsum(x)) ~ t, data = d, subset = seq_len(n), na.action = na.omit)),
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

      ## Naive estimates for each fitting window
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
      i1 <- match(pn1, tmb_data$X_info$par, 0L)

      ## Assign means over windows
      init_split$beta[i1] <- colMeans(Y_init[pn1])
    }

  ## If user specifies a full parameter vector
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

  ## Replace unused parameter objects with NA_real_
  if (!has_random(tmb_data)) {
    init_split[c("b", "theta")] <- list(NA_real_)
  }

  ## Retain all naive estimates when debugging
  if (is.null(init) && debug) {
    attr(init_split, "Y_init") <- Y_init[names(frame_par)]
  }
  init_split
}

#' Create TMB infrastructure
#'
#' Gathers necessary components of a call to [TMB::MakeADFun()].
#'
#' @param frame_ts,frame_par
#'   Elements of the list output of [make_frames()].
#' @inheritParams egf
#'
#' @return
#' A list with elements `data`, `parameters`, `map`, `random`,
#' `DLL`, and `silent`.
#'
#' @seealso [make_tmb_data()], [make_tmb_parameters()]
#' @keywords internal
make_tmb_args <- function(frame_ts, frame_par,
                          curve, distr, excess, weekday,
                          sparse_X, init, debug, trace) {
  tmb_data <- make_tmb_data(
    frame_ts = frame_ts,
    frame_par = frame_par,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    sparse_X = sparse_X,
    trace = trace
  )
  tmb_parameters <- make_tmb_parameters(
    tmb_data = tmb_data,
    frame_ts = frame_ts,
    frame_par = frame_par,
    curve = curve,
    init = init,
    debug = debug
  )
  if (has_random(tmb_data)) {
    ## Nothing to fix, since all parameter objects are used
    tmb_map <- list()
    ## Declare that `b` contains random effects
    tmb_random <- "b"
  } else {
    ## Fix `b` and `theta` (in this case to NA_real_),
    ## since only `beta` is used
    tmb_map <- list(b = factor(NA), theta = factor(NA))
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
#' @keywords internal
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

#' Split `ADREPORT()`ed variables
#'
#' Extracts reported variables from an `"sdreport"` object
#' and returns them in a list.
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
#' @keywords internal
split_sdreport <- function(sdreport) {
  ssdr <- summary(sdreport, select = "report")
  colnames(ssdr) <- c("estimate", "se")
  lapply(split(as.data.frame(ssdr), rownames(ssdr)), `row.names<-`, NULL)
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

#' Construct combined model frame
#'
#' Joins in a single data frame all mixed effects model frames
#' and further variables specified via [egf()] argument `append`.
#'
#' @param object An `"egf"` object returned by [egf()].
#'
#' @details
#' If a variable name occurs in multiple mixed effects model frames,
#' then only one instance is retained. Except in unusual cases,
#' all instances of a variable name are identical, and no information
#' is lost.
#'
#' Since the data frames being combined each correspond rowwise to
#' `object$endpoints`, so does the result.
#'
#' @return
#' A data frame combining (in the sense of [cbind()]) all data frames
#' listed in `object$frame_par` and the data frame `object$frame_append`.
#'
#' @keywords internal
make_combined <- function(object) {
  stop_if_not(
    inherits(object, "egf"),
    m = "`object` must inherit from class \"egf\"."
  )
  frames <- c(unname(object$frame_par), list(object$frame_append))
  combined <- do.call(cbind, frames)
  combined[!duplicated(names(combined))]
}

do_wald <- function(estimate, se, level) {
  q <- qchisq(level, df = 1)
  n <- length(estimate)
  lu <- estimate + rep.int(sqrt(q) * c(-1, 1), c(n, n)) * se
  dim(lu) <- c(n, 2L)
  colnames(lu) <- c("lower", "upper")
  lu
}

apply_inverse_link <- function(x, g) {
  if (is.data.frame(x)) {
    f <- function(x, s) {
      x[] <- lapply(x, get_inverse_link(s))
      x
    }
  } else {
    f <- function(x, s) {
      get_inverse_link(s)(x)
    }
  }
  x_split <- split(x, g, drop = TRUE)
  fx_split <- Map(f, x = x_split, s = get_link_string(names(x_split)))
  unsplit(fx_split, g, drop = TRUE)
}
