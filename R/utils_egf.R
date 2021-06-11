#' Get model parameter names
#'
#' Retrieves the names used internally for nonlinear and dispersion
#' model parameters.
#'
#' @param object
#'   An \code{"\link{egf_model}"}, \code{"\link{egf}"},
#'   or \code{"\link[=make_tmb_data]{tmb_data}"} object
#'   specifying a model. Otherwise, \code{\link{NULL}}.
#' @param link
#'   A \link{logical} scalar. If \code{TRUE}, then \code{"name"}
#'   is replaced with \code{"log(name)"} or \code{"logit(name)"}
#'   depending on the link function used internally for the parameter.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{character} vector giving a subset of \code{"r"}, \code{"alpha"},
#' \code{"c0"}, \code{"tinfl"}, \code{"K"}, \code{"p"}, \code{"a"},
#' \code{"b"}, \code{"nbdisp"}, and \code{\link{paste0}("w", 1:6)},
#' possibly modified according to \code{link}.
#'
#' If \code{object = \link{NULL}}, then the complete set is returned.
#'
#' @examples
#' get_par_names(NULL, link = FALSE)
#' get_par_names(NULL, link = TRUE)
#'
#' model <- egf_model(
#'   curve = "logistic",
#'   excess = FALSE,
#'   family = "nbinom",
#'   day_of_week = 0L
#' )
#' get_par_names(model, link = FALSE)
#'
#' @export
get_par_names <- function(object, ...) {
  UseMethod("get_par_names", object)
}

#' @rdname get_par_names
#' @export
get_par_names.default <- function(object, link = TRUE, ...) {
  stop_if_not(
    is.null(object),
    m = "This method accepts `object = NULL` only."
  )
  pn <- c("r", "alpha", "c0", "tinfl", "K",
          "p", "a", "b", "nbdisp", paste0("w", 1:6))
  if (link) {
    return(string_add_link(pn))
  }
  pn
}

#' @rdname get_par_names
#' @export
get_par_names.egf_model <- function(object, link = TRUE, ...) {
  pn <- switch(object$curve,
    exponential    = c("r", "c0"),
    logistic       = c("r", "tinfl", "K"),
    richards       = c("r", "tinfl", "K", "a"),
    subexponential = c("alpha", "c0", "p"),
    gompertz       = c("alpha", "c0", "K")
  )
  if (object$excess) {
    pn <- c(pn, "b")
  }
  if (object$family == "nbinom") {
    pn <- c(pn, "nbdisp")
  }
  if (object$day_of_week > 0L) {
    pn <- c(pn, paste0("w", 1:6))
  }
  if (link) {
    return(string_add_link(pn))
  }
  pn
}

#' @rdname get_par_names
#' @export
get_par_names.egf <- function(object, link = TRUE, ...) {
  get_par_names(object$model, link = link)
}

#' @rdname get_par_names
#' @export
get_par_names.tmb_data <- function(object, link = TRUE, ...) {
  pn <- levels(object$X_info$par)
  if (link) {
    return(pn)
  }
  string_remove_link(pn)
}

#' Manipulate model parameter names
#'
#' Utilities for modifying strings \code{s} and
#' \code{fs = \link{sprintf}("\%s(\%s)", f, s)},
#' where \code{s} is the name used internally for
#' a nonlinear or dispersion model parameter and
#' \code{f} is the name (either \code{"log"} or \code{"logit"})
#' of the link function used internally for that parameter.
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
#' is returned for invalid vector elements (see Examples).
#'
#' @examples
#' # string_get_link("r")
#' # string_add_link("r")
#' # string_remove_link("log(r)")
#' # string_extract_link("log(r)")
#' # string_extract_link("log(r)", "invalid", "input")
#'
#' @name string_get_link
#' @keywords internal
NULL

#' @rdname string_get_link
string_get_link <- function(s) {
  ok <- s %in% get_par_names(NULL, link = FALSE)
  s[ok] <- replace(rep_len("log", sum(ok)), s[ok] == "p", "logit")
  s[!ok] <- NA_character_
  s
}

#' @rdname string_get_link
string_add_link <- function(s) {
  ok <- s %in% get_par_names(NULL, link = FALSE)
  s[ok] <- sprintf("%s(%s)", string_get_link(s[ok]), s[ok])
  s[!ok] <- NA_character_
  s
}

#' @rdname string_get_link
string_remove_link <- function(fs) {
  ok <- fs %in% get_par_names(NULL, link = TRUE)
  fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\2", fs[ok])
  fs[!ok] <- NA_character_
  fs
}

#' @rdname string_get_link
string_extract_link <- function(fs) {
  ok <- fs %in% get_par_names(NULL, link = TRUE)
  fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\1", fs[ok])
  fs[!ok] <- NA_character_
  fs
}

#' Get link and inverse link functions
#'
#' Retrieve the link function named by a string, or its inverse.
#'
#' @param f
#'   A character string naming a link function,
#'   either \code{"log"} or \code{"logit"}.
#' @param inverse
#'   A logical scalar. If \code{TRUE}, then the inverse is returned.
#'
#' @return
#' A \link{function}.
#'
#' @examples
#' # match_link("log")
#' # match_link("log", inverse = TRUE)
#' # match_link("logit")
#' # match_link("logit", inverse = TRUE)
#'
#' @name match_link
#' @keywords internal
#' @importFrom stats plogis qlogis
match_link <- function(f, inverse = FALSE) {
  if (inverse) {
    switch(f,
           log = exp,
           logit = function(q) plogis(q),
           stop("Link not implemented.")
    )
  } else {
    switch(f,
           log = log,
           logit = function(p) qlogis(p),
           stop("Link not implemented.")
    )
  }
}

#' Construct model frames
#'
#' Constructs time series and mixed effects model frames to be used
#' by \code{\link{egf}}, while completing a battery of checks on the
#' supplied formulae and data.
#'
#' @inheritParams egf
#'
#' @return
#' A list with elements \code{frame}, \code{frame_par}, \code{frame_append},
#' and \code{endpoints}. See descriptions under \code{\link{egf}}.
#'
#' @keywords internal
#' @importFrom stats terms model.frame na.pass as.formula complete.cases
make_frames <- function(formula, formula_par,
                        data, data_par,
                        subset, subset_par,
                        na_action, na_action_par,
                        endpoints, origin,
                        model, init, append) {
  stop_if_not(
    vapply(list(data, data_par, endpoints), inherits, FALSE, c("data.frame", "list", "environment")),
    m = "`data`, `data_par`, and `endpoints` must be data frames, lists, or environments."
  )
  stop_if_not(
    inherits(origin, "Date"),
    length(origin) == 1L,
    !is.na(origin),
    m = "`origin` must be a Date vector of length 1."
  )


  ### Construction step ########################################################

  ### Time series model frame

  ## Check that `formula` has the form `x ~ time` or `x ~ time | ts`
  stop_if_not(
    inherits(formula, "formula"),
    m = "`formula` must be a formula."
  )
  formula <- simplify_terms(formula)
  a <- attributes(terms(formula))
  stop_if_not(
    a$response > 0L,
    a$intercept == 1L,
    is.null(a$offset),
    length(a$term.labels) == 1L,
    m = paste0(
      "To be parsed correctly, `formula` must have a response,\n",
      "an intercept, no offsets, and exactly one term.\n",
      "In particular, after expansion of all formula operators,\n",
      "`formula` should have the form `x ~ time` or `x ~ time | ts."
    )
  )

  ## Construct and evaluate call to `data.frame()`
  l <- list(quote(data.frame), time = formula[[3L]], x = formula[[2L]])
  is_bar <- function(x) is.call(x) && x[[1L]] == as.name("|")
  if (ib <- is_bar(formula[[3L]])) {
    l[c("time", "ts")] <- formula[[3L]][2:3]
  }
  frame <- eval(as.call(l), envir = data, enclos = environment(formula))
  if (!ib) {
    ## Assume that there is only one time series
    frame$ts <- rep_len(factor(1), nrow(frame))
  }

  ## Preserve deparsed `data.frame()` arguments for error messages
  nf <- lapply(l[-1L], deparse)


  ### Mixed effects model frames

  ## Nonlinear and dispersion model parameter names
  pn <- get_par_names(model, link = TRUE)
  p <- length(pn)

  ## If `formula_par` is a formula of the form `~terms`,
  ## then recycle to a list of length `p`
  if (inherits(formula_par, "formula") && length(formula_par) == 2L) {
    formula_par <- simplify_terms(formula_par)
    formula_par <- rep_len(list(formula_par), p)
    names(formula_par) <- pn
    ## Otherwise, test for a valid list of formulae
    ## of the form `par ~ terms`, and complete with
    ## the default `~1` to obtain a list of length `p`
  } else {
    stop_if_not(
      is.list(formula_par),
      length(formula_par) > 0L,
      vapply(formula_par, inherits, FALSE, "formula"),
      lengths(formula_par) == 3L,
      (nfp <- vapply(formula_par, function(x) deparse(x[[2L]]), "")) %in% pn,
      m = paste0(
        "`formula_par` must be a formula of the form `~terms`\n",
        "or a list of formulae of the form `par ~ terms`, with\n",
        "`deparse(par)` an element of `get_par_names(model, link = TRUE)`."
      )
    )
    formula_par <- lapply(formula_par, `[`, -2L)
    names(formula_par) <- nfp
    formula_par[setdiff(pn, names(formula_par))] <- list(~1)
    formula_par <- lapply(formula_par[pn], simplify_terms)
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
        "Default initial values for linear coefficients are not\n",
        "reliable for fixed effects models with an intercept.\n",
        "Consider setting `init` explicitly or including an intercept."
      )
    )
  }

  ## Replace `|` operators with `+` operators in formulae
  ## to enable use of `model.frame()` machinery
  formula_par_no_bars <- lapply(formula_par, gsub_bar_plus)
  frame_par <- lapply(formula_par_no_bars, model.frame,
                      data = data_par,
                      na.action = na.pass,
                      drop.unused.levels = TRUE
  )

  ## Test that all model frames have the same number of rows.
  ## For formulae without variables (e.g., `~1`), `model.frame()`
  ## returns data frames of zero length with zero rows.
  ## These can be excluded from the test and modified separately.
  nr <- vapply(frame_par, nrow, 0L)
  N <- max(nr)
  len <- lengths(frame_par)
  stop_if_not(
    nr[len > 0L] == N,
    m = paste0(
      "Mixed effects model frames constructed from `formula_par`\n",
      "and `data_par` must have the same number of rows."
    )
  )
  frame_par[len == 0L] <- list(data.frame(row.names = seq_len(N)))


  ### Table of fitting windows

  ## Construct and evaluate call to `data.frame()`
  l <- list(quote(data.frame), start = quote(start), end = quote(end))
  if (ib) {
    l$ts <- formula[[3L]][[3L]]
  }
  endpoints <- eval(as.call(l), envir = endpoints, enclos = environment(formula))
  if (!ib) {
    ## Assume that there is only one time series
    endpoints$ts <- rep_len(factor(1), nrow(endpoints))
  }

  ## Check for rowwise correspondence with mixed effects model frames
  stop_if_not(
    nrow(endpoints) == N,
    m = "`endpoints` and mixed effects model frames must have the same number of rows."
  )

  ## Preserve deparsed `data.frame()` arguments for error messages
  ne <- lapply(l[-1L], deparse1)


  ### Appended stuff

  if (is.data.frame(data_par) && nrow(data_par) == N) {
    if (is.null(append)) {
      append <- seq_along(data_par)
    } else {
      append <- append_to_index(append, data = data_par, enclos = baseenv())
    }
    frame_append <- data_par[append]
  } else {
    frame_append <- data.frame(row.names = seq_len(N))
  }


  ### Subsetting step ##########################################################

  if (is.data.frame(data)) {
    subset <- subset_to_index(subset, data = data, enclos = parent.frame())
    frame <- frame[subset, , drop = FALSE]
  }
  if (is.data.frame(data_par)) {
    subset_par <- subset_to_index(subset_par, data = data_par, enclos = parent.frame())
    endpoints <- endpoints[subset_par, , drop = FALSE]
    frame_par <- lapply(frame_par, `[`, subset_par, , drop = FALSE)
    frame_append <- frame_append[subset_par, , drop = FALSE]
  }


  ### Validation step ##########################################################

  ### Time series model frame

  stop_if_not(
    nrow(frame) > 0L,
    m = "Time series model frame has 0 rows (after subsetting with `subset`)."
  )

  ## Time variable
  if (inherits(frame$time, "Date")) {
    frame$time <- julian(frame$time, origin = origin)
  }
  stop_if_not(
    is.numeric(frame$time),
    is.finite(frame$time),
    m = sprintf("`%s` must be a finite numeric or Date vector.", nf$time)
  )
  if (model$day_of_week > 0L && is.double(frame$time)) {
    ## Day of week effect estimation requires integer time points,
    ## corresponding to midnights in local time
    stop_if_not(
      all.equal(frame$time, z <- round(frame$time)),
      m = sprintf("model$day_of_week > 0: `%s` must be an integer or Date vector.", nf$time)
    )
    frame$time <- z
  }

  ## Incidence variable
  stop_if_not(
    is.numeric(frame$x),
    frame$x[!is.na(frame$x)] >= 0,
    m = sprintf("`%s` must be a non-negative numeric vector.", nf$x)
  )
  if (is.double(frame$x)) {
    if (any(is_NaN_or_Inf <- is.nan(frame$x) | is.infinite(frame$x))) {
      warning(sprintf("NaN, Inf, and -Inf in `%s` replaced with NA.", nf$x))
      frame$x[is_NaN_or_Inf] <- NA
    }
    if (!isTRUE(all.equal(frame$x, z <- round(frame$x)))) {
      warning(sprintf("Non-integer numeric elements of `%s` rounded to nearest integer.", nf$x))
    }
    frame$x <- z
  }

  ## Grouping variable
  if (ib) {
    stop_if_not(
      is.factor(frame$ts),
      m = sprintf("`%s` must be a factor.", nf$ts)
    )
    frame$ts <- factor(frame$ts, exclude = NA)
    stop_if_not(
      !anyNA(frame$ts),
      m = sprintf("`%s` must not have missing values.", nf$ts)
    )
  }


  ### Mixed effects model frames

  stop_if_not(
    nrow(frame_par[[1L]]) > 0L,
    m = "Mixed effects model frames have 0 rows (after subsetting with `subset_par`)."
  )

  ## Test data frames for Inf and -Inf
  any_infinite <- function(d) {
    i <- vapply(d, is.double, FALSE)
    any(vapply(d[i], function(x) any(is.infinite(x)), FALSE))
  }
  stop_if_not(
    !vapply(frame_par, any_infinite, FALSE),
    m = "Numeric `formula_par` variables must not contain Inf or -Inf."
  )

  ## Require that `x` and `g` in random effects terms `(x | g)`
  ## are expressions involving only numeric vectors and factors,
  ## respectively
  get_var_names <- function(x) {
    vapply(attr(terms(as.formula(call("~", x))), "variables")[-1L], deparse, "")
  }
  is_bar_ok <- function(bar, data) {
    v <- lapply(bar[-1L], get_var_names)
    all(vapply(data[v[[1L]]], is.numeric, FALSE),
        vapply(data[v[[2L]]], is.factor,  FALSE))
  }
  all_bars_ok <- function(formula, data) {
    bars <- split_effects(formula)$random
    all(vapply(bars, is_bar_ok, FALSE, data = data))
  }
  stop_if_not(
    mapply(all_bars_ok, formula = formula_par, data = frame_par),
    m = paste0(
      "`formula_par` variables on left and right hand sides\n",
      "of `|` must be numeric vectors and factors, respectively."
    )
  )


  ### Table of fitting windows

  ## Endpoints
  if (inherits(endpoints$start, "Date")) {
    endpoints$start <- julian(endpoints$start, origin = origin)
  }
  if (inherits(endpoints$end, "Date")) {
    endpoints$end <- julian(endpoints$end, origin = origin)
  }
  stop_if_not(
    is.numeric(endpoints$start),
    is.numeric(endpoints$end),
    is.finite(endpoints$start),
    is.finite(endpoints$end),
    m = sprintf("In `endpoints`: `%s` and `%s` must be finite numeric or Date vectors.", ne$start, ne$end)
  )
  stop_if_not(
    endpoints$start < endpoints$end,
    m = sprintf("In `endpoints`: `%s` must precede `%s`.", ne$start, ne$end)
  )

  ## Grouping variable
  if (ib) {
    stop_if_not(
      is.factor(endpoints$ts),
      m = sprintf("In `endpoints`: `%s` must be a factor.", ne$ts)
    )
    endpoints$ts <- factor(endpoints$ts, exclude = NA)
    stop_if_not(
      !anyNA(endpoints$ts),
      m = sprintf("In `endpoints`: `%s` must not have missing values.", ne$ts)
    )
    sdl <- setdiff(levels(endpoints$ts), levels(frame$ts))
    stop_if_not(
      length(sdl) == 0L,
      m = paste0(
        sprintf("These levels of `%s` are missing corresponding time series data:\n", ne$ts),
        paste(sdl, collapse = "\n")
      )
    )
    stop_if_not(
      tapply(frame$time, frame$ts, function(x) all(diff(x) > 0)),
      m = sprintf("`%s` must be increasing in each level of `%s`.", nf$time, nf$ts)
    )
  } else {
    stop_if_not(
      all(diff(frame$time) > 0),
      m = sprintf("`%s` must be increasing.", nf$time)
    )
  }


  ### Another subsetting step ##################################################

  ## Discard fitting windows without observations
  ## on all mixed effects variables
  if (any(len > 0L)) {
    cc <- do.call(complete.cases, frame_par[len > 0L])
    if (!all(cc)) {
      if (na_action_par == "fail") {
        stop("na_action_par = \"fail\": `formula_par` variables must not have missing values.")
      }
      if (!any(cc)) {
        stop("At least one fitting window must have complete mixed effects data.")
      }
      warning("Fitting windows with incomplete mixed effects data discarded.")
      endpoints <- endpoints[cc, , drop = FALSE]
      frame_par <- lapply(frame_par, `[`, cc, , drop = FALSE)
      frame_append <- frame_append[cc, , drop = FALSE]
      endpoints$ts <- factor(endpoints$ts)
    }
  }

  ## Discard time series without fitting windows
  frame$ts <- factor(frame$ts, levels = levels(endpoints$ts))
  frame <- frame[!is.na(frame$ts), , drop = FALSE]

  ## Order everything by time series and chronologically within time series
  o1 <- do.call(order, unname(frame[c("ts", "time")]))
  frame <- frame[o1, , drop = FALSE]
  o2 <- do.call(order, unname(endpoints[c("ts", "start", "end")]))
  endpoints <- endpoints[o2, , drop = FALSE]
  frame_par <- lapply(frame_par, `[`, o2, , drop = FALSE)
  frame_append <- frame_append[o2, , drop = FALSE]


  ### Labeling step ############################################################

  ## Enumerate fitting windows as currently ordered
  N <- nrow(endpoints)
  endpoints$window <- gl(N, 1L, labels = sprintf("window_%0*d", nchar(N), seq_len(N)))

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
      m = "In `endpoints`: Intervals [start, end] must be disjoint within time series."
    )
    f <- function(t0, t1) which(time >= t0 & time <= t1)
    il <- Map(f, t0 = endpoints$start, t1 = endpoints$end)
    nl <- lengths(il)
    stop_if_not(
      nl >= 2L,
      m = "In `endpoints`: Intervals [start, end] must contain at least two time points."
    )
    window[unlist(il, FALSE, FALSE)] <- rep.int(endpoints$window, nl)
    window
  }
  window_split <- Map(make_window_segment,
                      time = split(frame$time, frame$ts),
                      endpoints = split(endpoints, endpoints$ts)
  )
  frame$window <- unsplit(window_split, frame$ts)

  ## Contract endpoints in table to range of time points between
  endpoints$start <- c(tapply(frame$time, frame$window, min))
  endpoints$end <- c(tapply(frame$time, frame$window, max))


  ### Another validation step ##################################################

  if (na_action == "fail") {
    stop_if_not(
      !tapply(frame$x, frame$window, function(x) anyNA(x[-1L])),
      m = sprintf("na_action = \"fail\": `%s` has missing values in at least one fitting window.", nf$x)
    )
  } else {
    stop_if_not(
      tapply(frame$x, frame$window, function(x) sum(!is.na(x[-1L])) > 0L),
      m = sprintf("`%s` has zero observations in at least one fitting window.", nf$x)
    )
  }
  if (model$day_of_week > 0L) {
    stop_if_not(
      tapply(frame$time, frame$window, function(x) all(diff(x) == 1)),
      n = sprintf("model$day_of_week > 0: `%s` must have 1-day spacing in all fitting windows.", nf$time)
    )
  }


  ### Cleaning up ##############################################################

  frame <- frame[c("ts", "window", "time", "x")]
  endpoints <- endpoints[c("ts", "window", "start", "end")]

  frame_par <- lapply(frame_par, droplevels)
  frame_append <- droplevels(frame_append)

  row.names(frame) <- NULL
  frame_par <- lapply(frame_par, `row.names<-`, NULL)
  row.names(frame_append) <- NULL
  row.names(endpoints) <- NULL

  attr(frame, "terms") <- terms(formula)
  frame_par <- Map(`attr<-`, frame_par, "terms", lapply(formula_par, terms))

  attr(frame, "origin") <- origin
  attr(endpoints, "origin") <- origin

  list(
    frame = frame,
    frame_par = frame_par,
    frame_append = frame_append,
    endpoints = endpoints
  )
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
#'   of the form \code{(tt | g)}. Here, \code{tt} is an expression
#'   composed of potentially many terms, while \code{g}
#'   is an unevaluated factor, possibly an interaction.
#' @param frame
#'   A \link[=model.frame]{model frame} listing the variables used
#'   in \code{x}.
#' @param sparse
#'   A \link{logical} scalar. If \code{TRUE}, then the design matrix
#'   is returned in sparse format.
#'
#' @details
#' \code{make_X(x, frame, sparse)} constructs an \code{X} matrix
#' by evaluating \code{\link{model.matrix}(x, frame)} or
#' \code{\link[Matrix]{sparse.model.matrix}(x, frame)}
#' (depending on \code{sparse}) and deleting from the result
#' columns containing only zeros.
#'
#' \code{make_Z(x, frame)} constructs a \code{Z} matrix
#' following the steps outlined in \code{\link{vignette}("lmer", "lme4")}.
#' It uses \code{make_X} to construct the so-called raw model matrix
#' from \code{x[[2L]]}. The result is always sparse.
#'
#' @return
#' A matrix with \link{attributes}:
#' \item{contrasts}{
#'   See \code{\link{model.matrix}}.
#'   \code{contrasts} is absent if the expression \code{tt}
#'   in \code{x = ~tt} or \code{x = (tt | g)} does not involve
#'   \link{factor}s.
#' }
#' \item{assign}{
#'   See \code{\link{model.matrix}}.
#'   \code{assign} indexes \code{\link{labels}(\link{terms}(~tt))}
#'   for the expression `tt` in \code{x = ~tt} or \code{x = (tt | g)}.
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
#' A utility to link each column of combined design matrices
#' (and similarly each element of combined parameter vectors)
#' to the corresponding nonlinear and dispersion model parameter,
#' mixed effects formula term, grouping variable, and group level.
#'
#' @param xl
#'   A named \link{list} of \link{formula}e of the form \code{~tt},
#'   or a named list of \link{call}s to binary operator \code{`|`}
#'   of the form \code{(tt | g)}. \code{\link{names}(xl)} must specify
#'   nonlinear or dispersion model parameters.
#' @param ml
#'   A \link{list} of \code{X} or \code{Z} matrices obtained by
#'   applying \code{\link{make_X}} or \code{\link{make_Z}} to the
#'   elements of \code{xl}.
#'
#' @return
#' A \link[=data.frame]{data frame} with rows corresponding to
#' columns of the combined design matrix \code{do.call(cbind, ml)},
#' and variables:
#' \item{par}{
#'   Nonlinear or dispersion model parameter.
#' }
#' \item{term}{
#'   Deparsed term from right hand side of \code{`~`} or left hand side
#'   of \code{`|`}. Otherwise, \code{"(Intercept)"}.}
#' \item{group}{
#'   Deparsed term from right hand side of \code{`|`}.
#'   \code{\link{NA}} if not applicable.
#' }
#' \item{level}{
#'   Level of the interaction specified by \code{group}.
#'   \code{\link{NA}} if not applicable.
#' }
#' \item{colname}{
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
    group <- rep.int(vapply(bar_rhs, deparse, ""), ml_ncol)
    level <- unlist(lapply(ml, attr, "group"), FALSE, FALSE)
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
  ), FALSE, FALSE)
  par <- rep.int(names(xl), ml_ncol)
  colname <- unlist(lapply(ml, colnames), FALSE, FALSE)
  data.frame(par, term, group, level, colname, row.names = NULL, stringsAsFactors = TRUE)
}

#' Create TMB infrastructure
#'
#' Gathers necessary components of a call to \code{\link[TMB]{MakeADFun}}.
#'
#' @param frame,frame_par
#'   Elements of the \link{list} output of \code{\link{make_frames}}.
#' @inheritParams egf
#'
#' @return
#' A \link{list} with elements \code{data}, \code{parameters}, \code{map},
#' \code{random}, \code{profile}, \code{DLL}, and \code{silent}.
#'
#' @seealso \code{\link{make_tmb_data}}, \code{\link{make_tmb_parameters}}
#' @keywords internal
make_tmb_args <- function(frame, frame_par, model, control, do_fit, init) {
  tmb_data <- make_tmb_data(
    frame = frame,
    frame_par = frame_par,
    model = model,
    control = control
  )
  tmb_parameters <- make_tmb_parameters(
    tmb_data = tmb_data,
    frame = frame,
    frame_par = frame_par,
    model = model,
    do_fit = do_fit,
    init = init
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
    profile = if (control$profile) "beta",
    DLL = "epigrowthfit",
    silent = (control$trace == 0L)
  )
}

#' Construct data objects for C++ template
#'
#' Gathers in a \link{list} data objects to be passed to
#' the package's C++ template via \pkg{TMB}'s \code{DATA_*} macros.
#' See also \code{\link[TMB]{MakeADFun}}.
#'
#' @inheritParams make_tmb_args
#'
#' @return
#' [Below,
#' \code{n = sum(!is.na(frame$window))}
#' is the total number of time points belonging a fitting window,
#' \code{N = nlevels(frame$window)}
#' is the number of fitting windows, and
#' \code{p = length(frame_par)}
#' is the number of nonlinear and dispersion model parameters.]
#'
#' A \link{list} inheriting from \link{class} \code{"tmb_data"},
#' with elements:
#' \item{t}{
#'   A \link[=double]{numeric} vector of length \code{n} giving time
#'   as a number of days since the earliest time point in the current
#'   fitting window.
#' }
#' \item{t_seg_len}{
#'   An \link{integer} vector of length \code{N} specifying the
#'   length of each fitting window as a number of time points.
#' }
#' \item{x}{
#'   A \link[=double]{numeric} vector of length \code{n-N} giving
#'   incidence in each fitting window. \code{x[i]} in window \code{k}
#'   is the number of cases observed from time \code{t[k+i-1]}
#'   to time \code{t[k+i]}.
#' }
#' \item{day_of_week_on_day0}{
#'   If \code{model$day_of_week > 0}, then an \link{integer} vector
#'   of length \code{N} indicating the first day of week in each
#'   fitting window, with value \code{i} in \code{0:6} mapping to
#'   the day of week \code{i} days after the reference day specified
#'   by \code{model$day_of_week}. Otherwise, an integer vector of
#'   the same length filled with \code{-1}.
#' }
#' \item{Yo}{
#'   The \link[=model.offset]{offset} matrix in dense format,
#'   with \code{N} rows and \code{p} columns.
#' }
#' \item{X}{
#'   The fixed effects \link[=model.matrix]{design} matrix
#'   in \link[Matrix:sparseMatrix]{sparse} or dense format,
#'   depending on \code{control$sparse_X}, with \code{N} rows.
#' }
#' \item{Z}{
#'   The random effects design matrix in
#'   \link[Matrix:sparseMatrix]{sparse} format, with \code{N} rows.
#'   If there are no random effects, then \code{Z} is an empty sparse matrix.
#' }
#' \item{Xs, Xd}{
#'   If \code{control$sparse_X = TRUE}, then \code{Xs = X} and \code{Xd}
#'   is an empty dense matrix. Otherwise, \code{Xd = X} and \code{Xs} is
#'   an empty sparse matrix.
#' }
#' \item{X_info, Z_info}{
#'   \link[=data.frame]{Data frame}s with \code{\link{ncol}(X)} and
#'   \code{\link{ncol}(Z)} rows, respectively. Row \code{j} describes
#'   the coefficient associated with column \code{j} of \code{X} or
#'   \code{Z}. \code{Z_info} has additional \link{factor}s \code{cor}
#'   and \code{vec} \link{split}ting coefficients by relation to a
#'   common covariance matrix and random vector, respectively.
#' }
#' \item{beta_seg_len, b_seg_len}{
#'   \link[=integer]{Integer} vectors of length \code{p} counting the
#'   columns of \code{X} and \code{Z}, respectively, that relate to a
#'   common nonlinear or dispersion model parameter.
#' }
#' \item{beta_seg_index, b_seg_index}{
#'   \link[=integer]{Integer} vectors of length \code{\link{ncol}(X)}
#'   and \code{\link{ncol}(Z)}, respectively, with values in \code{0:(p-1)}.
#'   These split the columns of \code{X} and \code{Z} by relation to
#'   a common nonlinear or dispersion model parameter.
#' }
#' \item{block_rows, block_cols}{
#'   \link[=integer]{Integer} vectors together giving the dimensions
#'   of each block of the random effects matrix.
#' }
#' \item{curve_flag, excess_flag, family_flag, day_of_week_flag, trace_flag, sparse_X_flag}{
#'   \link[=integer]{Integer} flags referencing so-named elements
#'   of \code{model} and \code{control}.
#' }
#' \item{`predict_flag`}{
#'   An \link{integer} flag set equal to 0 so that prediction code is not run.
#' }
#' Additional \link{integer} elements of the form \code{j_link_parameter}
#' (e.g., \code{j_log_r}) give the 0-index of nonlinear and dispersion
#' model parameter names (e.g., \code{"log(r)"}) in \code{names(frame_par)}.
#' The value \code{-1} is used for parameters not belonging to the model
#' being fit.
#'
#' @keywords internal
#' @importFrom stats formula terms model.matrix model.offset
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix sparse.model.matrix KhatriRao
make_tmb_data <- function(frame, frame_par, model, control, do_fit, init) {
  ## Nonlinear and dispersion model parameter names
  pn <- names(frame_par)
  p <- length(frame_par)

  ## Time series model frame excluding rows not associated
  ## with a fitting window
  origin <- attr(frame, "origin")
  frame <- frame[!is.na(frame$window), , drop = FALSE]

  ## Fitting window lengths as numbers of time points
  t_seg_len <- tabulate(frame$window)
  N <- length(t_seg_len)

  ## Time as number of days since earliest time point
  firsts <- which(!duplicated(frame$window))
  t <- frame$time - rep.int(frame$time[firsts], t_seg_len)

  ## Incidence without unused first elements
  x <- frame$x[-firsts]

  ## Day of week during 1-day interval starting at earliest time point
  ## NB: When `day_of_week > 0L`, time points are integer numbers of
  ##     days since midnight on a reference date, so this 1-day interval
  ##     indeed occurs on a unique day of week
  if (model$day_of_week > 0L) {
    ## Date during that interval
    ## NB: Reason for adding 1 is that the earliest time points are
    ##     23:59:59 on Dates `origin + time[firsts]`, so the 1-day
    ##     intervals that we seek occur on the following dates
    Date_on_day0 <- origin + frame$time[firsts] + 1
    ## Day of week during that interval, coded as an integer `i` in `0:6`
    ## NB: `i` maps to the day of week `i` days after the reference day,
    ##     which is the day `day_of_week` days after an arbitrary Saturday,
    ##     in this case `.Date(2)`
    day_of_week_on_day0 <- as.integer(julian(Date_on_day0, origin = .Date(2 + model$day_of_week)) %% 7)
  } else {
    day_of_week_on_day0 <- rep_len(-1L, N)
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
  Xl <- Map(make_X, x = fixed, frame = frame_par, sparse = control$sparse_X)
  X <- do.call(cbind, Xl)

  ## Random effects design matrices
  Zl <- Map(make_Z, x = do.call(c, unname(random)), frame = rep.int(frame_par, lengths(random)))
  Z <- do.call(cbind, Zl)

  ## Nonlinear or dispersion model parameter, formula term,
  ## group (g in (terms | g)), and group level associated
  ## with each column of each matrix
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
  ## convenient arrangement in random effects blocks,
  ## given column-major order
  o <- do.call(order, unname(Z_info[c("cor", "vec", "par")]))
  Z <- Z[, o, drop = FALSE]
  Z_info <- Z_info[o, , drop = FALSE]

  ## Number of coefficients related to each nonlinear and
  ## dispersion model parameter
  beta_seg_len <- c(table(X_info$par))
  b_seg_len <- c(table(Z_info$par))

  ## Coefficients factored by relation to a nonlinear or
  ## dispersion model parameter
  beta_seg_index <- as.integer(X_info$par) - 1L
  b_seg_index <- as.integer(Z_info$par) - 1L

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

  l1 <- list(
    t = t,
    t_seg_len = t_seg_len,
    x = x,
    day_of_week_on_day0 = day_of_week_on_day0,
    Yo = Yo,
    X = X,
    Z = if (is.null(Z)) empty_sparse_matrix else Z,
    Xs = if (control$sparse_X) X else empty_sparse_matrix,
    Xd = if (control$sparse_X) empty_dense_matrix else X,
    X_info = X_info,
    Z_info = Z_info,
    beta_seg_len = beta_seg_len,
    b_seg_len = b_seg_len,
    beta_seg_index = beta_seg_index,
    b_seg_index = b_seg_index,
    block_rows = block_rows,
    block_cols = block_cols,
    curve_flag = get_flag("curve", model$curve),
    excess_flag = as.integer(model$excess),
    family_flag = get_flag("family", model$family),
    day_of_week_flag = as.integer(model$day_of_week > 0L),
    trace_flag = control$trace,
    sparse_X_flag = as.integer(control$sparse_X),
    predict_flag = 0L
  )
  pn0 <- get_par_names(NULL, link = TRUE)
  l2 <- as.list(match(pn0, pn, 0L) - 1L)
  names(l2) <- sub("^(log|logit)\\((.*)\\)$", "j_\\1_\\2", pn0)
  out <- c(l1, l2)
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
#' When \code{init = NULL}, naive estimates of nonlinear model parameters
#' are obtained for each fitting window as follows:
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
#'   \code{exp(log_c0)}, where \code{log_c0} is the intercept
#'   of a linear model fit to \code{\link{log1p}(\link{cumsum}(x))}.
#' }
#' \item{\code{tinfl}}{
#'   \code{max(t)}, assuming that the fitting window ends near
#'   the time of a peak in interval incidence.
#' }
#' \item{\code{K}}{
#'   \code{2*sum(x)}, assuming that the fitting window ends near
#'   the time of a peak in interval incidence _and_ that interval
#'   incidence is roughly symmetric about the peak.
#' }
#' \item{\code{p}}{0.8}
#' \item{\code{a, b, nbdisp, w[123456]}}{1}
#' }
#' The naive estimates are log- or logit-transformed (all but
#' \code{p} use a log link), and, for each nonlinear or dispersion
#' model parameter, the initial value of the \code{"(Intercept)"}
#' coefficient in its fixed effects model (if one exists) is taken
#' to be the mean of the link scale estimates across fitting windows.
#' The remaining elements of \code{beta} take initial value 0.
#' All elements of \code{b} and \code{theta} take initial value 0,
#' so that, initially, all random vectors follow a standard normal
#' distribution.
#'
#' @return
#' A \link{list} with elements:
#' \item{beta}{
#'   A \link[=double]{numeric} vector of length
#'   \code{sum(tmb_data$beta_seg_len)} listing
#'   initial values for the fixed effects coefficients.
#' }
#' \item{b}{
#'   If there are random effects, then a \link[=double]{numeric}
#'   vector of length \code{sum(tmb_data$b_seg_len)} listing
#'   initial values for the (unit variance scale) random effects
#'   coefficients. If there are no random effects,
#'   then \code{b = \link{NA_real_}}.
#' }
#' \item{theta}{
#'   If there are random effects, then a numeric vector of length
#'   \code{sum(f(tmb_data$block_rows))}, where \code{f(n) = n*(n+1)/2},
#'   specifying initial covariance matrices for all random vectors.
#'   Let \code{u} be a segment of \code{theta} of length \code{f(n)},
#'   corresponding to a random vector \code{x} of length \code{n}.
#'   Then \code{u[1:n]} lists log standard deviations of the elements
#'   of \code{x}, and \code{u[-(1:n)]} lists (in row-major order)
#'   the subdiagonal elements of the unit lower triangular matrix
#'   \code{L} in the Cholesky factorization
#'   \code{S = A \link[=matmult]{\%*\%} \link{t}(A)},
#'   \code{A = D^-0.5 \link[=matmult]{\%*\%} L},
#'   \code{D = \link{diag}(diag(L \link[=matmult]{\%*\%} \link{t}(L)))}
#'   of the correlation matrix \code{S} of \code{x}.
#'   If there are no random effects, then \code{theta = \link{NA_real_}}.
#' }
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis terms
make_tmb_parameters <- function(tmb_data, frame, frame_par, model, do_fit, init) {
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
      window <- frame$window[!is.na(frame$window)]
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
      alpha <- switch(model$curve,
                      subexponential = r * c0^(1 - p),
                      gompertz = r / log(K / c0),
                      NA_real_
      )
      Y_init <- data.frame(r, alpha, c0, tinfl, K, p)
      Y_init[c("a", "b", "nbdisp", paste0("w", 1:6))] <- 1

      ## Link transform
      Y_init[] <- Map(function(x, s) match_link(s)(x),
                      x = Y_init,
                      s = string_get_link(names(Y_init))
      )
      names(Y_init) <- string_add_link(names(Y_init))

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
  if (is.null(init) && !do_fit) {
    attr(init_split, "Y_init") <- as.matrix(Y_init[names(frame_par)])
  }
  init_split
}

#' Check for random effects
#'
#' Determine whether an object specifies a random effects model.
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
#' Defines a wrapper function on top of \code{\link[TMB]{MakeADFun}}-generated
#' functions \code{fn} and \code{gr}, so that function and gradient evaluations
#' can retry inner optimization using fallback methods in event of failure using
#' the default method (usually \code{\link[TMB]{newton}}).
#'
#' @param f
#'   A \link{function} to be patched, either \code{fn} or \code{gr}
#'   from a \code{\link[TMB]{MakeADFun}}-generated \link{list} object.
#' @param order
#'   An integer, either 0 or 1,
#'   identifying \code{f} as \code{fn} (0) or \code{gr} (1).
#' @param inner_optimizer
#'   An \code{"\link{egf_inner_optimizer}"} object listing
#'   inner optimization methods to be tried in turn.
#'
#' @return
#' A patched version of \link{function} \code{f}.
#'
#' @keywords internal
patch_fn_gr <- function(f, order, inner_optimizer) {
  stop_if_not(
    is.numeric(order),
    length(order) == 1L,
    order %in% 0:1,
    m = "`order` must be 0 or 1."
  )
  f_name <- c("fn", "gr")[1L + order]
  e <- environment(f)
  function(x = e$last.par[-e$random], ...) {
    oim <- e$inner.method
    oic <- e$inner.control
    on.exit({
      e$inner.method <- oim
      e$inner.control <- oic
    })
    n <- length(x)^order
    ok <- function(v) is.numeric(v) && length(v) == n && all(is.finite(v))
    for (io in inner_optimizer) {
      e$inner.method <- io$method
      e$inner.control <- io$control
      v <- f(x, ...)
      if (ok(v)) {
        return(v)
      }
    }
    warning(sprintf("Unable to evaluate `%s(x)`, returning NaN.", f_name))
    NaN
    # warning(sprintf("Unable to evaluate `%s(x)`, returning `rep_len(NaN, %d)`.", f_name, n))
    # rep_len(NaN, n)
  }
}
