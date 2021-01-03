#' Get flags
#'
#' Returns the underlying integer value of an enumerator in the
#' C++ template.
#'
#' @param type,enum
#'   Character strings giving the name of an enumerated type
#'   and a name of an enumerator of that type, respectively.
#'
#' @details
#' This function is kept synchronized with the C++ template
#' using script `utils/update_flags.R`.
#'
#' @return
#' An integer.
#'
#' @examples
#' epigrowthfit:::get_flag("curve", "logistic")
#' epigrowthfit:::get_flag("distr", "nbinom")
#'
#' @keywords internal
get_flag <- function(type, enum) {
  curve_names <- c("exponential", "subexponential", "gompertz", "logistic", "richards") # GREP_FLAG
  distr_names <- c("pois", "nbinom") # GREP_FLAG
  switch(type,
    curve = match(enum, curve_names) - 1L,
    distr = match(enum, distr_names) - 1L,
    NA_integer_
  )
}

#' Get parameter names
#'
#' Returns the names used internally for parameters of the specified
#' incidence model.
#'
#' @param curve,distr,excess
#'   Character strings specifying an incidence model. See [egf()].
#'   Alternatively, `curve` can be an `"egf"` object returned by
#'   [egf()] or a `"tmb_data"` object returned by [make_tmb_data()],
#'   in which case `distr` and `excess` are ignored.
#' @param link
#'   A logical scalar. If `TRUE`, then a prefix indicating a link
#'   function (either `"log_"` or `"logit_"`) is prepended to each
#'   parameter name.
#' @param s
#'   A character vector listing parameter names with or without
#'   prefixes.
#'
#' @details
#' `add_link_string()` (internal) adds prefixes to parameter names.
#' `remove_link_string()` (internal) is its inverse.
#'
#' @return
#' `get_par_names(curve, distr, excess, link = FALSE)` returns the
#' subset of
#' `c("r", "alpha", "c0", "tinfl", "K", "p", "a", "b", "nbdisp")`
#' relevant to `curve`, `distr`, and `excess`.
#'
#' @examples
#' pn  <- get_par_names(curve = "exponential", distr = "pois",
#'                      excess = FALSE, link = FALSE)
#' lpn <- get_par_names(curve = "exponential", distr = "pois",
#'                      excess = FALSE, link = TRUE)
#' identical(pn,  epigrowthfit:::remove_link_string(lpn))
#' identical(lpn, epigrowthfit:::add_link_string(pn))
#'
#' @export
get_par_names <- function(curve = NULL, distr = NULL, excess = NULL,
                          link = TRUE) {
  UseMethod("get_par_names", curve)
}

#' @export
get_par_names.default <- function(curve = NULL, distr = NULL, excess = NULL,
                                  link = TRUE) {
  if (is.null(curve) && is.null(distr) && is.null(excess)) {
    pn <- c("r", "alpha", "c0", "tinfl", "K", "p", "a", "b", "nbdisp")
  } else {
    pn <- character(0)
    if (!is.null(curve)) {
      a <- switch(curve,
        exponential    = c("r", "c0"),
        subexponential = c("alpha", "c0", "p"),
        gompertz       = c("alpha", "c0", "K"),
        logistic       = c("r", "tinfl", "K"),
        richards       = c("r", "tinfl", "K", "a"),
      )
      pn <- c(pn, a)
    }
    if (!is.null(distr) && distr == "nbinom") {
      pn <- c(pn, "nbdisp")
    }
    if (!is.null(excess) && excess) {
      pn <- c(pn, "b")
    }
  }
  if (link) {
    add_link_string(pn)
  } else {
    pn
  }
}

#' @export
get_par_names.tmb_data <- function(curve, distr = NULL, excess = NULL,
                                   link = TRUE) {
  if (link) {
    colnames(curve$rid)
  } else {
    remove_link_string(colnames(curve$rid))
  }
}

#' @export
get_par_names.egf <- function(curve, distr = NULL, excess = NULL,
                              link = TRUE) {
  get_par_names(curve$tmb_args$data)
}

#' @rdname get_par_names
#' @keywords internal
add_link_string <- function(s) {
  if (is.null(s)) {
    return(NULL)
  }
  ok <- s %in% get_par_names(link = FALSE)
  if (!any(ok)) {
    return(s)
  }
  fmt <- ifelse(s[ok] == "p", "logit_%s", "log_%s")
  s[ok] <- sprintf(fmt, s[ok])
  s
}

#' @rdname get_par_names
#' @keywords internal
remove_link_string <- function(s) {
  if (is.null(s)) {
    return(NULL)
  }
  ok <- s %in% get_par_names(link = TRUE)
  if (!any(ok)) {
    return(s)
  }
  s[ok] <- sub("^(logit|log)_", "", s[ok])
  s
}

#' Validate model formulae
#'
#' Check that formulae passed to [egf()] are specified correctly
#' by deparsing and comparing against an accepted pattern.
#'
#' @inheritParams egf
#'
#' @return
#' `check_formula()` returns `formula` as is.
#'
#' `check_fixed()` and `check_random()` return lists of `p` formulae
#' with names `get_par_names(curve, distr, excess, link = TRUE)`.
#' If `fixed` and `random` are already subsets of such a list, then
#' those subsets are preserved, and the set difference is made up of
#' `~1` and `NULL`, respectively. If `fixed` and `random` are formulae
#' (rather than lists of formulae) or if `random` is `NULL`, then they
#' are recycled to length `p`.
#'
#' @details
#' In R 4.0.3, formulae are deparsed with one or zero spaces on each
#' side of the tilde operator `~` (one if there is a left-hand side,
#' zero if not), zero spaces around `( ) : /`, and one space around
#' `+ - * |`. `check_formula()` and friends unsafely assume that this
#' behavior is maintained across R versions, so unexpected breakage
#' is possible.
#'
#' FIXME: Refactor using [stats::terms()] instead of deparse-regex
#' to avoid breakage. Ideally, allow for arithmetic, nonsyntactic
#' variable names, and spurious parentheses, which are currently not
#' tolerated (see Examples).
#'
#' @examples
#' curve  <- "exponential"
#' distr  <- "pois"
#' excess <- FALSE
#' get_par_names(curve, distr, excess, link = TRUE)
#' ## [1] "log_r" "log_c0"
#'
#' cfm <- epigrowthfit:::check_formula
#' cfx <- function(fixed) {
#'   epigrowthfit:::check_fixed(fixed, curve, distr, excess)
#' }
#' crd <- function(random) {
#'   epigrowthfit:::check_random(random, curve, distr, excess)
#' }
#'
#' ## Time series
#' cfm(y ~ x)
#'
#' ## Fixed effects
#' cfx(list(log_r = ~z, log_c0 = ~z))
#' cfx(list(log_r = ~z0, log_c0 = ~z1:z2))
#' cfx(~z) # same as first line
#' cfx(list(r = ~z)) # `c0` defaults to `~1` (global mean)
#'
#' ## Random effects
#' crd(list(log_r = ~(1|z), log_c0 = ~(1|z)))
#' crd(list(log_r = ~(1|z0), log_c0 = ~(1|z1:z2)))
#' crd(~(1|z)) # same as first line
#' crd(list(r = ~(1|z))) # `c0` defaults to NULL (no RE)
#' crd(list(r = ~(1|z0), c0 = ~(1|z1/z2))) # RE component allowed >1 term,
#'                                         # `z1/z2` expands to 2 terms
#'
#' ## Not run
#' \dontrun{
#' cfm(I(y1 + y2) ~ x)  # error: arithmetic
#' cfm("0y" ~ x)        # error: nonsyntactic variable name
#' cfm((y) ~ x)         # error: spurious parentheses
#' cfx(list(r = ~z, c = ~z))       # error: misnamed parameter
#' cfx(list(r = ~z0, c0 = ~z1*z2)) # error: FE component restricted to 1 term,
#'                                 #        `z1*z2` expands to 3 terms
#' cfx(NULL)                       # error: for null FE model use `~1`
#' crd(~1)                         # error: RE terms need >1 group
#' }
#'
#' @name check_formula
#' @keywords internal
check_formula <- function(formula) {
  stop_if_not(
    inherits(formula, "formula"),
    length(formula) == 3L,
    grepl("^[[:alpha:].]{1}[[:alnum:]._]* ~ [[:alpha:].]{1}[[:alnum:]._]*$", deparse(formula)),
    m = "`formula` must be a formula of the form `y ~ x`."
  )
  formula
}

#' @rdname check_formula
#' @keywords internal
check_fixed <- function(fixed, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (inherits(fixed, "formula")) {
    fixed <- rep(list(fixed), p)
    names(fixed) <- pn
  }
  stop_if_not(
    inherits(fixed, "list"),
    length(fixed) >= 1L,
    vapply(fixed, inherits, logical(1L), "formula"),
    m = "`fixed` must be a formula or a named list of formulae."
  )
  ## Tolerate, e.g., "r" instead of "log_r" but not both
  names(fixed) <- add_link_string(remove_link_string(names(fixed)))
  stop_if_not(
    !is.null(names(fixed)),
    names(fixed) %in% pn,
    !any(duplicated(names(fixed))),
    m = paste0(
      "If `fixed` is a list, then `names(fixed)` must be a subset\n",
      "of `get_par_names(curve, distr, excess)`."
    )
  )
  rhs <- vapply(fixed, function(x) deparse(x[[2L]]), character(1L))
  stop_if_not(
    lengths(fixed) == 2L,
    grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", rhs),
    m = "`fixed` formulae must be `~1` or have the form `~f1:...:fn`."
  )

  ## Fill out and order the list
  fixed[setdiff(pn, names(fixed))] <- list(~1)
  fixed[pn]
}

#' @rdname check_formula
#' @keywords internal
check_random <- function(random, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (is.null(random) || inherits(random, "formula")) {
    random <- rep(list(random), p)
    names(random) <- pn
  }
  stop_if_not(
    inherits(random, "list"),
    length(random) >= 1L,
    vapply(random, inherits, logical(1L), c("NULL", "formula")),
    m = "`random` must be NULL, a formula, or a named list of formulae."
  )
  ## Tolerate, e.g., "r" instead of "log_r" but not both
  names(random) <- add_link_string(remove_link_string(names(random)))
  stop_if_not(
    !is.null(names(random)),
    names(random) %in% pn,
    !any(duplicated(names(random))),
    m = paste0(
      "If `random` is a list, then `names(random)` must be a subset\n",
      "of `get_par_names(curve, distr, excess)`."
    )
  )
  is_formula <- vapply(random, inherits, logical(1L), "formula")
  if (any(is_formula)) {
    rhs <- vapply(random[is_formula], function(x) deparse(x[[2L]]), character(1L))
    stop_if_not(
      lengths(random[is_formula]) == 2L,
      grepl("^(( (\\+|-) )?\\(1 \\| [[:alnum:]._]+((:[[:alnum:]._]+)*|(/[[:alnum:]._]+)*|(( \\* )[[:alnum:]._]+)*)\\))+$", rhs),
      m = paste0(
        "`random` formulae must have the form `~rhs`,\n",
        "with `rhs` a sum of one or more terms of the form\n",
        "`(1 | r1:...:rk)`, `(1 | r1/.../rm)`, or `(1 | r1 * ... * rn)`."
      )
    )
  }

  ## Fill out and order the list
  random[setdiff(pn, names(random))] <- list(NULL)
  random[pn]
}

#' Construct a model frame
#'
#' Constructs a model frame to be used by [egf()], after performing
#' myriad checks on the arguments.
#'
#' @inheritParams egf
#'
#' @details
#' Rows of `data` are permuted by `order(index)` so that rows belonging
#' to the same fitting window are contiguous in the model frame.
#' The permuted rows indexed by `!is.na(index[order(index)])` form the
#' model frame. The remaining rows are preserved in attribute `extra`.
#'
#' @return
#' A data frame containing the variables named in `formula`, `fixed`,
#' and `random`, with attributes:
#' \item{`extra`}{
#'   A data frame preserving data not belonging to a fitting window
#'   (see Details).
#'   Plot methods may use all of the data in `rbind(frame, extra)`.
#' }
#' \item{`index`}{
#'   A factor of length `nrow(frame)` such that `split(frame, index)`
#'   splits the frame by fitting window. Not necessarily identical to
#'   the argument of the same name (see Details).
#' }
#' \item{`date_name`, `cases_name`}{
#'   Characters strings naming the date and incidence variables in the
#'   data frame.
#' }
#' \item{`fixed_term_labels`, `random_term_labels`}{
#'   Named lists of character vectors naming mixed effects model terms
#'   for each parameter.
#' }
#'
#' @examples
#' r <- log(2) / 10
#' mu <- diff(10 * exp(r * 0:99))
#'
#' date <- rep(as.Date(0:99, origin = "2020-01-01"), 4L)
#' cases <- c(replicate(4L, c(NA, rpois(length(mu), mu))))
#' continent <- rep(c("asia", "europe"), rep(200L, 2L))
#' country <- rep(c("china", "japan", "france", "germany"), rep(100L, 4L))
#' wave <- rep(rep(c(NA, 1, NA, 2, NA), rep(20L, 5L)), 4L)
#'
#' data <- data.frame(date, cases, continent, country, wave)
#' data[3:5] <- lapply(data[3:5], factor)
#'
#' x <- c(NA, 1, NA, 2, NA,
#'        NA, 3, NA, 4, NA,
#'        NA, 5, NA, 6, NA,
#'        NA, 7, NA, 8, NA)
#' index <- factor(rep(x, rep(20L, 20L)))
#'
#' formula <- cases ~ date
#' fixed <- list(
#'   r  = ~1,
#'   c0 = ~wave
#' )
#' random <- list(
#'   r  = ~(1 | continent/country) + (1 | wave),
#'   c0 = ~(1 | continent/country)
#' )
#'
#' curve <- "exponential"
#' distr <- "pois"
#' excess <- FALSE
#' na_action <- "pass"
#' date_format <- "%Y-%m-%d" # unused in this case because
#'                           # inherits(date, "Date") = TRUE
#'
#' frame <- epigrowthfit:::make_frame(formula, fixed, random,
#'                                    data, index,
#'                                    curve, distr, excess,
#'                                    na_action, date_format)
#'
#' @keywords internal
make_frame <- function(formula, fixed, random, data, index,
                       curve, distr, excess, na_action, date_format) {
  ## Check that model formulae are specified correctly
  formula <- check_formula(formula)
  fixed   <- check_fixed(fixed, curve, distr, excess)
  random  <- check_random(random, curve, distr, excess)

  ## Check that formula variables can be found
  stop_if_not(
    inherits(data, c("data.frame", "list", "environment")),
    m = "`data` must be a data frame, list, or environment."
  )
  if (is.environment(data)) {
    data <- as.list(data)
  }
  dn <- all.vars(formula[[3L]]) # date
  cn <- all.vars(formula[[2L]]) # cases
  fn <- unique(unlist(lapply(c(fixed, random), all.vars))) # factors
  an <- unique(c(dn, cn, fn))
  found <- an %in% names(data)
  if (any(!found)) {
    stop(
      "Variables not found in `data`:\n",
      paste(an[!found], collapse = ", ")
    )
  }

  ## Check that time series is specified correctly
  if (is.character(data[[dn]])) {
    data[[dn]] <- try(as.Date(data[[dn]], tryFormats = date_format), silent = TRUE)
  }
  stop_if_not(
    inherits(data[[dn]], "Date"),
    m = sprintf("`%s` must be of class \"Date\" or so coercible with\n`as.Date(%s, tryFormats = date_format)`.", dn, dn)
  )
  stop_if_not(
    length(data[[dn]]) > 0L,
    m = sprintf("`%s` must have nonzero length.", dn)
  )
  stop_if_not(
    !anyNA(data[[dn]]),
    m = sprintf("`%s` must not have missing values.", dn)
  )
  n <- length(data[[dn]])
  stop_if_not(
    is.numeric(data[[cn]]),
    length(data[[cn]]) == n,
    m = sprintf("`%s` must be numeric and have length `length(%s)`.", cn, dn)
  )
  stop_if_not(
    all(data[[cn]] >= 0, na.rm = TRUE),
    !is.infinite(data[[cn]]),
    all(data[[cn]] %% 1 == 0, na.rm = TRUE), # FIXME: isTRUE(all.equal())?
    m = sprintf("Elements of `%s` must be non-negative integers or NA.", cn)
  )
  data[[cn]] <- as.integer(data[[cn]])

  ## Check that grouping variables are specified correctly
  if (length(fn) > 0L) {
    stop_if_not(
      vapply(data[fn], is.factor, logical(1L)),
      lengths(data[fn]) == n,
      m = sprintf("Grouping variables must be factors of length `length(%s)`.", dn)
    )
    stop_if_not(
      !vapply(data[fn], anyNA, logical(1L)),
      m = "Grouping variables must not have missing values."
    )
    data[fn] <- lapply(data[fn], droplevels, exclude = NA)
  }

  ## Check that fitting windows are specified correctly
  frame <- data.frame(data[an])
  if (is.null(index)) {
    index <- factor(integer(n))
  } else {
    stop_if_not(
      is.factor(index),
      length(index) == n,
      m = sprintf("`index` must be a factor of length `length(%s)`.", dn)
    )
    index <- droplevels(index, exclude = NA)
    stop_if_not(
      nlevels(index) > 0L,
      m = "`index` must have at least one nonempty level."
    )
  }
  frame_split <- split(frame, index)
  date_split <- lapply(frame_split, "[[", dn)
  stop_if_not(
    lengths(date_split) >= 7L,
    m = sprintf("`%s` must have length 7 or greater in each level of `index`.", dn)
  )
  stop_if_not(
    vapply(date_split, function(x) all(diff(x) > 0), logical(1L)),
    m = sprintf("`%s` must be increasing in each level of `index`.", dn)
  )
  cases_split <- lapply(frame_split, "[[", cn)
  if (na_action == "fail") {
    stop_if_not(
      vapply(cases_split, function(x) !anyNA(x[-1L]), logical(1L)),
      m = sprintf("There are fitting windows with missing values in `%s`.", cn)
    )
  } else {
    stop_if_not(
      vapply(cases_split, function(x) sum(!is.na(x)) >= 7L, logical(1L)),
      m = sprintf("There are fitting windows with insufficient data\n(fewer than 7 observations) in `%s`.", cn)
    )
  }
  if (length(fn) > 0L) {
    factors_split <- lapply(frame_split, function(x) droplevels(x[fn]))
    stop_if_not(
      vapply(factors_split, function(x) all(vapply(x, nlevels, integer(1L)) == 1L), logical(1L)),
      m = paste0(
        "Grouping variables must have exactly one level\n",
        "in each level of `index`."
      )
    )
    stop_if_not(
      vapply(frame[fn], nlevels, integer(1L)) > 1L,
      m = paste0(
        "Grouping variables must not have the same level\n",
        "in every level of `index`."
      )
    )
  }

  ## Order by index
  ord <- order(index)
  index <- index[ord]
  frame <- frame[ord, ]
  row.names(frame) <- NULL

  ## Subset rows belonging to a fitting window
  ## and preserve the difference as an attribute
  ## for use by plot methods
  keep <- !is.na(index)
  structure(droplevels(frame[keep, , drop = FALSE]),
    extra = droplevels(frame[!keep, , drop = FALSE]),
    index = index[keep],
    date_name = dn,
    cases_name = cn,
    fixed_term_labels = lapply(fixed, get_term_labels),
    random_term_labels = lapply(random, get_term_labels)
  )
}

#' Get term labels from a formula
#'
#' A wrapper for `labels(terms())` with special handling of certain
#' formulae.
#'
#' @param f A formula.
#'
#' @return
#' If `f` has no variable terms, as in `f = ~1`, then `"(1)"`.
#' Otherwise, the result of `labels(terms(f))` *after* replacing
#' terms of `f` of the form `(1 | rhs_of_bar)` with `rhs_of_bar`.
#'
#' @examples
#' epigrowthfit:::get_term_labels(~w + x + y:z)
#' epigrowthfit:::get_term_labels(~1)
#' epigrowthfit:::get_term_labels(~w*x + y/z) # crosses, nests expanded
#' epigrowthfit:::get_term_labels(~w*x + (1 | y/z)) # bars ignored
#'
#' @keywords internal
#' @importFrom stats terms reformulate
get_term_labels <- function(f) {
  if (is.null(f)) {
    return(NULL)
  }
  tl <- labels(terms(f))
  if (length(tl) == 0L) {
    return("(1)")
  }
  has_bar <- grepl("^1 \\| ", tl)
  if (any(has_bar)) {
    tl[has_bar] <- sub("^1 \\| ", "", tl[has_bar])
    tl <- labels(terms(reformulate(tl)))
  }
  tl
}

#' Evaluate a term label in a data frame
#'
#' Constructs the interaction specified by a term label
#' (a character string of the form `"f1:...:fn"`) from
#' factors in a data frame.
#'
#' @param term_label
#'   A character string giving a colon-separated list of names
#'   of factors in `frame`.
#' @param frame
#'   A data frame containing factors.
#'
#' @return
#' A factor of length `nrow(frame)` giving the interaction
#' specified by `term_label`. Unused levels are dropped.
#'
#' In the special case `term_label = "(1)"`, the result is
#' `factor(rep("(1)", nrow(frame)))`.
#'
#' @examples
#' f <- function() factor(sample(5L, 20L, replace = TRUE))
#' d <- data.frame(x = f(), y = f())
#' epigrowthfit:::get_factor("x:y", frame = d)
#' epigrowthfit:::get_factor("(1)", frame = d)
#'
#' @keywords internal
get_factor <- function(term_label, frame) {
  if (is.null(term_label)) {
    return(NULL)
  }
  if (term_label == "(1)") {
    return(factor(rep("(1)", nrow(frame))))
  }
  s <- unique(unlist(strsplit(term_label, ":")))
  interaction(frame[s], drop = TRUE, sep = ":")
}

#' Construct a design matrix from a factor
#'
#' A wrapper for
#' [stats::model.matrix()] and [Matrix::sparse.model.matrix()]
#' operating on factors instead of formulae and data frames.
#'
#' @param x
#'   A factor.
#' @param sparse
#'   A logical scalar. If TRUE, then matrix is returned in sparse
#'   format.
#' @param intercept
#'   A logical scalar. If TRUE, then zero-sum contrasts are used.
#'   If FALSE, then no contrasts are used. See [stats::contr.sum()].
#'
#' @details
#' Factors `x` with `nlevels(x) < 2` are tolerated. This is not true
#' of [stats::model.matrix()] or [Matrix::sparse.model.matrix()].
#'
#' @return
#' A matrix with `length(x)-d` rows and `nlevels(x)` columns,
#' where `d` is `sum(is.na(x))` if `NA` is not among `levels(x)`
#' and 0 otherwise.
#'
#' @examples
#' x <- factor(sample(5L, 20L, replace = TRUE))
#' epigrowthfit:::factor_to_matrix(x, sparse = TRUE, intercept = TRUE)
#' epigrowthfit:::factor_to_matrix(x, sparse = TRUE, intercept = FALSE)
#' epigrowthfit:::factor_to_matrix(x, sparse = FALSE, intercept = TRUE)
#' epigrowthfit:::factor_to_matrix(x, sparse = FALSE, intercept = FALSE)
#'
#' @keywords internal
#' @importFrom Matrix sparseMatrix sparse.model.matrix
#' @importFrom stats model.matrix
factor_to_matrix <- function(x, sparse, intercept) {
  if (is.null(x)) {
    return(NULL)
  }
  d <- if (!anyNA(levels(x))) sum(is.na(x)) else 0L
  m <- length(x) - d
  n <- nlevels(x)
  if (m == 0L || n < 2L) {
    if (sparse) {
      X <- sparseMatrix(seq_len(m * n), rep(1L, m * n), x = rep(1L, m * n), dims = c(m, n))
    } else {
      X <- matrix(1L, nrow = m, ncol = n)
    }
  } else {
    fn <- if (sparse) "sparse.model.matrix" else "model.matrix"
    f <- get(fn)
    X <- f(
      object = if (intercept) ~x else ~-1 + x,
      data = list(x = x),
      contrasts.arg = list(x = "contr.sum")
    )
    dimnames(X) <- list(NULL, NULL)
  }
  X
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
#' \item{`index`}{
#'   A copy of `attr(frame, "index")`. See [make_frame()].
#' }
#' \item{`t`}{
#'   An integer vector of length `nrow(frame)` giving time as a number
#'   of days since the earliest date in the fitting window.
#' }
#' \item{`x`}{
#'   An integer vector of length `nrow(frame)`. Within each fitting
#'   window, the first element is `NA` and the remaining elements
#'   are such that `x[i]` is the number of cases observed from `t[i-1]`
#'   to `t[i]`.
#' }
#' \item{`wl`}{
#'   An integer vector of length `nlevels(index)` giving the length
#'   of each fitting window.
#' }
#' \item{`Xs`, `Xd`}{
#'   The fixed effects design matrix in sparse and dense format.
#'   If `sparse_X = TRUE`, then `Xd` is empty. If `sparse_X = FALSE`,
#'   then `Xs` is empty. Uses zero-sum contrasts.
#' }
#' \item{`Z`}{
#'   The random effects design matrix in sparse format, without
#'   contrasts.
#' }
#' \item{`ffr`, `rfr`}{
#'   Data frames containing the factors referenced by the elements of
#'   `attr(frame, "fixed_term_labels")` and
#'   `attr(frame, "random_term_labels")`. See [make_frame()].
#' }
#' \item{`fnl`, `rnl`}{
#'   Integer vectors listing the number of levels of each factor in
#'   `ffr` and `rfr`.
#' }
#' \item{`fid`, `rid`}{
#'   Indicator matrices with term labels for row names and parameter
#'   names (with prefixes) for column names. Element `[i,j]` is equal
#'   to 1 if and only if term `i` appears in the mixed effects model
#'   (fixed/random component) for parameter `j`.
#' }
#' \item{`curve_flag`, `distr_flag`, `excess_flag`, `sparse_X_flag`}{
#'   Integer flags referencing `curve`, `distr`, `excess`, and
#'   `sparse_X`.
#' }
#' \item{`predict_flag`}{
#'   An integer flag set equal to 0 so that prediction code is not run.
#' }
#' Additional integer elements of the form `j_link_parameter`
#' (e.g., `j_log_r`) give the index of parameter names (with prefixes)
#' in `colnames(rid)`. Indexing starts at 0, and the value -1 indicates
#' that the parameter does not appear in the incidence model being fit.
#'
#' @examples
#' example("make_frame", package = "epigrowthfit")
#' tmb_data <-
#'   epigrowthfit:::make_tmb_data(frame, curve, distr, excess, sparse_X = FALSE)
#'
#' @keywords internal
#' @importFrom Matrix sparseMatrix
make_tmb_data <- function(frame, curve, distr, excess, sparse_X) {
  a <- attributes(frame)

  ## Number of parameters
  p <- length(a$fixed_term_labels)

  ## Parameter names (with prefixes)
  pn <- get_par_names(curve, distr, excess, link = TRUE) # model
  pn0 <- get_par_names(link = TRUE) # all

  ## Compute time as number of days since start of fitting window
  date_to_integer <- function(x) as.integer(x - x[1L])
  dn <- a$date_name
  t <- unlist(lapply(unname(split(frame[[dn]], a$index)), date_to_integer))

  ## Set unused first observations to NA (avoids confusion)
  cn <- a$cases_name
  x <- frame[[cn]]
  x[!duplicated(a$index)] <- NA_integer_

  ## Find lengths of fitting windows
  wl <- as.vector(table(a$index))

  ## Get character vectors containing unique term labels
  ftl_incl_dupl <- unlist(a$fixed_term_labels, use.names = FALSE) # length `p`
  rtl_incl_dupl <- unlist(a$random_term_labels, use.names = FALSE)
  ftl <- unique(ftl_incl_dupl)
  rtl <- unique(rtl_incl_dupl)
  names(ftl) <- ftl
  names(rtl) <- rtl

  ## Evaluate term labels in `frame` and construct a data frame from
  ## the result
  ffr <- data.frame(lapply(ftl, get_factor, frame), check.names = FALSE)
  rfr <- data.frame(lapply(rtl, get_factor, frame), check.names = FALSE)

  ## Get number of levels for each factor
  fnl <- vapply(ffr, nlevels, integer(1L))[ftl_incl_dupl] # length `p`
  rnl <- vapply(rfr, nlevels, integer(1L))

  ## Get design matrix for each factor
  fmat <- lapply(ffr, factor_to_matrix, sparse = sparse_X, intercept = TRUE)[ftl_incl_dupl] # length `p`
  rmat <- lapply(rfr, factor_to_matrix, sparse = TRUE, intercept = FALSE)

  ## Combine design matrices columnwise
  if (sparse_X) {
    Xs <- do.call(cbind, fmat)
    Xd <- matrix(integer(0), nrow = 0L, ncol = ncol(Xs))
  } else {
    Xd <- do.call(cbind, fmat)
    Xs <- sparseMatrix(i = integer(0L), j = integer(0L), x = integer(0L), dims = c(0L, ncol(Xd)))
  }
  if (is.null(rtl_incl_dupl)) {
    Z <- sparseMatrix(integer(0L), integer(0L), x = integer(0L), dims = c(nrow(frame), 0L))
  } else {
    Z <- do.call(cbind, rmat)
  }

  ## Construct indicator matrices such that 1 at index (i,j)
  ## indicates that parameter j follows model using term i
  fid <- matrix(0L, nrow = length(ftl), ncol = p, dimnames = list(ftl, pn))
  for (i in seq_len(nrow(fid))) { # FIXME: outer(FUN = "==")?
    for (j in seq_len(ncol(fid))) {
      fid[i, j] <- 1L * (ftl[[i]] == a$fixed_term_labels[[j]])
    }
  }
  if (is.null(rtl)) {
    rid <- matrix(integer(0L), ncol = p, dimnames = list(NULL, pn))
  } else {
    rid <- matrix(0L, nrow = length(rtl), ncol = p, dimnames = list(rtl, pn))
    for (i in seq_len(nrow(rid))) { # FIXME: outer(FUN = "%in%")?
      for (j in seq_len(ncol(rid))) {
        rid[i, j] <- 1L * (rtl[[i]] %in% a$random_term_labels[[j]])
      }
    }
  }

  curve_list <- c("exponential", "subexponential", "gompertz",
                  "logistic", "richards")
  distr_list <- c("pois", "nbinom")

  l1 <- list(
    index = a$index,
    t = t,
    x = x,
    wl = wl,
    Xs = Xs,
    Xd = Xd,
    Z = Z,
    ffr = ffr,
    rfr = rfr,
    fnl = fnl,
    rnl = rnl,
    fid = fid,
    rid = rid,
    curve_flag = get_flag("curve", curve),
    distr_flag = get_flag("distr", distr),
    excess_flag = 1L * excess,
    sparse_X_flag = 1L * sparse_X,
    predict_flag = 0L
  )
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
#' \item{`a`}{1}
#' \item{`b`}{1}
#' \item{`nbdisp`}{1}
#' }
#' The naive estimates are log-transformed (all but `p`) or
#' logit-transformed (`p` only), and coefficients under zero-sum
#' contrasts are obtained from the within-group means to form
#' segments of parameter object `beta`.
#'
#' When `par_init = NULL` and there are random effects, `log_sd_b`
#' and `b` are for simplicity initialized as zero vectors.
#'
#' @return
#' A list with elements:
#' \item{`beta`}{
#'   A numeric vector of length `sum(tmb_data$fnl)` listing initial
#'   values for the fixed effects coefficients, namely between-group
#'   means and within-group offsets.
#' }
#' \item{`log_sd_b`}{
#'   If there are random effects, then a numeric vector of length
#'   `sum(tmb_data$rid)` listing initial values for the log standard
#'   deviations of the within-group random effects. Otherwise, `NA_real_`.
#' }
#' \item{`b`}{
#'   If there are random effects, then a numeric vector of length
#'   `sum(tmb_data$rnl * rowSums(tmb_data$rid))` listing initial
#'   values for the within-group random effects. Otherwise, `NA_real_`.
#' }
#'
#' @examples
#' example("make_tmb_data", package = "epigrowthfit")
#' tmb_parameters <-
#'   epigrowthfit:::make_tmb_parameters(tmb_data, curve, distr, excess,
#'                                      par_init = NULL)
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit
make_tmb_parameters <- function(tmb_data, curve, distr, excess, par_init) {
  re <- has_random(tmb_data)

  ## If user did not specify a parameter vector
  if (is.null(par_init)) {
    ts_split  <- split(data.frame(tmb_data[c("t", "x")]), tmb_data$index)
    ffr_split <- split(tmb_data$ffr, tmb_data$index)

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
    get_log_K     <- function(d) log(2) + log(sum(d$x[-1L], na.rm = TRUE))

    log_r_log_c0 <- vapply(ts_split, get_log_r_log_c0, numeric(2L))
    log_tinfl    <- vapply(ts_split, get_log_tinfl,    numeric(1L))
    log_K        <- vapply(ts_split, get_log_K,        numeric(1L))
    log_r  <- log_r_log_c0[1L, ]
    log_c0 <- log_r_log_c0[2L, ]
    p <- 0.8
    log_alpha <- switch(curve,
      subexponential = log_r + (1 - p) * log_c0,
      gompertz       = log_r - log(log_K - log_c0),
      NA
    )

    pfr <- data.frame(
      log_r      = log_r,
      log_alpha  = log_alpha,
      log_c0     = log_c0,
      log_tinfl  = log_tinfl,
      log_K      = log_K,
      logit_p    = log(p) - log1p(-p),
      log_a      = 0,
      log_b      = 0,
      log_nbdisp = 0
    )
    ffr <- do.call(rbind, lapply(ffr_split, "[", 1L, , drop = FALSE))

    pn <- get_par_names(curve, distr, excess, link = TRUE)
    beta <- unlist(lapply(pn, function(s) {
      tl <- names(tmb_data$fnl)[match(s, pn)]
      ## within-group means
      m <- vapply(split(pfr[[s]], ffr[[tl]]), mean, numeric(1L))
      ## mean of means
      mm <- mean(m)
      unname(c(mm, m[-length(m)] - mm))
    }))
    if (re) {
      log_sd_b <- rep(0, sum(tmb_data$rid))
      b <- rep(0, sum(tmb_data$rnl * rowSums(tmb_data$rid)))
    } else {
      log_sd_b <- NA_real_
      b <- NA_real_
    }

  ## Otherwise
  } else {
    p <- c(
      beta     = sum(tmb_data$fnl),
      log_sd_b = if (re) sum(tmb_data$rid),
      b        = if (re) sum(tmb_data$rnl * rowSums(tmb_data$rid))
    )
    stop_if_not(
      is.numeric(par_init),
      length(par_init) == sum(p),
      is.finite(par_init),
      m = sprintf("`par_init` must be a finite numeric vector of length %d.", sum(p))
    )

    par_init <- unname(par_init)
    beta <- par_init[seq_len(p[1L])]
    if (re) {
      log_sd_b <- par_init[p[1L] + seq_len(p[2L])]
      b <- par_init[p[1L] + p[2L] + seq_len(p[3L])]
    } else {
      log_sd_b <- NA_real_
      b <- NA_real_
    }
  }

  list(beta = beta, log_sd_b = log_sd_b, b = b)
}

#' Create TMB infrastructure
#'
#' Gathers necessary components of a call to [TMB::MakeADFun()].
#'
#' @param frame
#'   A model frame returned by [make_frame()].
#' @param par_init
#'   A user-specified full parameter vector, or otherwise `NULL`.
#'   See [egf()].
#' @inheritParams egf
#'
#' @return
#' A list with elements `data`, `parameters`, `map`, `random`,
#' `DLL`, and `silent`.
#'
#' @examples
#' example("make_frame", package = "epigrowthfit")
#' tmb_args <-
#'   epigrowthfit:::make_tmb_args(frame, curve, distr, excess,
#'                                sparse_X = FALSE, par_init = NULL)
#' tmb_out <- do.call(TMB::MakeADFun, tmb_args)
#'
#' @seealso [make_tmb_data()], [make_tmb_parameters()]
#' @keywords internal
make_tmb_args <- function(frame, curve, distr, excess,
                          sparse_X, par_init) {
  data <- make_tmb_data(frame, curve, distr, excess, sparse_X)
  parameters <- make_tmb_parameters(data, curve, distr, excess, par_init)
  if (has_random(data)) {
    map <- list()
    random <- "b"
  } else {
    map <- list(log_sd_b = factor(NA), b = factor(NA))
    random <- NULL
  }
  list(
    data = data,
    parameters = parameters,
    map = map,
    random = random,
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
#' @keywords internal
has_random <- function(object) {
  UseMethod("has_random", object)
}

#' @export
has_random.tmb_data <- function(object) {
  nrow(object$rid) > 0L
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
#' @examples
#' example("make_tmb_args", package = "epigrowthfit")
#' m <- c("nlminb", "nlm", "Nelder-Mead", "BFGS")
#' s <- c("par", "estimate", "par", "par")
#' f <- function(m, s) epigrowthfit:::optim_tmb_out(tmb_out, method = m)[[s]]
#' ## Not run
#' \dontrun{
#' pp <- mapply(f, m = m, s = s) # needs 1-2 min
#' }
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
#' s <- c("beta", "log_sd_b", "b")
#' n <- 10L
#' par <- numeric(n)
#' names(par) <- sample(s, n, replace = TRUE)
#' epigrowthfit:::rename_par(par)
#'
#' @keywords internal
#' @importFrom stats setNames
rename_par <- function(par) {
  for (s in c("beta", "log_sd_b", "b")) {
    i <- which(names(par) == s)
    names(par)[i] <- sprintf("%s[%d]", s, seq_along(i))
  }
  par
}

#' Undo zero-sum contrasts
#'
#' Multiplies each segment of a vector of fixed effects
#' coefficients in front by the appropriate matrix of contrasts
#' to recover within-group means.
#'
#' @param beta
#'   The segment of a TMB-generated parameter vector with names
#'   `"beta"` or, more generally, a numeric vector.
#' @param fnl
#'   An integer vector giving the length of each segment of `beta`,
#'   such that `sum(fnl) = length(beta)`. See [make_tmb_data()].
#'
#' @return
#' A named numeric vector of length `length(beta)`.
#'
#' @examples
#' fnl <- 1:4
#' beta <- sample(5L, sum(fnl), replace = TRUE)
#' epigrowthfit:::decontrast_beta(beta, fnl)
#'
#' @keywords internal
decontrast_beta <- function(beta, fnl) {
  names(beta) <- sprintf("beta[%d]", seq_along(beta))
  beta_split <- unname(split(beta, rep(seq_along(fnl), fnl)))

  extract_index <- function(s) {
    gsub("^beta\\[([0-9]+)\\]$", "\\1", s)
  }
  decontrast_beta_segment <- function(x) {
    l <- length(x)
    if (l == 1L) {
      return(x)
    }
    nx <- names(x)
    y <- c(x[1L] + x[-1L], x[1L] - sum(x[-1L])) # matrix multiply
    names(y)[-l] <- sprintf("%s+%s", nx[1L], nx[-1L])
    if (l == 2L) {
      names(y)[l] <- sprintf("%s-%s", nx[1L], nx[2L])
    } else {
      names(y)[l] <- sprintf("%s-sum(beta[%s:%s])",
                             nx[1L],
                             extract_index(nx[2L]),
                             extract_index(nx[l]))
    }
    y
  }
  unlist(lapply(beta_split, decontrast_beta_segment))
}

#' Describe elements of a parameter vector
#'
#' @param par
#'   A TMB-generated parameter vector, renamed using [rename_par()],
#'   so that elements have names of the form
#'   `"beta[i]"`, `"log_sd_b[j]"`, or `"b[k]"`.
#' @inheritParams make_tmb_parameters
#'
#' @return
#' A data frame with one row per element of `beta`, `log_sd_b`, `b`,
#' and `u = decontrast_beta(beta, tmb_data$fnl)` (in that order) and
#' six variables:
#' \item{`vector`}{
#'   A factor with levels `beta`, `log_sd_b`, `b`, and
#'   `.decontrast(beta)`, allowing the data frame to be split by vector.
#' }
#' \item{`name`}{
#'   A character vector concatenating `names(beta)`, `names(log_sd_b)`,
#'   `names(b)`, and `names(u)`.
#' }
#' \item{`estimate`}{
#'   A numeric vector concatenating `beta`, `log_sd_b`, `b`, and `u`.
#' }
#' \item{`response`}{
#'   A factor naming the relevant response variable, i.e., the response
#'   variable of the generalized linear mixed effects model to which
#'   the parameter belongs. For example, when fitting an exponential
#'   model, `levels(response)` will include `"log_r"` and `"log_c0"`.
#' }
#' \item{`term`}{
#'   A factor naming the relevant term of the mixed effects model
#'   formula.
#' }
#' \item{`level`}{
#'   A character vector. For elements of `u` and `b` (within-group
#'   means and within-group random effects), the appropriate element
#'   of `levels(term_label)` (i.e., the group label) is given. For
#'   elements of `beta` (between-group means and within-group offsets),
#'   `".mean(term_label)"` or `".offset(group_label)"` is given. For
#'   elements of `log_sd_b` (log standard deviations of within-group
#'   random effects), `".log_sd_b(term_label)"` is given.
#' }
#'
#' @keywords internal
get_par_info <- function(par, tmb_data) {
  with(tmb_data[c("ffr", "rfr", "fnl", "rnl", "rid")], {
    re <- has_random(tmb_data)
    pn <- get_par_names(tmb_data, link = TRUE)

    beta_no_contrasts <- decontrast_beta(
      beta = par[grep("beta", names(par))],
      fnl  = fnl
    )
    f <- function(x, s) {
      c(sprintf(".mean(%s)", s),
        sprintf(".offset(%s)", levels(x)[-nlevels(x)]))
    }

    d_vector <- c(
      sub("\\[[0-9]+\\]$", "", names(par)),
      rep(".decontrast(beta)", length(beta_no_contrasts))
    )
    d_name <- c(names(par), names(beta_no_contrasts))
    d_estimate <- unname(c(par, beta_no_contrasts))
    d_response <- c(
      rep(pn, fnl),
      if (re) rep(pn, nrow(rid))[t(rid) > 0L],
      if (re) rep(rep(pn, nrow(rid)), t(rnl * rid)),
      rep(pn, fnl)
    )
    d_term <- c(
      rep(names(fnl), fnl),
      if (re) rep(rownames(rid), rowSums(rid)),
      if (re) rep(rownames(rid), rowSums(rnl * rid)),
      rep(names(fnl), fnl)
    )
    d_level <- c(
      unlist(Map(f, ffr[names(fnl)], names(fnl))),
      if (re) rep(sprintf(".log_sd_b(%s)", rownames(rid)), rowSums(rid)),
      if (re) unlist(Map(rep, lapply(rfr[rownames(rid)], levels), rowSums(rid))),
      unlist(lapply(ffr[names(fnl)], levels))
    )

    data.frame(
      vector   = factor(d_vector, levels = unique(d_vector)),
      name     = d_name,
      estimate = d_estimate,
      response = factor(d_response, levels = unique(d_response)),
      term     = factor(d_term, levels = unique(d_term)),
      level    = d_level,
      stringsAsFactors = FALSE
    )
  })
}

#' Split up `ADREPORT()`ed variables
#'
#' Extracts reported variables from an `"sdreport"` object and
#' returns them in a list.
#'
#' @param sdreport
#'   An `"sdreport"` object returned by [TMB::sdreport()].
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
  lapply(split(as.data.frame(ssdr), rownames(ssdr)), "row.names<-", NULL)
}



