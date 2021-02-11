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
#'   A logical scalar. If `TRUE`, then a prefix indicating the
#'   link function used internally (either `"log_"` or `"logit_"`)
#'   is prepended to each parameter name.
#'
#' @return
#' If `link = FALSE`, then the default method returned the subset
#' of `c("r", "alpha", "c0", "tinfl", "K", "p", "a", "b", "nbdisp")`
#' relevant to `curve`, `distr`, and `excess`.
#'
#' @examples
#' pn  <- get_par_names(curve = "exponential", distr = "pois",
#'                      excess = FALSE, link = FALSE)
#' lpn <- get_par_names(curve = "exponential", distr = "pois",
#'                      excess = FALSE, link = TRUE)
#' identical(pn, sub("^(log|logit)_", "", lpn))
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
get_par_names.tmb_data <- function(curve, distr, excess, link = TRUE) {
  if (link) {
    colnames(curve$rid)
  } else {
    remove_link_string(colnames(curve$rid))
  }
}

#' @export
get_par_names.egf <- function(curve, distr, excess, link = TRUE) {
  get_par_names(curve$tmb_args$data, link = link)
}

#' Manipulate parameter names
#'
#' Add, remove, and extract prefixes from parameter names.
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
#' Retrieve the link or inverse link function corresponding to a string.
#'
#' @param s A character string, either `"log"` or `"logit"`.
#'
#' @return
#' If `s = "log"`, then `log` or its inverse `exp`.
#' If `s = "logit"`, then `function(p) qlogis(p)`
#' or its inverse `function(q) plogis(q)`.
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

######################## REFACTOR



check_formula_ts <- function(formula_ts) {
  stop_if_not(
    inherits(formula_ts, "formula"),
    m = "`formula_ts` must be a formula."
  )
  stop_if_not(
    length(formula_ts) == 3L,
    m = "`formula_ts` must have a response."
  )
  formula_ts[[2L]] <- deparen(formula_ts[[2L]])
  stop_if_not(
    is.name(formula_ts[[2L]]),
    m = "Left hand side of `formula_ts` must be a name."
  )
  formula_ts[[3L]] <- deparen(formula_ts[[3L]])
  if (is.name(formula_ts[[3L]])) {
    return(formula_ts)
  }
  if (is.call(formula_ts[[3L]]) && formula_ts[[3L]][[1L]] == as.name("|")) {
    formula_ts[[3L]][[2L]] <- deparen(formula_ts[[3L]][[2L]])
    stop_if_not(
      is.name(formula_ts[[3L]][[2L]]),
      m = "Left hand side of `|` in `formula_ts`\nmust be a name."
    )
    formula_ts[[3L]][[3L]] <- deparen(formula_ts[[3L]][[3L]])
    stop_if_not(
      is_interaction(formula_ts[[3L]][[3L]]),
      m = "Right hand side of `|` in `formula_ts`\nmust be a name or interaction."
    )
    formula_ts[[3L]][[3L]] <- expand_terms(formula_ts[[3L]][[3L]])
    return(formula_ts)
  }
  stop("Right hand side of `formula_ts`\nmust be a name or a call to `|`.")
}

check_formula_glmm <- function(formula_glmm, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (inherits(formula_glmm, "formula") && length(formula_glmm) == 2L) {
    formula_glmm[[3L]] <- formula_glmm[[2L]]
    formula_glmm <- Map(
      f = function(s, x) `[[<-`(x, 2L, as.name(s)),
      s = pn,
      x = rep.int(list(formula_glmm), p)
    )
  } else if (inherits(formula_glmm, "list") &&
             length(formula_glmm) > 0L &&
             all(vapply(formula_glmm, inherits, FALSE, "formula")) &&
             all(lengths(formula_glmm) == 3L)) {
    formula_glmm <- lapply(formula_glmm, function(x) `[[<-`(x, 2L, deparen(x[[2L]])))
    lhs <- lapply(formula_glmm, `[[`, 2L)
    stop_if_not(
      vapply(lhs, is.name, FALSE),
      m = "Left hand side of each `formula_glmm` formula\nmust be a name."
    )
    lhs_as_character <- vapply(lhs, all.vars, "")
    stop_if_not(
      lhs_as_character %in% pn,
      !duplicated(lhs_as_character),
      m = paste0(
        "Name on left hand side of each `formula_glmm` formula\n",
        "must be an element of:\n",
        "`get_par_names(curve, distr, excess, link = TRUE)`."
      )
    )
    names(formula_glmm) <- lhs_as_character
    make_default_formula <- function(s) {
      as.formula(call("~", as.name(s), 1))
    }
    ss <- setdiff(pn, lhs_as_character)
    formula_glmm[ss] <- lapply(ss, make_default_formula)
    formula_glmm <- formula_glmm[pn]
  } else {
    stop("`formula_glmm` must be a formula of the form `~rhs`\n",
         "or a list of formulae of the form `par ~ rhs`.")
  }
  lapply(formula_glmm, function(x) {
    ftrt <- split_effects(x)
    formula_glmm[[3L]] <- unsplit_terms(do.call(c, ftrt))
    structure(formula_glmm,
      fixed_terms = ftrt$fixed_terms,
      random_terms = ftrt$random_terms
    )
  })
}

split_effects <- function(x) {
  if (!(inherits(x, "formula") || is.call(x) || is.name(x) || is.numeric(x))) {
    stop("`x` must be a formula, call, name, or number.")
  }
  orig_terms <- split_terms(x)
  has_bar <- vapply(orig_terms, function(x) is.call(x) && x[[1L]] == as.name("|"), FALSE)

  if (all(has_bar)) {
    fixed_terms <- list()
  } else {
    fixed_terms <- split_terms(expand_terms(unsplit_terms(orig_terms[!has_bar])))
  }

  if (any(has_bar)) {
    lhs <- lapply(orig_terms[has_bar], function(x) deparen(x[[2L]]))
    stop_if_not(
      vapply(lhs, function(x) is.name(x) || x == 1, FALSE),
      m = "Left hand side of each `|` in `formula_glmm`\nmust be 1 or a name."
    )
    rhs <- lapply(orig_terms[has_bar], function(x) expand_terms(x[[3L]]))
    rhs_split <- lapply(rhs, split_terms)
    stop_if_not(
      vapply(do.call(c, rhs_split), is_interaction, FALSE),
      m = "Right hand side of each `|` in `formula_glmm`\nmust expand to a sum of names and interactions."
    )
    random_terms <- Map(
      f = function(x, y) call("|", x, y),
      x = rep.int(lhs, lengths(rhs_split)),
      y = do.call(c, rhs_split)
    )
  } else {
    random_terms <- list()
  }

  list(fixed_terms = fixed_terms, random_terms = random_terms)
}




#' Construct a model frame
#'
#' Constructs a model frame to be used by [egf()],
#' while performing myriad checks on the arguments.
#'
#' @inheritParams egf
#'
#' @details
#' Rows of `data` belonging to time series without fitting windows
#' are discarded. Rows of `data[-discarded, ]` belonging to fitting
#' windows are returned. The remaining rows of `data[-discarded, ]`
#' are preserved in attribute `extra`. These rows are kept because
#' it is desirable for [plot.egf()] to display entire time series,
#' not only segments to which a curve has been fit.
#'
#' @return
#' A data frame containing the variables named in `formula_ts` and
#' `formula_glmm`, with attributes:
#' \item{`extra`}{
#'   A data frame preserving unused data (see Details).
#' }
#' \item{`window`}{
#'   A factor of length `nrow(frame)` such that `split(frame, window)`
#'   splits the data frame by fitting window. Not usually identical
#'   to the so-named argument (see Details).
#' }
#' \item{`formula_ts`, `formula_glmm`}{
#'   Processed versions of the so-named arguments.
#' }
#'
#' @examples
#' r <- log(2) / 10
#' c0 <- 10
#' time <- 0:99
#' mu <- diff(c0 * exp(r * time))
#'
#' date <- rep(as.Date(time, origin = "2020-01-01"), 4)
#' cases <- c(replicate(4, c(NA, rpois(length(mu), mu))))
#' continent <- factor(rep(c("asia", "europe"), each = 200))
#' country <- factor(rep(c("china", "japan", "france", "germany"), each = 100))
#' wave <- factor(rep(rep(c(1, 2), each = 50), 4))
#'
#' data <- data.frame(date, cases, continent, country, wave)
#'
#' x <- c(NA, 1, NA, 2, NA,
#'        NA, 3, NA, 4, NA,
#'        NA, 5, NA, 6, NA,
#'        NA, 7, NA, 8, NA)
#' index <- factor(rep(x, each = 20))
#'
#' curve <- "exponential"
#' distr <- "pois"
#' excess <- FALSE
#' get_par_names(curve, distr, excess, link = FALSE)
#' ## [1] "r" "c0"
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
#' group_by <- ~country
#'
#' frame <- epigrowthfit:::make_frame(
#'   formula, fixed, random, group_by,
#'   data, index,
#'   curve, distr, excess,
#'   na_action = "pass"
#' )
#'
#' @keywords internal
#' @importFrom stats complete.cases
make_frame <- function(formula_ts, formula_glmm, data, window,
                       curve, distr, excess, na_action) {
  ## Check model formulae
  formula_ts   <- check_formula_ts(formula_ts)
  formula_glmm <- check_formula_glmm(formula_glmm, curve, distr, excess)
  a <- attributes(formula_glmm)

  ## Check that formula variables can be found
  stop_if_not(
    inherits(data, c("data.frame", "list", "environment")),
    m = "`data` must be a data frame, list, or environment."
  )
  if (is.environment(data)) {
    data <- as.list(data)
  }
  av <- unique(c(all.vars(formula_ts), all.vars(formula_glmm)))
  found <- (av %in% names(data))
  stop_if_not(
    all(found),
    m = paste0(
      "Formula variables not found in `data`:\n",
      paste(sprintf("`%s`", av[!found]), collapse = ", ")
    )
  )

  ## Construct model frame
  stop_if_not(
    vapply(data[av], is.atomic, FALSE),
    m = "Formula variables must be atomic."
  )
  if (!is.data.frame(data)) {
    stop_if_not(
      diff(range(lengths(data[av]))) == 0L,
      m = "Formula variables must have a common length."
    )
    for (s in av) {
      dim(data[[s]]) <- NULL
    }
  }
  frame <- as.data.frame(data[av])
  N <- nrow(frame)
  min_window_length <- 7L
  stop_if_not(
    N >= min_window_length,
    m = sprintf("Model frame must have at least %d rows.", min_window_length)
  )

  ## Check data types
  if (is.call(formula_ts[[3L]])) {
    gv <- all.vars(formula_ts[[3L]][[3L]])
  } else {
    gv <- character(0L)
  }
  dv <- setdiff(all.vars(formula_ts[[3L]]), gv)
  stop_if_not(
    inherits(frame[[dv]], "Date"),
    m = sprintf("`%s` must inherit from class \"Date\".", dv)
  )
  cv <- all.vars(formula_ts[[2L]])
  stop_if_not(
    is.numeric(frame[[cv]]),
    all(frame[[cv]] >= 0, na.rm = TRUE),
    m = sprintf("`%s` must be a non-negative numeric vector.", cv)
  )
  rv_lhs <- unique(unlist(lapply(a$random_terms, function(x) all.vars(x[[2L]]))))
  stop_if_not(
    vapply(frame[rv_lhs], is.numeric, FALSE),
    m = "Formula variables on left hand side of `|`\nmust be numeric vectors."
  )
  rv_rhs <- unique(unlist(lapply(a$random_terms, function(x) all.vars(x[[3L]]))))
  stop_if_not(
    vapply(frame[c(rv_rhs, cv)], is.factor, FALSE),
    m = "Formula variables on right hand side of `|`\nmust be factors."
  )
  fv <- unique(unlist(lapply(a$fixed_terms, all.vars)))
  stop_if_not(
    vapply(frame[fv], function(x) is.factor(x) || is.numeric(x), FALSE),
    m = "Fixed effects formula variables\nmust be factors or numeric vectors."
  )

  ## Drop unused and NA factor levels
  facv <- an[vapply(frame[an], is.factor, FALSE)]
  frame[facv] <- lapply(frame[facv], droplevels, exclude = NA)

  ## Check missing values
  stop_if_not(
    !vapply(frame[setdiff(av, cv)], anyNA, FALSE),
    m = sprintf("Formula variables other than incidence (`%s`)\ncannot have missing values.", cv)
  )

  ## Check fitting windows
  stop_if_not(
    is.factor(window),
    length(window) == N,
    m = sprintf("`window` must be a factor of length `length(%s)`.", dv)
  )
  window <- droplevels(window, exclude = NA)
  stop_if_not(
    nlevels(window) > 0L,
    m = "`window` must have at least one nonempty level."
  )
  stop_if_not(
    tabulate(window) >= min_window_length,
    m = sprintf("Nonempty levels of `window` must have\nfrequency %d or greater.", min_window_length)
  )

  if (length(facv) > 0L) {
    ## Check factors
    stop_if_not(
      nlevels(interaction(cbind(window, frame[facv]), drop = TRUE)) == nlevels(window),
      m = "Factors must be constant in each level of `window`."
    )
    if (length(facv) > length(gv)) {
      frame_red <- droplevels(frame[!is.na(window) & !duplicated(window), setdiff(facv, gv), drop = FALSE])
      stop_if_not(
        vapply(frame_red, nlevels, 0L) > 1L,
        m = "Factors in `formula_glmm` must not be constant\nacross all levels of `window`."
      )
    }
  }

  if (length(gv) > 0L) {
    ## Discard rows belonging to time series without fitting windows
    frame_red <- droplevels(frame[!is.na(window) & !duplicated(window), gv, drop = FALSE])
    frame[gv] <- Map(factor,
      x = frame[gv],
      levels = lapply(frame_red[gv], levels)
    )
    keep <- complete.cases(frame[gv])
    frame <- droplevels(frame[keep, , drop = FALSE])
    index <- index[keep]
    N <- nrow(frame)

    ## Order rows by time series
    ts <- interaction(frame[gv], drop = TRUE)
    ord <- order(ts)
    frame <- frame[ord, ]
    index <- index[ord]
    ts <- ts[ord]

    ## Check date
    date_split <- split(frame[[dv]], ts)
    stop_if_not(
      vapply(date_split, function(x) all(diff(x) > 0), FALSE),
      m = sprintf("`%s` must be increasing in each level of `%s`.", dv, paste(gv, collapse = ":"))
    )
  } else {
    ## Check date
    stop_if_not(
      all(diff(frame[[dv]]) > 0),
      m = sprintf("`%s` must be increasing.", dv)
    )
  }

  ## Check incidence
  if (!is.integer(frame[[cv]])) {
    which_is_infinite <- which(!is.finite(frame[[cv]]))
    if (length(which_is_infinite) > 0L) {
      warning(sprintf("Non-finite numeric elements of `%s` replaced with NA.", cv))
      frame[[cv]][which_is_infinite] <- NA
    }
    which_is_non_integer <- which(frame[[cn]] %% 1 != 0)
    if (length(which_is_non_integer) > 0L) {
      warning(sprintf("Non-integer numeric elements of `%s` truncated.", cn))
      frame[[cn]][which_is_non_integer] <- trunc(frame[[cn]][which_is_non_integer])
    }
  }
  cases_split <- split(frame[[cn]], window)
  if (na_action == "fail") {
    stop_if_not(
      !vapply(cases_split, function(x) anyNA(x[-1L]), FALSE),
      m = sprintf("na_action = \"fail\": `%s` has missing values\nin at least one fitting window.", cv)
    )
  } else {
    stop_if_not(
      vapply(cases_split, function(x) sum(!is.na(x)) >= min_window_length, FALSE),
      m = sprintf("`%s` has insufficient data (fewer than %d observations)\nin at least one fitting window.", cv, min_window_length)
    )
  }

  ## Check that fitting windows are contiguous
  stop_if_not(
    vapply(split(seq_len(N), window), function(i) all(diff(i) == 1L), FALSE),
    m = "Fitting windows must be contiguous."
  )

  ## Clean up
  row.names(frame) <- NULL

  # Variable/level order is a dependency of multiple methods
  # (for better or for worse), so edit with care
  frame <- frame[c(dv, cv, setdiff(av, c(dv, cv)))] # put time series first, everything else after
  window <- factor(window, levels = unique(window), exclude = NA) # order levels by occurrence

  ## Subset rows belonging to a fitting window
  ## and preserve the difference as an attribute
  ## for use by plot.egf()
  keep <- !is.na(index)
  structure(droplevels(frame[keep, ]),
    extra = droplevels(frame[!keep, ]),
    window = window[keep],
    formula_ts = formula_ts,
    formula_glmm = formula_glmm
  )
}

#' Get term labels from a formula
#'
#' A wrapper for `labels(terms())` with special handling of certain
#' formulae.
#'
#' @param f A formula of the form `~rhs`.
#'
#' @return
#' If `f` has no variable terms, as in `f = ~1`, then `"1"`.
#' Otherwise, the result of `labels(terms(f))` *after* replacing
#' terms of `f` of the form `(1 | rhs_of_bar)` with `rhs_of_bar`.
#'
#' @examples
#' get_term_labels(~w + x + y:z)
#' get_term_labels(~1)
#' get_term_labels(~w*x + y/z) # crosses, nests expanded
#' get_term_labels(~w*x + (1 | y/z)) # bars ignored
#'
#' @noRd
#' @importFrom stats terms reformulate
get_term_labels <- function(f) {
  if (is.null(f)) {
    return(NULL)
  }
  tl <- labels(terms(f))
  if (length(tl) == 0L) {
    return("1")
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
#' In the special case `term_label = "1"`, the result is
#' `rep(factor(1), nrow(frame))`.
#'
#' @examples
#' f <- function() factor(sample(5L, 20L, replace = TRUE))
#' d <- data.frame(x = f(), y = f())
#' get_factor("x:y", frame = d)
#' get_factor("1", frame = d)
#'
#' @noRd
get_factor <- function(term_label, frame) {
  if (is.null(term_label)) {
    return(NULL)
  }
  if (term_label == "1") {
    return(rep(factor(1), nrow(frame)))
  }
  s <- unique(strsplit(term_label, ":")[[1L]])
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
#'   A logical scalar. If TRUE, then the matrix is returned in
#'   sparse format.
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
#' where `d = 0` if `levels(x)` contains `NA` and `d = sum(is.na(x))`
#' otherwise.
#'
#' @examples
#' x <- factor(sample(5L, 20L, replace = TRUE))
#' factor_to_matrix(x, sparse = TRUE, intercept = TRUE)
#' factor_to_matrix(x, sparse = TRUE, intercept = FALSE)
#' factor_to_matrix(x, sparse = FALSE, intercept = TRUE)
#' factor_to_matrix(x, sparse = FALSE, intercept = FALSE)
#'
#' @noRd
#' @importFrom Matrix sparseMatrix sparse.model.matrix
#' @importFrom stats model.matrix
factor_to_matrix <- function(x, sparse, intercept) {
  if (is.null(x)) {
    return(NULL)
  }
  d <- if (anyNA(levels(x))) 0L else sum(is.na(x))
  m <- length(x) - d
  n <- nlevels(x)
  if (m == 0L || n < 2L) {
    mn <- m * n
    if (sparse) {
      X <- sparseMatrix(i = seq_len(mn), j = rep.int(1L, mn),
                        x = rep.int(1L, mn), dims = c(m, n))
    } else {
      X <- rep.int(1L, mn)
      dim(X) <- c(m, n)
    }
  } else {
    f <- if (sparse) sparse.model.matrix else model.matrix
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
#'   An integer vector of length `nrow(frame)` giving time as
#'   a number of days since the earliest date in the relevant
#'   fitting window.
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
#'   Data frames containing the factors referenced by the elements
#'   of `frame` attributes `fixed_term_labels` and `random_term_labels`.
#'   See [make_frame()].
#' }
#' \item{`fnl`, `rnl`}{
#'   Integer vectors listing the number of levels of each factor in
#'   `ffr` and `rfr`.
#' }
#' \item{`fid`, `rid`}{
#'   Indicator matrices with term labels for row names and parameter
#'   names (with prefixes) for column names. Element `[i, j]` is equal
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
#' (e.g., `j_log_r`) give the index of response variables
#' (e.g., `"log_r"`) in `colnames(rid)`. Indexing starts at 0,
#' and the value -1 indicates that the parameter does not
#' belong to the incidence model being fit.
#'
#' @examples
#' example("make_frame", package = "epigrowthfit")
#' tmb_data <- epigrowthfit:::make_tmb_data(
#'   frame,
#'   curve, distr, excess,
#'   sparse_X = FALSE
#' )
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

  ## Compute time as number of days since start of fitting window.
  ## No issues with permutation here as `index` is ordered by level.
  date_to_integer <- function(x) days(x, since = x[1L])
  t <- unlist(lapply(split(frame[[1L]], a$index), date_to_integer), use.names = FALSE)

  ## Set unused first observations to NA. They are unused regardless,
  ## and this makes that explicit.
  x <- frame[[2L]]
  x[!duplicated(a$index)] <- NA_integer_

  ## Find lengths of fitting windows
  wl <- tabulate(a$index)

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
  fnl <- vapply(ffr, nlevels, 0L)[ftl_incl_dupl] # length `p`
  rnl <- vapply(rfr, nlevels, 0L)

  ## Get design matrix for each factor
  fmat <- lapply(ffr, factor_to_matrix, sparse = sparse_X, intercept = TRUE)[ftl_incl_dupl] # length `p`
  rmat <- lapply(rfr, factor_to_matrix, sparse = TRUE, intercept = FALSE)

  ## Combine design matrices columnwise
  if (sparse_X) {
    Xs <- do.call(cbind, fmat)
    Xd <- integer(0L)
    dim(Xd) <- c(0L, ncol(Xs))
  } else {
    Xd <- do.call(cbind, fmat)
    Xs <- sparseMatrix(i = integer(0L), j = integer(0L),
                       x = integer(0L), dims = c(0L, ncol(Xd)))
  }
  if (is.null(rtl_incl_dupl)) {
    Z <- sparseMatrix(i = integer(0L), j = integer(0L),
                      x = integer(0L), dims = c(nrow(frame), 0L))
  } else {
    Z <- do.call(cbind, rmat)
  }

  ## Construct indicator matrices such that 1 at index (i,j)
  ## indicates that parameter j follows model using term i
  fid <- matrix(0L, nrow = length(ftl), ncol = p, dimnames = list(ftl, pn))
  for (i in seq_len(nrow(fid))) { # FIXME: outer(FUN = `==`)?
    for (j in seq_len(ncol(fid))) {
      fid[i, j] <- 1L * (ftl[[i]] == a$fixed_term_labels[[j]])
    }
  }
  if (is.null(rtl)) {
    rid <- matrix(integer(0L), ncol = p, dimnames = list(NULL, pn))
  } else {
    rid <- matrix(0L, nrow = length(rtl), ncol = p, dimnames = list(rtl, pn))
    for (i in seq_len(nrow(rid))) { # FIXME: outer(FUN = `%in%`)?
      for (j in seq_len(ncol(rid))) {
        rid[i, j] <- 1L * (rtl[[i]] %in% a$random_term_labels[[j]])
      }
    }
  }

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
#' @param tmb_data A `"tmb_data"` object returned by [make_tmb_data()].
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
#'   A numeric vector of length `sum(tmb_data$fnl)` listing initial
#'   values for the fixed effects coefficients, each of which represents
#'   the mean of a set of within-group means or a within-group offset
#'   from that mean.
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
#' tmb_parameters <- epigrowthfit:::make_tmb_parameters(
#'   tmb_data,
#'   curve, distr, excess,
#'   par_init = NULL
#' )
#'
#' @keywords internal
#' @importFrom stats coef lm na.omit qlogis
make_tmb_parameters <- function(tmb_data, curve, distr, excess,
                                par_init) {
  any_re <- has_random(tmb_data)

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

    log_r_log_c0 <- vapply(ts_split, get_log_r_log_c0, rep(0, 2L))
    log_tinfl    <- vapply(ts_split, get_log_tinfl,    0)
    log_K        <- vapply(ts_split, get_log_K,        0)
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
      logit_p    = qlogis(p),
      log_a      = 0,
      log_b      = 0,
      log_nbdisp = 0
    )
    ffr <- do.call(rbind, lapply(ffr_split, `[`, 1L, , drop = FALSE))

    pn <- get_par_names(curve, distr, excess, link = TRUE)
    beta <- unlist(lapply(pn, function(s) {
      tl <- names(tmb_data$fnl)[match(s, pn)]
      ## within-group means
      m <- vapply(split(pfr[[s]], ffr[[tl]]), mean, 0, USE.NAMES = FALSE)
      ## mean of means
      mm <- mean(m)
      c(mm, m[-length(m)] - mm)
    }))
    if (any_re) {
      log_sd_b <- rep.int(0, sum(tmb_data$rid))
      b <- rep.int(0, sum(tmb_data$rnl * rowSums(tmb_data$rid)))
    } else {
      log_sd_b <- NA_real_
      b <- NA_real_
    }

  ## Otherwise
  } else {
    p <- c(
      beta     = sum(tmb_data$fnl),
      log_sd_b = if (any_re) sum(tmb_data$rid),
      b        = if (any_re) sum(tmb_data$rnl * rowSums(tmb_data$rid))
    )
    stop_if_not(
      is.numeric(par_init),
      length(par_init) == sum(p),
      is.finite(par_init),
      m = sprintf("`par_init` must be a finite numeric vector of length %d.", sum(p))
    )
    if (!is.null(names(par_init))) {
      names(par_init) <- NULL
    }
    beta <- par_init[seq_len(p[1L])]
    if (any_re) {
      log_sd_b <- par_init[seq.int(to = sum(p[1:2]), length.out = p[2L])]
      b <- par_init[seq.int(to = sum(p), length.out = p[3L])]
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
#' tmb_args <- epigrowthfit:::make_tmb_args(
#'   frame,
#'   curve, distr, excess,
#'   sparse_X = FALSE,
#'   par_init = NULL
#' )
#' tmb_out <- do.call(TMB::MakeADFun, tmb_args)
#'
#' @seealso [make_tmb_data()], [make_tmb_parameters()]
#' @keywords internal
make_tmb_args <- function(frame, curve, distr, excess,
                          sparse_X, par_init) {
  data <- make_tmb_data(
    frame = frame,
    curve = curve,
    distr = distr,
    excess = excess,
    sparse_X = sparse_X
  )
  parameters <- make_tmb_parameters(
    tmb_data = data,
    curve = curve,
    distr = distr,
    excess = excess,
    par_init = par_init
  )
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
#' @noRd
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
#' f <- function(m, s) optim_tmb_out(tmb_out, method = m)[[s]]
#' ## Not run
#' \dontrun{
#' pp <- Map(f,
#'   m = c("nlminb", "nlm", "Nelder-Mead", "BFGS")
#'   s = c("par", "estimate", "par", "par")
#' )
#' }
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
