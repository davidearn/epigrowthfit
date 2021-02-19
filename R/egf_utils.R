#' Get parameter names
#'
#' Returns the names used internally for parameters of the specified
#' incidence model.
#'
#' @param curve,distr,excess,weekday
#'   Character or logical flags specifying an incidence model.
#'   See [egf()]. Alternatively, `curve` can be an `"egf"` object
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
#' pn  <- get_par_names(curve = "exponential", distr = "pois",
#'                      excess = FALSE, weekday = FALSE, link = FALSE)
#' lpn <- get_par_names(curve = "exponential", distr = "pois",
#'                      excess = FALSE, weekday = FALSE, link = TRUE)
#' identical(pn, sub("^(log|logit)_", "", lpn))
#'
#' @export
get_par_names <- function(curve = NULL,
                          distr = NULL,
                          excess = NULL,
                          weekday = NULL,
                          link = TRUE) {
  UseMethod("get_par_names", curve)
}

#' @export
get_par_names.default <- function(curve = NULL,
                                  distr = NULL,
                                  excess = NULL,
                                  weekday = NULL,
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

#' Validate time series formula
#'
#' Checks that [egf()] argument `formula_ts` is a valid formula
#' of the form `y ~ x` or `y ~ x | group`.
#'
#' @inheritParams egf
#'
#' @return
#' `formula_ts`, simplified if possible.
#'
#' @noRd
check_formula_ts <- function(formula_ts) {
  stop_if_not(
    inherits(formula_ts, "formula"),
    m = "`formula_ts` must be a formula."
  )
  stop_if_not(
    length(formula_ts) == 3L,
    m = "`formula_ts` must have a response."
  )
  formula_ts <- deparen(formula_ts)
  lhs <- formula_ts[[2L]]
  stop_if_not(
    is.name(lhs),
    m = "Left hand side of `formula_ts` must be a name."
  )
  rhs <- formula_ts[[3L]]
  if (is.name(rhs)) {
    return(formula_ts)
  }
  stop_if_not(
    is.call(rhs) && rhs[[1L]] == as.name("|"),
    m = "Right hand side of `formula_ts`\nmust be a name or a call to `|`."
  )
  stop_if_not(
    is.name(rhs[[2L]]),
    m = "Left hand side of `|` in `formula_ts` must be a name."
  )
  stop_if_not(
    is_interaction(rhs[[3L]]),
    m = "Right hand side of `|` in `formula_ts`\nmust be a name or an interaction of names."
  )
  rhs[[3L]] <- expand_terms(rhs[[3L]])
  formula_ts[[3L]] <- rhs
  formula_ts
}

#' Validate mixed effects formulae
#'
#' Checks that [egf()] argument `formula_glmm` is a valid formula
#' of the form `~terms` or a valid list of formulae of the form
#' `par ~ terms`.
#'
#' @inheritParams egf
#'
#' @details
#' If `formula_glmm` is a formula, then it is recycled to a list.
#' If it is an incomplete list of formulae, then missing formulae
#' are assigned `par ~ 1`. Terms using formula operators such as
#' `"*"` are expanded in the result.
#'
#' @return
#' A list of formulae of the form `par ~ terms`, one for each
#' nonlinear model parameter. Each has attributes `fixed_terms`
#' and `random_terms`, which split the fixed and random effects
#' components of `terms` at `"+"`.
#'
#' @noRd
check_formula_glmm <- function(formula_glmm, curve, distr, excess, weekday) {
  pn <- get_par_names(curve = curve, distr = distr,
                      excess = excess, weekday = weekday, link = TRUE)
  make_default_formula <- function(s) as.formula(call("~", as.name(s), 1))

  if (inherits(formula_glmm, "formula") &&
      length(formula_glmm) == 2L) {
    formula_glmm <- check_formula_terms(formula_glmm)
    formula_glmm[[3L]] <- formula_glmm[[2L]]
    formula_glmm <- Map(
      f = function(s, x) {x[[2L]] <- as.name(s); x},
      s = pn,
      x = list(formula_glmm)
    )
  } else if (inherits(formula_glmm, "list") &&
             length(formula_glmm) > 0L &&
             all(vapply(formula_glmm, inherits, FALSE, "formula")) &&
             all(lengths(formula_glmm) == 3L)) {
    formula_glmm <- lapply(formula_glmm, deparen)
    lhs <- lapply(formula_glmm, `[[`, 2L)
    lhs_as_character <- vapply(lhs, deparse, "")
    stop_if_not(
      vapply(lhs, is.name, FALSE),
      lhs_as_character %in% pn,
      !duplicated(lhs_as_character),
      m = paste0(
        "Left hand side of each `formula_glmm` formula\n",
        "must be a name in:\n",
        "`get_par_names(curve, distr, excess, weekday, link = TRUE)`."
      )
    )
    names(formula_glmm) <- lhs_as_character
    fill <- setdiff(pn, lhs_as_character)
    formula_glmm[fill] <- lapply(fill, make_default_formula)
    formula_glmm <- lapply(formula_glmm[pn], check_formula_glmm_terms)
  } else {
    stop("`formula_glmm` must be a formula of the form `~terms`\n",
         "or a list of formulae of the form `par ~ terms`.")
  }

  formula_glmm
}

check_formula_glmm_terms <- function(x) {
  stop_if_not(
    inherits(x, "formula"),
    m = "`x` must be a formula."
  )
  x[[length(x)]] <- expand_terms(x)
  terms <- split_terms(x)
  has_bar <- vapply(terms, function(x) is.call(x) && x[[1L]] == as.name("|"), FALSE)
  stop_if_not(
    vapply(terms[!has_bar], is_interaction, FALSE),
    m = paste0(
      "Fixed effects formulae must expand to\n",
      "a sum of names and interactions of names."
    )
  )
  if (any(has_bar)) {
    lhs <- lapply(terms[has_bar], `[[`, 2L)
    stop_if_not(
      vapply(lhs, `==`, FALSE, 1) | vapply(lhs, is_interaction, FALSE),
      m = paste0(
        "Left hand side of `|` in `formula_glmm`\n",
        "must expand to 1 or a sum of names and\n",
        "interactions of names."
      )
    )
    rhs <- lapply(terms[has_bar], `[[`, 3L)
    stop_if_not(
      vapply(rhs, is_interaction, FALSE),
      m = paste0(
        "Right hand side of `|` in `formula_glmm`\n",
        "must expand to a sum of names and\n",
        "interactions of names."
      )
    )
  }

  structure(x,
    fixed_terms = terms[!has_bar],
    random_terms = terms[has_bar]
  )
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
#' are preserved in attribute `extra`.
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
                       curve, distr, excess, weekday, na_action) {
  ## Check model formulae
  formula_ts <- check_formula_ts(formula_ts)
  formula_glmm <- check_formula_glmm(formula_glmm, curve, distr, excess)
  ft <- do.call(c, lapply(formula_glmm, attr, "fixed_terms"))
  rt <- do.call(c, lapply(formula_glmm, attr, "random_terms"))

  ## Check that formula variables can be found
  stop_if_not(
    inherits(data, c("data.frame", "list", "environment")),
    m = "`data` must be a data frame, list, or environment."
  )
  if (is.environment(data)) {
    data <- as.list(data)
  }
  av <- unique(unlist(lapply(c(list(formula_ts), ft, rt), all.vars)))
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
  min_window_length <- 1L + length(get_par_names(curve, distr, excess, weekday))
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
  fv <- unique(unlist(lapply(ft, all.vars)))
  stop_if_not(
    vapply(frame[fv], function(x) is.factor(x) || is.numeric(x), FALSE),
    m = "Fixed effects formula variables\nmust be factors or numeric vectors."
  )
  rv_lhs <- unique(unlist(lapply(rt, function(x) all.vars(x[[2L]]))))
  stop_if_not(
    vapply(frame[rv_lhs], is.numeric, FALSE),
    m = "Formula variables on left hand side of `|`\nmust be numeric vectors."
  )
  rv_rhs <- unique(unlist(lapply(rt, function(x) all.vars(x[[3L]]))))
  stop_if_not(
    vapply(frame[c(rv_rhs, gv)], is.factor, FALSE),
    m = "Formula variables on right hand side of `|`\nmust be factors."
  )

  ## Drop unused and NA factor levels
  facv <- av[vapply(frame[av], is.factor, FALSE)]
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

  ## Variable/level order is a dependency of multiple methods
  ## (for better or for worse), so edit with care
  frame <- frame[c(dv, cv, setdiff(av, c(dv, cv)))] # order variables with time series first, everything else after
  window <- factor(window, levels = unique(window), exclude = NA) # order levels by occurrence

  ## Clean up
  row.names(frame) <- NULL

  ## Condense, e.g., `x + y` to `x:y` if `x` and `y` are factors
  condense_fixed_terms <- function(x) {
    fixed_terms <- attr(x, "fixed_terms")
    if (length(fixed_terms) == 1L && fixed_terms[[1L]] == 1) {
      return(x)
    }
    av <- lapply(fixed_terms, all.vars)
    av_numeric <- lapply(av, function(s) sort(s[vapply(frame[s], is.numeric, FALSE)]))
    av_factor <- Map(setdiff, av, av_numeric)
    l <- tapply(
      X = av_factor,
      INDEX = vapply(av_numeric, paste, "", collapse = ":"),
      FUN = function(l) paste(unique(unlist(l)), collapse = ":"),
      simplify = FALSE
    )
    l <- Map(
      f = function(s1, s2) strsplit(paste(s1, s2, sep = ":"), ":")[[1L]],
      s1 = names(l),
      s2 = l,
      USE.NAMES = FALSE
    )
    l <- lapply(l, function(s) {
      tt <- as.name(s[1L])
      for (i in seq_along(s)[-1L]) {
        tt <- call(":", tt, as.name(s[i]))
      }
      tt
    })
    if (all(lengths(av_numeric) > 0L)) {
      l <- c(list(1), l)
    }
    structure(unsplit_terms(c(l, a$random_terms)),
      fixed_terms = l,
      random_terms = a$random_terms
    )
  }
  formula_glmm <- lapply(formula_glmm, condense_fixed_terms)

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

#' Construct a design matrix from a formula term
#'
#' A wrapper for
#' [stats::model.matrix()] and [Matrix::sparse.model.matrix()]
#' operating on irreducible formula terms.
#'
#' @param x
#'   A formula term, one of
#'   (1) 1, as in `(~1)[[2L]]`;
#'   (2) an interaction of one or more numeric vectors and factors,
#'   as in `x` and `x:y`; and
#'   (3) a call to `|` with 1 or an interaction of one or more
#'   numeric variables on the left hand side and an interaction
#'   of one or more factors on the right hand side, as in
#'   `1 | f` and `x | f:g`.
#' @param frame
#'   A data frame containing the variables named in `x`,
#'   which must not have missing values, as enforced by
#'   [make_frame()].
#' @param sparse_X
#'   A logical scalar. If `TRUE`, then the matrix is returned in
#'   sparse format. If `x` is a call to `|`, then the matrix is
#'   returned in sparse format regardless of `sparse`.
#'
#' @return
#' A matrix with `nrow(frame)` rows and `n` columns,
#' where `n` is the number of nonempty groups if `x`
#' has a grouping component and 1 otherwise.
#'
#' @noRd
#' @importFrom Matrix sparseMatrix sparse.model.matrix
#' @importFrom stats model.matrix
term_to_matrix <- function(x, frame, sparse_X) {
  if (is.call(x) && x[[1L]] == as.name("|")) {
    group <- interaction(frame[all.vars(x[[3L]])], drop = TRUE)
    Z <- sparse.model.matrix(~-1 + group, data = list(group = group))
    if (x[[2L]] != 1) {
      Z <- Z * Reduce(`*`, frame[all.vars(x[[2L]])])
    }
    attr(Z, "assign") <- attr(Z, "contrasts") <- NULL
    attr(Z, "variables") <- list(numeric = all.vars(x[[2L]]), factor = all.vars(x[[3L]]))
    attr(Z, "group") <- group
    dimnames(Z) <- list(NULL, NULL)
    return(Z)
  }

  av <- all.vars(x)
  av_is_numeric <- vapply(frame[av], is.numeric, FALSE)
  if (all(av_is_numeric)) {
    X <- if (x == 1) rep.int(1L, nrow(frame)) else Reduce(`*`, frame[av])
    if (sparse_X) {
      X <- sparseMatrix(i = seq_along(X), j = rep.int(1L, length(X)),
                        x = X, dims = c(length(X), 1L))
    } else {
      dim(X) <- c(length(X), 1L)
    }
    attr(X, "variables") <- list(numeric = av, factor = character(0L))
    attr(X, "group") <- rep(factor(1), nrow(X))
    return(X)
  }

  mm <- if (sparse_X) sparse.model.matrix else model.matrix
  f <- if (any(av_is_numeric)) ~-1 + group else ~group
  group <- interaction(frame[av[!av_is_numeric]], drop = TRUE)
  X <- mm(f, data = list(group = group), contrasts.arg = list(group = "contr.sum"))
  if (any(av_is_numeric)) {
    X <- X * Reduce(`*`, frame[av[av_is_numeric]])
  }
  attr(X, "assign") <- attr(X, "contrasts") <- NULL
  attr(X, "variables") <- list(numeric = av[av_is_numeric], factor = av[!av_is_numeric])
  attr(X, "group") <- group
  dimnames(X) <- list(NULL, NULL)
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
#'   are such that `x[i]` is the number of cases observed from the
#'   end of day `t[i-1]` to the end of day `t[i]`.
#' }
#' \item{`wl`}{
#'   An integer vector of length `nlevels(index)` giving the length
#'   of each fitting window.
#' }
#' \item{`dow0`}{
#'   An integer vector giving the day-of-week of the earliest date
#'   in each fitting window, with Sunday mapping to 0, Monday mapping
#'   to 1, and so on.
#' }
#' \item{`Xs`, `Xd`}{
#'   The fixed effects design matrix in sparse and dense format.
#'   If `sparse_X = TRUE`, then `Xd` is empty. If `sparse_X = FALSE`,
#'   then `Xs` is empty. Uses zero-sum contrasts.
#' }
#' \item{`Z`}{
#'   The random effects design matrix in sparse format.
#'   Does not use contrasts.
#' }
#' \item{`fvar`, `rvar`}{
#'   Character vectors naming the independent variable
#'   in each fixed or random effects term.
#' }
#' \item{`fgrp`, `rgrp`}{
#'   Data frames containing the grouping variable corresponding
#'   to each fixed or random effects term.
#' }
#' \item{`fnc`, `rnc`}{
#'   Integer vectors listing the number of coefficients corresponding
#'   to each fixed or random effects term.
#' }
#' \item{`fid`, `rid`}{
#'   Indicator matrices with term labels for row names and parameter
#'   labels (with prefixes) for column names. Element `[i, j]` is equal
#'   to 1 if and only if term `i` appears in the mixed effects model
#'   (fixed/random component) for parameter `j`.
#' }
#' \item{`curve_flag`, `distr_flag`, `excess_flag`, `weekday_flag`, `sparse_X_flag`}{
#'   Integer flags referencing `curve`, `distr`, `excess`, `weekday`,
#'   and `sparse_X`.
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
#'   frame = frame,
#'   curve = curve,
#'   distr = distr,
#'   excess = excess,
#'   weekday = weekday,
#'   sparse_X = FALSE
#' )
#'
#' @keywords internal
#' @importFrom Matrix sparseMatrix
make_tmb_data <- function(frame, curve, distr, excess, weekday, sparse_X) {
  a <- attributes(frame)

  ## Parameter names (with prefixes)
  pn <- get_par_names(curve = curve, distr = distr,
                      excess = excess, weekday = TRUE, link = TRUE) # model
  pn0 <- get_par_names(link = TRUE) # all
  p <- length(pn)

  ## Compute time as number of days since start of fitting window.
  ## No issues with permutation as `window` is ordered by level.
  date_to_integer <- function(x) as.integer(julian(x, origin = x[1L]))
  t <- unlist(tapply(frame[[1L]], a$window, date_to_integer, simplify = FALSE), use.names = FALSE)

  ## Set unused first observations to NA. They are unused regardless,
  ## and this makes that explicit.
  x <- unlist(tapply(frame[[2L]], a$window, `[<-`, 1L, NA, simplify = FALSE), use.names = FALSE)

  ## Get lengths of fitting windows
  wl <- tabulate(a$window)

  ## Get day-of-week of first date, mapping Sunday to 0, etc.
  d0 <- do.call(c, tapply(frame[[1L]], a$window, `[`, 1L, simplify = FALSE))
  dow0 <- as.integer(julian(d0, origin = as.Date("1995-03-05"))) %% 7L # "1995-03-05" was a Sunday

  ## Get unique mixed effects terms
  ftl <- lapply(a$formula_glmm, attr, "fixed_terms")
  rtl <- lapply(a$formula_glmm, attr, "random_terms")
  ftl <- Map(`names<-`, ftl, lapply(ftl, function(l) vapply(l, deparse, "")))
  rtl <- Map(`names<-`, rtl, lapply(rtl, function(l) vapply(l, deparse, "")))
  ft_incl_dupl <- do.call(c, unname(ftl))
  rt_incl_dupl <- do.call(c, unname(rtl))
  ft <- ft_incl_dupl[!duplicated(names(ft_incl_dupl))]
  rt <- rt_incl_dupl[!duplicated(names(rt_incl_dupl))]

  ## Get design matrix for each term
  fmat <- lapply(ft, term_to_matrix, frame = frame, sparse_X = sparse_X)
  rmat <- lapply(rt, term_to_matrix, frame = frame)

  ## Get independent variable name for each term
  get_varname <- function(m) {
    s <- attr(m, "variables")$numeric
    if (length(s) == 0L) {
      return("(Intercept)")
    }
    paste(s, collapse = ":")
  }
  fvar <- lapply(fmat, get_varname)
  rvar <- lapply(rmat, get_varname)

  ## Get grouping variable for each term
  fgrp <- as.data.frame(lapply(fmat, attr, "group"))
  rgrp <- as.data.frame(lapply(rmat, attr, "group"))

  ## Get number of coefficients for each term
  fnc <- vapply(fmat, ncol, 0L)
  rnc <- vapply(rmat, ncol, 0L)

  ## Combine design matrices columnwise
  if (sparse_X) {
    Xs <- do.call(cbind, fmat)
    Xd <- matrix(integer(0L), ncol = ncol(Xs))
  } else {
    Xd <- do.call(cbind, fmat)
    Xs <- sparseMatrix(i = integer(0L), j = integer(0L),
                       x = integer(0L), dims = c(0L, ncol(Xd)))
  }
  if (length(rt) > 0L) {
    Z <- do.call(cbind, rmat)
  } else {
    Z <- sparseMatrix(i = integer(0L), j = integer(0L),
                      x = integer(0L), dims = c(nrow(frame), 0L))
  }

  ## Construct indicator matrices such that 1 at index (i,j)
  ## indicates that parameter j follows model using term i
  fid <- matrix(0L, nrow = length(ft), ncol = p, dimnames = list(names(ft), pn))
  for (s1 in rownames(fid)) { # FIXME: outer(FUN = `%in%`)?
    for (s2 in colnames(fid)) {
      fid[s1, s2] <- 1L * (s1 %in% ftl[[s2]])
    }
  }
  if (length(rt) > 0L) {
    rid <- matrix(0L, nrow = length(rt), ncol = p, dimnames = list(names(rt), pn))
    for (s1 in rownames(rid)) { # FIXME: outer(FUN = `%in%`)?
      for (s2 in colnames(rid)) {
        rid[s1, s2] <- 1L * (s1 %in% rtl[[s2]])
      }
    }
  } else {
    rid <- matrix(integer(0L), ncol = p, dimnames = list(NULL, pn))
  }

  l1 <- list(
    window = a$window,
    t = t,
    x = x,
    wl = wl,
    dow0 = dow0,
    Xs = Xs,
    Xd = Xd,
    Z = Z,
    fvar = fvar,
    rvar = rvar,
    fgrp = fgrp,
    rgrp = rgrp,
    fnc = fnc,
    rnc = rnc,
    fid = fid,
    rid = rid,
    curve_flag = get_flag("curve", curve),
    distr_flag = get_flag("distr", distr),
    excess_flag = 1L * excess,
    weekday_flag = 1L * weekday,
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
#' @param frame
#'   A model frame returned by [make_frame()].
#' @inheritParams egf
#'
#' @return
#' A list with elements `data`, `parameters`, `map`, `random`,
#' `DLL`, and `silent`.
#'
#' @examples
#' example("make_frame", package = "epigrowthfit")
#' tmb_args <- epigrowthfit:::make_tmb_args(
#'   frame = frame,
#'   curve = curve,
#'   distr = distr,
#'   excess = excess,
#'   weekday = weekday,
#'   sparse_X = FALSE,
#'   par_init = NULL
#' )
#' tmb_out <- do.call(TMB::MakeADFun, tmb_args)
#'
#' @seealso [make_tmb_data()], [make_tmb_parameters()]
#' @keywords internal
make_tmb_args <- function(frame, curve, distr, excess, weekday,
                          sparse_X, par_init) {
  data <- make_tmb_data(
    frame = frame,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    sparse_X = sparse_X
  )
  parameters <- make_tmb_parameters(
    tmb_data = data,
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
