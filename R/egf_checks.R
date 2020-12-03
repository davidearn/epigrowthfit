check_formula <- function(formula) {
  check(formula,
    what = "formula",
    len = 3L,
    yes = function(x) grepl("^[[:alpha:].]{1}[[:alnum:]._]* ~ [[:alpha:].]{1}[[:alnum:]._]*$", deparse(x)),
    "`formula` must be a formula of the form `y ~ x`."
  )
  formula
}

check_fixed <- function(fixed, par_names) {
  p <- length(par_names)
  if (is.null(fixed)) {
    fixed <- rep(list(as.formula("~1", env = .GlobalEnv)), p)
    names(fixed) <- par_names
    return(fixed)
  }
  if (inherits(fixed, "formula")) {
    fixed <- rep(list(fixed), p)
    names(fixed) <- par_names
  }
  check(fixed,
    what = "list",
    len = c(1L, NA),
    yes = function(x) all(sapply(x, inherits, "formula")),
    "`fixed` must be `NULL`, a formula, or a named list of formula."
  )
  check(names(fixed),
    what = "character",
    opt = par_names,
    no = function(x) any(duplicated(x)),
    "If `fixed` is a list, then `names(fixed)` must be a subset\n",
    "of `get_par_names(curve, distr, include_baseline)`."
  )
  check(sapply(fixed, length),
    val = 2L,
    "`fixed` formula must have the form `~rhs`."
  )
  rhs <- sapply(fixed, function(x) deparse(x[[2]]))
  check(rhs,
    yes = function(x) all(grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", x)),
    "`fixed` formula must be `~1` or have the form `~f1:...:fn`."
  )

  ## Fill out and order the list
  fixed[setdiff(par_names, names(fixed))] <- list(as.formula("~1", env = .GlobalEnv))
  fixed[par_names]
}

check_random <- function(random, par_names) {
  p <- length(par_names)
  if (is.null(random)) {
    random <- rep(list(NULL), p)
    names(random) <- par_names
    return(random)
  }
  if (inherits(random, "formula")) {
    random <- rep(list(random), p)
    names(random) <- par_names
  }
  check(random,
    what = "list",
    len = c(1L, NA),
    yes = function(x) all(sapply(x, inherits, "formula")),
    "`random` must be `NULL`, a formula, or a named list of formula."
  )
  check(names(random),
    what = "character",
    opt = par_names,
    no = function(x) any(duplicated(x)),
    "If `random` is a list, then `names(random)` must be a subset\n",
    "of `get_par_names(curve, distr, include_baseline)`."
  )
  check(sapply(random, length),
    val = 2L,
    "`random` formula must have the form `~rhs`."
  )
  rhs <- sapply(random, function(x) deparse(x[[2]]))
  check(rhs,
    yes = function(x) all(grepl("^(( \\+|- )?\\(1 \\| [[:alnum:]._]+((:[[:alnum:]._]+)*|(/[[:alnum:]._]+)*|(( \\* )[[:alnum:]._]+)*)\\))+$", x)),
    "`random` formula must have the form `~rhs`,\n",
    "with `rhs` a sum of one or more terms of the form\n",
    "`(1 | r1:...:rk)`, `(1 | r1/.../rm)`, or `(1 | r1 * ... * rn)`."
  )

  ## Fill out and order the list
  random[setdiff(par_names, names(random))] <- list(NULL)
  random[par_names]
}

check_data <- function(formula, data, index,
                       fixed, random,
                       dfmt, na_action) {
  ## Check that formula variables can be found
  check(data,
    what = c("data.frame", "list", "environment"),
    "`data` must be a data frame, list, or environment."
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
    stop(call. = FALSE,
      "Variables not found in `data`:\n",
      paste(an[!found], collapse = ", ")
    )
  }

  ## Check that time series is specified correctly
  if (is.character(data[[dn]])) {
    data[[dn]] <- try(as.Date(data[[dn]], tryFormats = dfmt), silent = TRUE)
  }
  check(data[[dn]],
    what = "Date",
    sprintf("`%s` must be of class \"Date\" or so coercible with\n`as.Date(%s, tryFormats = dfmt)`.", dn, dn)
  )
  check(data[[dn]],
    len = c(1L, NA),
    sprintf("`%s` must have nonzero length.", dn)
  )
  check(data[[dn]],
    no = anyNA,
    sprintf("`%s` must not have missing values.", dn)
  )
  n <- length(data[[dn]])
  check(data[[cn]],
    what = "numeric",
    len = n,
    sprintf("`%s` must be numeric and have length `length(%s)`.", cn, dn)
  )
  check(data[[cn]],
    val = c(0, Inf),
    rel = c(">=", "<"),
    yes = function(x) all(x %% 1 == 0, na.rm = TRUE),
    sprintf("Elements of `%s` must be non-negative integers or NA.", cn)
  )
  data[[cn]] <- as.integer(data[[cn]])

  ## Check that grouping variables are specified correctly
  if (length(fn) > 0L) {
    check(data[fn],
      yes = function(x) all(sapply(x, is.factor) & sapply(x, length) == n),
      sprintf("Grouping variables must be factors of length `length(%s)`", dn)
    )
    data[fn] <- lapply(data[fn], droplevels, exclude = NA)
  }

  ## Check that fitting windows are specified correctly
  if (is.null(index)) {
    index <- factor(integer(n))
  } else {
    check(index,
      what = "factor",
      len = n,
      sprintf("`index` must be a factor of length `length(%s)`.", dn)
    )
    index <- droplevels(index, exclude = NA)
    check(levels(index),
      len = c(1L, NA),
      "`index` must have at least one nonempty level."
    )
    check(na.omit(unclass(index)),
      yes = function(x) sum(diff(x) != 0L) == nlevels(index) - 1L,
      "Elements of `index` of a given level must be contiguous."
    )
  }
  data <- data.frame(data[an])[!is.na(index), ]
  index <- index[!is.na(index)]
  data_split <- split(data, index)
  date_split <- lapply(data_split, "[[", dn)
  check(date_split,
    no = function(x) any(sapply(x, length) < 7L),
    sprintf("`%s` must have length 7 or greater in each level of `index`.", dn)
  )
  check(date_split,
    yes = function(x) all(unlist(lapply(x, diff)) > 0),
    sprintf("`%s` must be increasing in each level of `index`.", dn)
  )
  cases_split <- lapply(data_split, "[[", cn)
  cases_split <- lapply(cases_split, "[", -1L)
  if (na_action == "fail") {
    check(cases_split,
      no = function(x) anyNA(unlist(x)),
      sprintf("There are fitting windows with missing values in `%s`.", cn)
    )
  }
  check(cases_split,
    no = function(x) any(sapply(x, function(y) sum(!is.na(y)) < 6L)),
    sprintf("There are fitting windows with insufficient data in `%s`.", cn)
  )
  factors_split <- lapply(data_split, "[", fn)
  factors_split <- lapply(factors_split, droplevels, exclude = NA)
  check(do.call(cbind, factors_split),
    yes = function(x) all(sapply(x, nlevels) == 1L),
    no = anyNA,
    "Factors must be constant and without missing values\n",
    "in each level of `index`."
  )

  structure(data,
    index = index,
    date_name = dn,
    cases_name = cn,
    fixed_term_labels = lapply(fixed, get_term_labels),
    random_term_labels = lapply(random, get_term_labels)
  )
}
