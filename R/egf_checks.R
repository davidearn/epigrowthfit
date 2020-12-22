check_formula <- function(formula) {
  stop_if_not(
    inherits(formula, "formula"),
    length(formula) == 3L,
    grepl("^[[:alpha:].]{1}[[:alnum:]._]* ~ [[:alpha:].]{1}[[:alnum:]._]*$", deparse(formula)),
    m = "`formula` must be a formula of the form `y ~ x`."
  )
  formula
}

check_fixed <- function(fixed, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (is.null(fixed)) {
    fixed <- rep(list(~1), p)
    names(fixed) <- pn
    return(fixed)
  }
  if (inherits(fixed, "formula")) {
    fixed <- rep(list(fixed), p)
    names(fixed) <- pn
  }
  stop_if_not(
    inherits(fixed, "list"),
    length(fixed) >= 1L,
    vapply(fixed, inherits, logical(1L), "formula"),
    m = "`fixed` must be NULL, a formula, or a named list of formula."
  )
  names(fixed) <- add_link_string(remove_link_string(names(fixed)))
  stop_if_not(
    !is.null(names(fixed)),
    names(fixed) %in% pn,
    !any(duplicated(names(fixed))),
    m = paste0(
      "If `fixed` is a list, then `names(fixed)` must be a subset\n",
      "of `get_par_names(curve, distr, excess, link)`."
    )
  )
  rhs <- vapply(fixed, function(x) deparse(x[[2L]]), character(1L))
  stop_if_not(
    lengths(fixed) == 2L,
    grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", rhs),
    m = "`fixed` formula must be `~1` or have the form `~f1:...:fn`."
  )

  ## Fill out and order the list
  fixed[setdiff(pn, names(fixed))] <- list(~1)
  fixed[pn]
}

check_random <- function(random, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (is.null(random)) {
    random <- rep(list(NULL), p)
    names(random) <- pn
    return(random)
  }
  if (inherits(random, "formula")) {
    random <- rep(list(random), p)
    names(random) <- pn
  }
  stop_if_not(
    inherits(random, "list"),
    length(random) >= 1L,
    vapply(random, inherits, logical(1L), "formula"),
    m = "`random` must be NULL, a formula, or a named list of formula."
  )
  names(random) <- add_link_string(remove_link_string(names(random)))
  stop_if_not(
    !is.null(names(random)),
    names(random) %in% pn,
    !any(duplicated(names(random))),
    m = paste0(
      "If `random` is a list, then `names(random)` must be a subset\n",
      "of `get_par_names(curve, distr, excess, link)`."
    )
  )
  rhs <- vapply(random, function(x) deparse(x[[2L]]), character(1L))
  stop_if_not(
    lengths(random) == 2L,
    grepl("^(( \\+|- )?\\(1 \\| [[:alnum:]._]+((:[[:alnum:]._]+)*|(/[[:alnum:]._]+)*|(( \\* )[[:alnum:]._]+)*)\\))+$", rhs),
    m = paste0(
      "`random` formula must have the form `~rhs`,\n",
      "with `rhs` a sum of one or more terms of the form\n",
      "`(1 | r1:...:rk)`, `(1 | r1/.../rm)`, or `(1 | r1 * ... * rn)`."
    )
  )

  ## Fill out and order the list
  random[setdiff(pn, names(random))] <- list(NULL)
  random[pn]
}

check_data <- function(formula, fixed, random,
                       data, index,
                       date_format, na_action) {
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
      m = sprintf("Grouping variables must be factors of length `length(%s)`", dn)
    )
    data[fn] <- lapply(data[fn], droplevels, exclude = NA)
  }
  data <- data.frame(data[an])

  ## Check that fitting windows are specified correctly
  if (is.null(index)) {
    index <- factor(integer(n))
  } else {
    stop_if_not(
      is.atomic(index),
      length(index) == n,
      m = sprintf("`index` must be a factor or atomic vector of length `length(%s)`.", dn)
    )
    index <- droplevels(index, exclude = NA)
    stop_if_not(
      nlevels(index) > 0L,
      m = "`index` must have at least one nonempty level."
    )
    not_na <- !is.na(index)
    index <- index[not_na]
    data <- data[not_na, ]
    ord <- order(index)
    index <- index[ord]
    data <- data[ord, ]
  }
  data_split <- split(data, index)
  date_split <- lapply(data_split, "[[", dn)
  stop_if_not(
    lengths(date_split) >= 7L,
    m = sprintf("`%s` must have length 7 or greater in each level of `index`.", dn)
  )
  stop_if_not(
    vapply(date_split, function(x) all(diff(x) > 0), logical(1L)),
    m = sprintf("`%s` must be increasing in each level of `index`.", dn)
  )
  cases_split <- lapply(data_split, "[[", cn)
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
  factors_split <- lapply(data_split, "[", fn)
  factors_split <- lapply(factors_split, droplevels, exclude = NA)
  stop_if_not(
    vapply(factors_split, function(x) all(vapply(x, nlevels, integer(1L)) == 1L) && !anyNA(x), logical(1L)),
    m = paste0(
      "Grouping variables must be constant and without\n",
      "missing values in each level of `index`."
    )
  )
  row.names(data) <- NULL

  structure(data,
    index = index,
    date_name = dn,
    cases_name = cn,
    fixed_term_labels = lapply(fixed, get_term_labels),
    random_term_labels = lapply(random, get_term_labels)
  )
}
