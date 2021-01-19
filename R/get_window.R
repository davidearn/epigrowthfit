get_window <- function(formula, data, group_by = ~1, spar = 0.66) {
  ## Check that `data` makes sense
  stop_if_not(
    inherits(data, c("data.frame", "list", "environment")),
    # typeof(data) %in% c("list", "environment"), # FIXME: better/worse?
    m = "`data` must be a data frame, list, or environment."
  )
  if (is.environment(data)) {
    data <- as.list(data)
  }

  ## Check that formulae make sense
  check_formula(formula)
  stop_if_not(
    inherits(group_by, "formula"),
    length(group_by) == 2L,
    grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", deparse(group_by[[2L]])),
    all.vars(group_by) %in% names(data),
    m = paste0(
      "`group_by` must be a formula of the form\n",
      "`~1` or `~f1:...:fn`, with `f1`,...,`fn`\n",
      "naming factors in `data`."
    )
  )

  ## Check that formula variables can be found
  dn <- all.vars(formula[[3L]]) # date
  cn <- all.vars(formula[[2L]]) # cases
  fn <- all.vars(group_by) # factors
  an <- unique(c(dn, cn, fn)) # all
  found <- an %in% names(data)
  stop_if_not(
    all(found),
    m = paste0(
      "Variables not found in `data`:\n",
      paste(an[!found], collapse = ", ")
    )
  )

  ## Check that time series is specified correctly
  if (is.character(data[[dn]])) {
    data[[dn]] <- try(as.Date(data[[dn]], tryFormats = date_format), silent = TRUE)
  }
  stop_if_not(
    inherits(data[[dn]], "Date"),
    m = sprintf("`%s` must be of class \"Date\" or so coercible with\n`as.Date(%s, tryFormats = date_format)`.", dn, dn)
  )
  n <- length(data[[dn]])
  stop_if_not(
    n > 0L,
    m = sprintf("`%s` must have nonzero length.", dn)
  )
  stop_if_not(
    !anyNA(data[[dn]]),
    m = sprintf("`%s` must not have missing values.", dn)
  )
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
      vapply(data[fn], is.factor, FALSE),
      lengths(data[fn]) == n,
      m = sprintf("Grouping variables must be factors of length `length(%s)`.", dn)
    )
    data[fn] <- lapply(data[fn], droplevels, exclude = NA)
    stop_if_not(
      !vapply(data[fn], anyNA, FALSE),
      m = "Grouping variables must not have missing values."
    )
  }

  group <- get_factor(get_term_labels(group_by), as.data.frame(data[fn]))
  ts <- as.data.frame(data[c(dn, cn)])
  ts_split <- split(ts, group)

  date_split <- lapply(ts_split, "[[", dn)
  stop_if_not(
    lengths(date_split) >= 7L,
    m = sprintf("`%s` must have length 7 or greater in each group.", dn)
  )
  stop_if_not(
    vapply(date_split, function(x) all(diff(x) > 0), FALSE),
    m = sprintf("`%s` must be increasing in each group.", dn)
  )
  cases_split <- lapply(frame_split, "[[", cn)
  stop_if_not(
    vapply(cases_split, function(x) sum(!is.na(x)) >= 7L, FALSE),
    m = sprintf("There are groups with insufficient data\n(fewer than 7 observations) in `%s`.", cn)
  )



}
