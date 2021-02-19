check_formula_glmm <- function(formula_glmm, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (inherits(formula_glmm, "formula") &&
      length(formula_glmm) == 2L) {
    formula_glmm <- check_formula_glmm_terms(formula_glmm)
    formula_glmm[[3L]] <- deparen(formula_glmm[[2L]])
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
        "`get_par_names(curve, distr, excess, link = TRUE)`."
      )
    )
    names(formula_glmm) <- lhs_as_character
    make_default_formula <- function(s) {
      as.formula(call("~", as.name(s), 1))
    }
    ss <- setdiff(pn, lhs_as_character)
    formula_glmm[ss] <- lapply(ss, make_default_formula)
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
    vapply(terms[!has_bar], is_interaction, FALSE, strict = FALSE),
    m = paste0(
      "Fixed effects formulae must expand\n",
      "to a sum of names and interactions\n",
      "of names and calls to I()."
    )
  )
  if (any(has_bar)) {
    lhs <- lapply(terms[has_bar], `[[`, 2L)
    stop_if_not(
      vapply(lhs, `==`, FALSE, 1) | vapply(lhs, is_interaction, FALSE, strict = FALSE),
      m = paste0(
        "Left hand side of `|` in `formula_glmm`\n",
        "must expand to 1 or a sum of names and\n",
        "interactions of names and calls to I()."
      )
    )
    rhs <- lapply(terms[has_bar], `[[`, 3L)
    stop_if_not(
      vapply(rhs, is_interaction, FALSE, strict = TRUE),
      m = paste0(
        "Right hand side of `|` in `formula_glmm`\n",
        "must expand to a sum of names and interactions of names."
      )
    )
  }

  structure(x,
    fixed_terms = terms[!has_bar],
    random_terms = terms[has_bar]
  )
}
