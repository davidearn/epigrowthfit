check_formula_glmm <- function(formula_glmm, curve, distr, excess) {
  pn <- get_par_names(curve, distr, excess, link = TRUE)
  p <- length(pn)

  if (inherits(formula_glmm, "formula") && length(formula_glmm) == 2L) {
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
    formula_glmm <- formula_glmm[pn]
  } else {
    stop("`formula_glmm` must be a formula of the form `~terms`\n",
         "or a list of formulae of the form `par ~ terms`.")
  }




  lapply(formula_glmm, function(x) {
    ftrt <- split_effects(x)
    x[[3L]] <- unsplit_terms(do.call(c, ftrt))
    structure(x,
      fixed_terms = ftrt$fixed_terms,
      random_terms = ftrt$random_terms
    )
  })
}

split_effects <- function(x) {
  terms <- split_terms(expand_terms(x))
  has_bar <- vapply(terms, function(x) is.call(x) && x[[1L]] == as.name("|"), FALSE)
  fixed_terms <- terms[!has_bar]
  stop_if_not(
    vapply(fixed_terms, is_interaction, FALSE, strict = FALSE),
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
