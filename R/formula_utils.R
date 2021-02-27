#' Negate formula terms
#'
#' Negates formula terms, being careful to resolve `-(-x)` to `x`.
#'
#' @param x A call, name, or atomic scalar.
#'
#' @return
#' A call, name, or atomic scalar.
#'
#' @examples
#' x <- as.name(x)
#' minus_x <- call("-", x)
#' do_minus(x)
#' do_minus(minus_x)
#'
#' @noRd
do_minus <- function(x) {
  if (is.call(x) && x[[1L]] == as.name("-") && length(x) == 2L) {
    return(x[[2L]])
  }
  call("-", x)
}

#' Split and unsplit formula terms
#'
#' Split formula terms related by ```+``` and ```-```
#' operators with `split_terms()`. Perform the inverse
#' operation with `unsplit_terms()`.
#'
#' @param x A call, name, or atomic scalar.
#' @param l A list of calls, names, and atomic scalars.
#'
#' @details
#' If `x` is a formula, then `x[[length(x)]]` is split.
#' Hence `unsplit_terms(split_terms(x))` reproduces
#' `x[[length(x)]]` but not `x`.
#'
#' Extraneous `(` calls, as in the expression
#' `(x) + (y)` are stripped when splitting and
#' not replaced when unsplitting. In these cases,
#' `identical(x, unsplit_terms(split_terms(x)))`
#' may return `FALSE`.
#'
#' In contrast to the behaviour of [stats::terms()],
#' zeros (0) and ones (1) are not treated specially,
#' and expressions like `x * y - y` are not expanded
#' and simplified prior to splitting.
#'
#' @return
#' A list of calls, names, and atomic scalars.
#'
#' @examples
#' x <- str2lang("1 + a + b * c - c + (d | e) + (0 + f * g | h)")
#' l <- split_terms(x)
#' y <- unsplit_terms(l)
#' identical(x, y)
#' # [1] FALSE
#'
#' @name split_terms
NULL

#' @rdname split_terms
#' @export
split_terms <- function(x) {
  if (inherits(x, "formula")) {
    x <- x[[length(x)]]
  }
  if (is.name(x) || (is.atomic(x) && length(x) == 1L)) {
    return(list(x))
  }
  stop_if_not(
    is.call(x),
    m = "`x` must be a call, name, or atomic scalar."
  )
  if (x[[1L]] == as.name("(")) {
    return(split_terms(x[[2L]]))
  }
  if (x[[1L]] == as.name("+")) {
    if (length(x) == 2L) {
      return(split_terms(x[[2L]]))
    }
    return(c(split_terms(x[[2L]]), split_terms(x[[3L]])))
  }
  if (x[[1L]] == as.name("-")) {
    if (length(x) == 2L) {
      return(lapply(split_terms(x[[2L]]), do_minus))
    }
    return(c(split_terms(x[[2L]]), split_terms(x[-2L])))
  }
  list(x)
}

#' @rdname split_terms
#' @export
unsplit_terms <- function(l) {
  stop_if_not(
    inherits(l, "list"),
    m = "`l` must be a list."
  )
  if (length(l) == 0L) {
    return(NULL)
  }
  is_pm_call <- function(x) {
    is.call(x) && (x[[1L]] == as.name("+") || x[[1L]] == as.name("-")) && length(x) == 2L
  }
  x <- l[[1L]]
  for (i in seq_along(l)[-1L]) {
    if (is_pm_call(l[[i]])) {
      x <- as.call(list(l[[i]][[1L]], x, l[[i]][[2L]]))
    } else {
      x <- call("+", x, l[[i]])
    }
  }
  x
}

#' Is an object a call to ```|```?
#'
#' Test whether an object is a call to the ```|``` operator.
#'
#' @param x An \R object, typically a formula term.
#'
#' @return
#' `TRUE` or `FALSE`.
#'
#' @examples
#' x <- call("|", as.name("f"), as.name("g"))
#' x_plus_x <- call("+", x, x)
#' is_bar(x)
#' is_bar(x_plus_x)
#'
#' @noRd
is_bar <- function(x) {
  is.call(x) && x[[1L]] == as.name("|")
}

#' Expand and simplify formula terms
#'
#' Expand and simplify formulae and formula terms
#' (including those with calls to the ```|``` operator)
#' using `terms(simplify = TRUE)` machinery.
#'
#' @param x A call, name, or atomic scalar.
#'
#' @return
#' A call, name, or atomic scalar.
#'
#' @examples
#' expand_terms(~x * y)
#' expand_terms(~x * y + (1 | f/g))
#' expand_terms(~x * y + (1 | f/g) + (a | f) + (b | f:g))
#'
#' @export
#' @importFrom stats terms
expand_terms <- function(x) {
  if (inherits(x, "formula")) {
    x[-1L] <- lapply(x[-1L], expand_terms)
    return(x)
  }
  if (is.name(x) || (is.atomic(x) && length(x) == 1L)) {
    return(x)
  }
  stop_if_not(
    is.call(x),
    m = "`x` must be a call, name, or atomic scalar."
  )
  tl <- split_terms(x)
  tl_is_bar <- vapply(tl, is_bar, FALSE)
  if (all(tl_is_bar)) {
    no_bar <- NULL
  } else {
    no_bar <- terms(as.formula(call("~", unsplit_terms(tl[!tl_is_bar]))), simplify = TRUE)[[2L]]
    if (!any(tl_is_bar)) {
      return(no_bar)
    }
  }
  expand_bar <- function(x) {
    lhs <- expand_terms(x[[2L]])
    rhs <- split_terms(expand_terms(x[[3L]]))
    lapply(rhs, function(x) call("|", lhs, x))
  }
  merge_bars <- function(l) {
    lhs <- lapply(l, `[[`, 2L)
    rhs <- l[[1L]][[3L]]
    call("|", expand_terms(unsplit_terms(lhs)), rhs)
  }
  bar <- do.call(c, lapply(tl[tl_is_bar], expand_bar))
  rhs_deparsed <- vapply(bar, function(x) deparse(x[[3L]]), "")
  bar <- tapply(bar, rhs_deparsed, merge_bars, simplify = FALSE)
  unsplit_terms(c(no_bar, bar))
}

#' Replace ```|``` with ```+``` in formula terms
#'
#' Replaces the ```|``` operator in formulae and formula terms
#' with the ```+``` operator, enabling construction of model frames
#' from mixed effects formulae using [stats::model.frame()].
#'
#' @param x A call, name, or atomic scalar.
#'
#' @return
#' A call, name, or atomic scalar.
#'
#' @examples
#' gsub_bar_plus(~x + (1 | f))
#'
#' @export
gsub_bar_plus <- function(x) {
  if (inherits(x, "formula")) {
    x[-1L] <- lapply(x[-1L], gsub_bar_plus)
    return(x)
  }
  if (is.name(x) || (is.atomic(x) && length(x) == 1L)) {
    return(x)
  }
  stop_if_not(
    is.call(x),
    m = "`x` must be a call, name, or atomic scalar."
  )
  tl <- split_terms(x)
  tl_is_bar <- vapply(tl, is_bar, FALSE)
  tl[tl_is_bar] <- lapply(tl[tl_is_bar], `[[<-`, 1L, as.name("+"))
  unsplit_terms(tl)
}

split_effects <- function(x) {
  tl <- split_terms(x)
  tl_is_bar <- vapply(tl, is_bar, FALSE)
  fixed_rhs <- if (all(tl_is_bar)) 1 else unsplit_terms(tl[!tl_is_bar])
  list(
    fixed = as.formula(call("~", fixed_rhs), env = environment(x)),
    random = tl[tl_is_bar]
  )
}

split_interaction <- function(x) {
  if (is.call(x) && x[[1L]] == as.name(":")) {
    return(do.call(c, lapply(x[-1L], split_interaction)))
  }
  list(x)
}



#
# extract_level <- function(term_label, column_name) {
#   f <- function(s1, s1s2) {
#     ## Handle "(Intercept)" and numeric variable names
#     if (s1 == s1s2) {
#       return(s1)
#     }
#     ## Extract `factor_level` from `paste0(factor_name, factor_level)`
#     sub(sprintf("^%s", s1), "", s1s2)
#   }
#   split_string_at_colon <- function(s) {
#     ## `strsplit(split = ":")` is naive as there could be
#     ## colons in the variable name
#     vapply(split_interaction(str2lang(s)), deparse, "")
#   }
#   ss <- mapply(f,
#     s1 = split_string_at_colon(term_label),
#     s2 = split_string_at_colon(column_name)
#   )
#   paste(ss, collapse = ":")
# }
