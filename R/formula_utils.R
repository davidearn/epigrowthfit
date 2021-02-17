do_minus <- function(x) {
  if (is.call(x) && x[[1L]] == as.name("-") && length(x) == 2L) {
    return(x[[2L]])
  }
  call("-", x)
}

split_terms <- function(x) {
  if (inherits(x, "formula")) {
    x <- x[[length(x)]]
  }
  if (is.name(x) || is.numeric(x)) {
    return(list(x))
  }
  stop_if_not(
    is.call(x),
    m = "`x` must be a call, name, or number."
  )
  if (x[[1L]] == as.name("(")) {
    return(split_terms(x[[2L]]))
  }
  if (x[[1L]] == as.name("+")) {
    if (length(x) == 2L) {
      return(split_terms(x[[2L]]))
    } else {
      return(c(split_terms(x[[2L]]), split_terms(x[[3L]])))
    }
  }
  if (x[[1L]] == as.name("-")) {
    if (length(x) == 2L) {
      return(lapply(split_terms(x[[2L]]), do_minus))
    } else {
      return(c(split_terms(x[[2L]]), split_terms(x[-2L])))
    }
  }
  list(x)
}

unsplit_terms <- function(l) {
  stop_if_not(
    inherits(l, "list"),
    length(l) > 0L,
    m = "`l` must be a list of nonzero length."
  )
  stop_if_not(
    vapply(l, function(x) is.call(x) || is.name(x) || is.numeric(x), FALSE),
    m = "Elements of `l` must be calls, names, or numbers."
  )
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

#' @importFrom stats terms reformulate
expand_terms <- function(x) {
  if (inherits(x, "formula")) {
    x <- x[[length(x)]]
  }
  if (is.numeric(x) || is.name(x)) {
    return(x)
  }
  stop_if_not(
    is.call(x),
    m = "`x` must be a call, name, or number."
  )
  l <- split_terms(x)
  has_bar <- vapply(l, function(x) is.call(x) && x[[1L]] == as.name("|"), FALSE)
  if (any(has_bar)) {
    f <- function(bar) {
      lhs <- split_terms(expand_terms(x[[2L]]))
      rhs <- split_terms(expand_terms(x[[3L]]))
      Map(function(a, b) call("|", a, b),
        a = rep(lhs, each = length(rhs)),
        b = rep(rhs, times = length(lhs))
      )
    }
    l <- c(l[!has_bar], do.call(c, lapply(l[has_bar], f)))
  }
  reformulate(labels(terms(call("~", unsplit_terms(l)))))[[2L]]
}

deparen <- function(x) {
  if (is.name(x) || is.numeric(x)) {
    return(x)
  }
  stop_if_not(
    is.call(x),
    m = "`x` must be a call, name, or number."
  )
  ops <- c("~", "|", "+", "-", "^", "*", "/", ":", "%in%")
  if (deparse(x[[1L]]) %in% ops) {
    x[[2L]] <- deparen(x[[2L]])
    if (length(x) == 3L) {
      x[[3L]] <- deparen(x[[3L]])
    }
    return(x)
  }
  while (is.call(x) && x[[1L]] == as.name("(")) {
    x <- x[[2L]]
  }
  x
}

is_interaction <- function(x, strict = TRUE) {
  if (is.name(x)) {
    return(TRUE)
  }
  if (is.call(x)) {
    if (strict) {
      return(x[[1L]] == as.name(":") && is_interaction(x[[2L]], strict) && is.name(x[[3L]]))
    }
    if (x[[1L]] == as.name(":")) {
      return(is_interaction(x[[2L]], strict) && is_interaction(x[[3L]], strict))
    }
    if (x[[1L]] == as.name("(") {
      return(is_interaction(x[[2L]], strict))
    }
    if (x[[1L]] == as.name("I")) {
      return(TRUE)
    }
  }
  FALSE
}
