deparen <- function(x) {
  while (is.call(x) && x[[1L]] == as.name("(")) {
    x <- x[[2L]]
  }
  x
}

do_minus <- function(x) {
  if (is.call(x) && x[[1L]] == as.name("-") && length(x) == 2L) {
    return(x[[2L]])
  }
  call("-", x)
}

is_interaction <- function(x) {
  if (is.call(x)) {
    if (x[[1L]] == as.name("(")) {
      return(is_interaction(x[[2L]]))
    }
    if (x[[1L]] == as.name(":")) {
      return(is_interaction(x[[2L]]) && is_interaction(x[[3L]]))
    }
  }
  if (is.name(x)) {
    return(TRUE)
  }
  FALSE
}

split_terms <- function(x) {
  if (inherits(x, "formula")) {
    x <- x[[length(x)]]
  }
  if (is.call(x)) {
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
    return(list(x))
  }
  if (is.name(x) || is.numeric(x)) {
    return(list(x))
  }
  stop("`x` must be a call, name, or number.")
}

unsplit_terms <- function(x) {
  is_pm_call <- function(x) {
    is.call(x) && (x[[1L]] == as.name("+") || x[[1L]] == as.name("-")) && length(x) == 2L
  }
  y <- x[[1L]]
  for (i in seq_along(x)[-1L]) {
    if (is_pm_call(x[[i]])) {
      y <- as.call(list(x[[i]][[1L]], y, x[[i]][[2L]]))
    } else {
      y <- call("+", y, x[[i]])
    }
  }
  y
}

#' @importFrom stats terms reformulate
expand_terms <- function(x) {
  fml <- as.formula(call("~", x))
  tl <- labels(terms(fml))
  if (length(tl) == 0L) {
    return(1)
  }
  reformulate(tl)[[2L]]
}
