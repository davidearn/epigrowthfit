#' Nonstandard evaluation
#'
#' Utilities for obtaining index vectors from unevaluated \code{subset},
#' \code{append}, and \code{order} expressions and character vectors from
#' unevaluated \code{label} expressions.
#'
#' @param expr
#'   A \link{language} object to be evaluated in \code{data} or,
#'   in the case of \code{egf_eval_append},
#'   in \code{`\link{names<-}`(\link{as.list}(\link{seq_along}(data)), \link{names}(data))}.
#'   Alternatively, a non-language object to be used in place of
#'   the hypothetical result of evaluation.
#' @param data
#'   A \link[=data.frame]{data frame}.
#' @param enclos
#'   An \link{environment} to enclose \code{data}.
#'
#' @details
#' \code{subset} must evaluate to a valid index vector
#' for the rows of \code{data} (see \code{\link{[.data.frame}}).
#'
#' \code{order} must evaluate to a permutation
#' of \code{\link{seq_len}(\link{nrow}(data))}.
#'
#' \code{append} must evaluate to a valid index vector
#' for the variables in \code{data} (see \code{\link{[.data.frame}}).
#'
#' \code{label} must evaluate to an \R object coercible
#' via \code{\link{as.character}} to a \link{character}
#' vector of length 1 or length \code{\link{nrow}(data)}.
#'
#' @return
#' Let \code{res} be the result of evaluating \code{expr}.
#'
#' \code{egf_eval_subset} returns an increasing \link{integer} vector
#' indexing rows of \code{data}, obtained after some processing of
#' \code{res}.
#'
#' \code{egf_eval_order} returns \code{\link{as.integer}(res)}.
#'
#' \code{egf_eval_append} returns
#' \code{\link{match}(\link{names}(data[res]), \link{names}(data), 0L)},
#' with zeros (if any) deleted.
#'
#' \code{egf_eval_label} returns \code{\link{as.character}(res)},
#' repeated to length \code{\link{nrow}(data)}.
#'
#' @examples
#' ## year <- 2021L
#' ## data <- data.frame(
#' ##   month = sample(month.abb, 20L, replace = TRUE),
#' ##   day = sample(30L, 20L, replace = TRUE)
#' ## )
#' ##
#' ## subset <- quote(grepl("^J", month) & day < 16L)
#' ## egf_eval_subset(subset, data)
#' ##
#' ## order <- quote(order(month, day))
#' ## egf_eval_order(order, data)
#' ##
#' ## append <- quote(-day)
#' ## egf_eval_append(append, data)
#' ##
#' ## label <- quote(sprintf("%d-%02d-%02d", year, match(month, month.abb, 0L), day))
#' ## egf_eval_label(label, data)
#'
#' @name egf_eval
#' @keywords internal
NULL

#' @rdname egf_eval
egf_eval_subset <- function(expr, data, enclos = parent.frame()) {
  stopifnot(is.data.frame(data))
  n <- nrow(data)
  if (is.null(expr)) {
    return(seq_len(n))
  }
  if (is.language(expr)) {
    subset <- eval(expr, data, enclos)
  } else {
    subset <- expr
  }
  index <- seq_len(n)
  names(index) <- row.names(data)
  subset <- index[subset]
  sort(unique(subset[!is.na(subset)]))
}

#' @rdname egf_eval
egf_eval_order <- function(expr, data, enclos = parent.frame()) {
  stopifnot(is.data.frame(data))
  n <- nrow(data)
  if (is.null(expr)) {
    return(seq_len(n))
  }
  if (is.language(expr)) {
    order <- eval(expr, data, enclos)
  } else {
    order <- expr
  }
  stopifnot(
    is.numeric(order),
    length(order) == n,
    sort(order) == seq_len(n)
  )
  as.integer(order)
}

#' @rdname egf_eval
egf_eval_append <- function(expr, data, enclos = baseenv()) {
  stopifnot(is.data.frame(data))
  if (is.null(expr)) {
    return(integer(0L))
  }
  if (is.language(expr)) {
    l <- as.list(seq_along(data))
    names(l) <- names(data)
    append <- eval(expr, l, enclos)
  } else {
    append <- expr
  }
  append <- match(names(data[append]), names(data), 0L)
  append[append > 0L]
}

#' @rdname egf_eval
egf_eval_label <- function(expr, data, enclos = parent.frame()) {
  stopifnot(is.data.frame(data))
  if (is.null(expr)) {
    return(NULL)
  }
  if (is.language(expr)) {
    label <- eval(expr, data, enclos)
  } else {
    label <- expr
  }
  label <- as.character(label)
  stopifnot(length(label) %in% c(1L, n <- nrow(data)))
  rep_len(label, n)
}
