#' Nonstandard evaluation
#'
#' Utilities for obtaining index vectors from unevaluated \code{subset},
#' \code{append}, and \code{order} expressions and character vectors from
#' unevaluated \code{label} expressions.
#'
#' @param subset,order,label
#'   Expressions to be evaluated in \code{envir}.
#' @param append
#'   An expression to be evaluated in
#'   \code{`\link{names<-}`(\link{as.list}(\link{seq_along}(envir)), \link{names}(envir))}.
#' @param envir,enclos
#'   See \code{\link{eval}}. \code{\link{is.list}(envir)} must be \code{TRUE}
#'   for \code{eval_append}.
#' @param .subset,.order,.append,.label
#'   \link[=atomic]{Atomic} vectors to be used (if non-\code{\link{NULL}})
#'   in place of the result of evaluating the undotted argument.
#'
#' @details
#' \code{subset} must evaluate to a \link{logical} vector.
#' The result must have length \code{\link{nrow}(envir)}
#' if \code{envir} is a \link[=data.frame]{data frame}.
#'
#' \code{order} must evaluate to a permutation of \code{\link{seq_len}(n)}
#' for some \code{n}. \code{n} must be \code{nrow(envir)} if \code{envir}
#' is a data frame.
#'
#' \code{append} must evaluate to an \link{atomic} vector
#' indexing \code{envir}.
#'
#' \code{label} must evaluate to an \link{atomic} vector.
#' The result must have length of length 1 or length \code{nrow(envir)}
#' if \code{envir} is a data frame.
#'
#' @section Note:
#' \code{subset} and \code{append} are processed similarly to arguments
#' \code{subset} and \code{select} of function \code{\link{subset}}.
#' See \code{\link{subset}} for additional usage examples.
#'
#' @section Warning:
#' Nonstandard evaluation of \code{subset}, \code{order}, \code{append},
#' and \code{label} is intended only to make interactive use more convenient.
#' To avoid unexpected behaviour, especially when programming, use the dotted
#' arguments.
#'
#' @return
#' Let \code{r} be the result of evaluating the supplied expression.
#'
#' \code{eval_subset} returns \code{r},
#' with \code{\link{NA}} (if any) coerced to \code{FALSE}.
#'
#' \code{eval_order} returns \code{r} as is.
#'
#' \code{eval_append} returns
#' \code{\link{match}(\link{names}(envir[r]), \link{names}(envir), 0L)},
#' with zeros (if any) deleted.
#'
#' \code{eval_label} returns \code{\link{as.character}(r)},
#' repeated to length \code{\link{nrow}(envir)}
#' if \code{envir} is a \link[=data.frame]{data frame}.
#'
#' @examples
#' year <- 2021L
#' data <- data.frame(
#'   month = sample(month.abb, 20L, replace = TRUE),
#'   day = sample(30L, 20L, replace = TRUE)
#' )
#'
#' subset <- quote(grepl("^J", month) & day < 16L)
#' order <- quote(order(month, day))
#' append <- quote(-day)
#' label <- quote(sprintf("%d-%02d-%02d", year, match(month, month.abb, 0L), day))
#'
#' ## eval_subset(subset, data, parent.frame())
#' ## eval_order(order, data, parent.frame())
#' ## eval_append(append, data, parent.frame())
#' ## eval_label(label, data, parent.frame())
#'
#' @name nse
#' @keywords internal
NULL

#' @rdname nse
eval_subset <- function(subset, envir, enclos, .subset = NULL) {
  if (is.null(.subset)) {
    if (is.null(subset)) {
      return(if (is.data.frame(envir)) rep_len(TRUE, nrow(envir)) else NULL)
    }
    subset <- eval(subset, envir, enclos)
  } else {
    subset <- .subset
  }
  stopifnot(is.logical(subset))
  if (is.data.frame(envir)) {
    stopifnot(length(subset) == nrow(envir))
  }
  !is.na(subset) & subset
}

#' @rdname nse
eval_append <- function(append, envir, enclos, .append = NULL) {
  stopifnot(is.list(envir))
  if (is.null(.append)) {
    if (is.null(append)) {
      return(integer(0L))
    }
    l <- as.list(seq_along(envir))
    names(l) <- names(envir)
    append <- eval(append, l, enclos)
  } else {
    append <- .append
  }
  m <- match(names(envir[append]), names(envir), 0L)
  m[m > 0L]
}

#' @rdname nse
eval_order <- function(order, envir, enclos, .order = NULL) {
  if (is.null(.order)) {
    if (is.null(order)) {
      return(if (is.data.frame(envir)) seq_len(nrow(envir)) else NULL)
    }
    order <- eval(order, envir, enclos)
  } else {
    order <- .order
  }
  stopifnot(
    is.numeric(order),
    sort(order) == seq_along(order)
  )
  if (is.data.frame(envir)) {
    stopifnot(length(order) == nrow(envir))
  }
  order
}

#' @rdname nse
eval_label <- function(label, envir, enclos, .label = NULL) {
  if (is.null(.label)) {
    if (is.null(label)) {
      return(NULL)
    }
    label <- eval(label, envir, enclos)
  } else {
    label <- .label
  }
  stopifnot(is.atomic(label))
  label <- as.character(label)
  if (is.data.frame(envir)) {
    stopifnot(length(label) %in% c(1L, n <- nrow(envir)))
    return(rep_len(label, n))
  }
  label
}
