#' Nonstandard evaluation
#'
#' Utilities for obtaining index vectors from unevaluated \code{subset},
#' \code{append}, and \code{order} expressions and character vectors from
#' unevaluated \code{label} expressions (\code{xlab}, \code{ylab}, etc.).
#'
#' @param subset,order,label
#'   Expressions to be evaluated in \code{data}.
#' @param append
#'   An expression to be evaluated in
#'   \code{`\link{names}<-`(\link{as.list}(\link{seq_along}(data)), names(data))}.
#' @param data
#'   A \link[=data.frame]{data frame}.
#' @param enclos
#'   An \link{environment} to enclose \code{data}.
#' @param .subset,.order,.append,.label
#'   \link[=atomic]{Atomic} vectors to be used (if non-\code{\link{NULL}})
#'   in place of the result of evaluating the undotted argument.
#'
#' @details
#' \code{subset} must evaluate to a logical vector of length
#' \code{n = \link{nrow}(data)}.
#' \code{\link{NULL}} is equivalent to \code{\link{rep_len}(TRUE, n)}.
#'
#' \code{order} must evaluate to a permutation of \code{\link{seq_len}(n)}.
#' \code{\link{NULL}} is equivalent to \code{\link{seq_len}(n)}.
#'
#' \code{append} must evaluate to an \link{atomic} vector indexing \code{data}.
#' \code{\link{NULL}} is equivalent to \code{\link{integer}(0L)}.
#'
#' \code{label} must evaluate to an \link{atomic} vector of length 1
#' or \code{n}.
#' \code{\link{NULL}} is a no-op.
#'
#' @section Note:
#' Note that \code{subset} and \code{append} are processed similarly to
#' arguments \code{subset} and \code{select} of function \code{\link{subset}}.
#' See \code{\link{subset}} for additional usage examples.
#'
#' @section Warning:
#' Nonstandard evaluation of \code{subset}, \code{order}, \code{append},
#' and \code{label} is intended only to make interactive use more convenient.
#' To avoid unexpected behaviour, especially when programming, use the dotted
#' versions.
#'
#' @return
#' Let \code{r} be the result of evaluating the supplied expression or,
#' otherwise, the dotted argument.
#'
#' \code{subset_to_index} returns \code{r} with \code{FALSE} replacing
#' \code{\link{NA}} (if any).
#'
#' \code{order_to_index} returns \code{r} as is.
#'
#' \code{append_to_index} returns
#' \code{\link{match}(\link{names}(data[r]), names(data), 0L)},
#' with zeros (if any) deleted.
#'
#' \code{label_to_character} returns
#' \code{\link{rep_len}(\link{as.character}(r), \link{nrow}(data))}.
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
#' # subset_to_index(subset, data, parent.frame())
#' # order_to_index(order, data, parent.frame())
#' # append_to_index(append, data, parent.frame())
#' # label_to_character(label, data, parent.frame())
#'
#' @name nse
#' @keywords internal
NULL

#' @rdname nse
subset_to_index <- function(subset, data, enclos = parent.frame(), .subset = NULL) {
  n <- nrow(data)
  if (is.null(.subset)) {
    if (is.null(subset)) {
      return(seq_len(n))
    }
    r <- eval(subset, envir = data, enclos = enclos)
  } else {
    r <- .subset
  }
  stop_if_not(
    is.logical(r),
    length(r) == n,
    m = "`subset` must evaluate to a logical vector of length `nrow(data)`."
  )
  r[is.na(r)] <- FALSE
  r
}

#' @rdname nse
append_to_index <- function(append, data, enclos = baseenv(), .append = NULL) {
  if (is.null(.append)) {
    if (is.null(append)) {
      return(integer(0L))
    }
    l <- as.list(seq_along(data))
    names(l) <- names(data)
    r <- eval(append, envir = l, enclos = enclos)
  } else {
    r <- .append
  }
  m <- match(names(data[r]), names(data), 0L)
  m[m > 0L]
}

#' @rdname nse
order_to_index <- function(order, data, enclos = parent.frame(), .order = NULL) {
  n <- nrow(data)
  if (is.null(.order)) {
    if (is.null(order)) {
      return(seq_len(n))
    }
    r <- eval(order, data, enclos)
  } else {
    r <- .order
  }
  stop_if_not(
    is.numeric(r),
    length(r) == n,
    sort(r) == seq_len(n),
    m = "`order` must be a permutation of `seq_len(nrow(data))`."
  )
  r
}

#' @rdname nse
label_to_character <- function(label, data, enclos = parent.frame(), .label = NULL) {
  n <- nrow(data)
  if (is.null(.label)) {
    if (is.null(label)) {
      return(NULL)
    }
    r <- eval(label, data, enclos)
  } else {
    r <- .label
  }
  stop_if_not(
    is.atomic(r),
    any(length(r) == c(1L, n)),
    m = "`label` must be an atomic vector of length 1 or `nrow(data)`."
  )
  rep_len(as.character(r), n)
}
