# #' Recursively merge lists
# #'
# #' Recursively merge elements of a partially specified \link{list}
# #' into a completely specified list of defaults.
# #'
# #' @param x
# #'   A \link{list}.
# #' @param template
# #'   A \link{list} into which \code{x} is merged.
# #'
# #' @details
# #' Recursion proceeds as follows:
# #'
# #' If \code{x} is \code{\link{NULL}}, then \code{template} is replaced
# #' with \code{\link{NULL}}.
# #'
# #' Otherwise, if \code{x} and \code{template} are both \link{list}s
# #' or are both \link{atomic} vectors and \code{template} does not have
# #' \link{names}, then \code{x} is recycled to the length of \code{template}.
# #'
# #' Otherwise, if \code{x} and \code{template} are both \link{list}s
# #' and \code{template} has \link{names}, then elements of \code{x}
# #' are recursively merged into the so-named elements of \code{template}.
# #'
# #' Otherwise, if \code{x} and \code{template} are both \link{atomic} vectors
# #' and \code{template} has \link{names}, then elements of \code{x} replace
# #' the so-named elements of \code{template}.
# #'
# #' Otherwise, \code{x} is discarded and \code{template} is preserved as is.
# #'
# #' @return
# #' A \link{list}.
# #'
# #' @examples
# #' x <- list(
# #'   a = NULL,
# #'   b = 1,
# #'   c = list(
# #'     c1 = list(
# #'       c11 = 1,
# #'       c12 = "a"
# #'     ),
# #'     c2 = c(Jan = 12L, Dec = 1L),
# #'     c3 = FALSE
# #'   )
# #' )
# #' template <- list(
# #'   a = list(
# #'     a1 = 0,
# #'     a2 = "",
# #'     a3 = FALSE
# #'   ),
# #'   b = seq_len(10L),
# #'   c = list(
# #'     c1 = list(
# #'       c11 = 0,
# #'       c12 = "",
# #'       c13 = FALSE
# #'     ),
# #'     c2 = `names<-`(seq_len(12L), month.abb),
# #'     c3 = list(TRUE)
# #'   )
# #' )
# #' # rlmerge(x, template)
# #'
# #' @keywords internal
# rlmerge <- function(x, template) {
#   if (is.null(x)) {
#     return(NULL)
#   }
#   if ((is.list(template) && is.list(x)) ||
#       (is.atomic(template) && is.atomic(x))) {
#     if (is.null(tn <- names(template))) {
#       names(x) <- NULL
#       return(rep_len(x, length(template)))
#     }
#     if (!is.null(xn <- names(x))) {
#       s <- intersect(tn, xn)
#       if (is.list(template)) {
#         template[s] <- Map(rlmerge, x = x[s], template = template[s])
#       } else {
#         template[s] <- x[s]
#       }
#       return(template)
#     }
#   }
#   template
# }
#
# #' Recursively split data frames
# #'
# #' \link[=split]{Split}s a supplied \link[=data.frame]{data frame}
# #' on one or more variables in turn.
# #'
# #' @param x
# #'   A \link[=data.frame]{data frame}.
# #' @param by
# #'   A subset of \code{\link{seq_along}(x)} or \code{\link{names}(x)}
# #'   indicating (in order) variables on which to split \code{x}.
# #' @param drop
# #'   A \link{logical} flag passed to \code{\link{split}}.
# #'
# #' @return
# #' A recursive \link{list} of \link[=data.frame]{data frame}s with
# #' \code{\link{length}(by)} "layers", unless \code{\link{length}(by) = 0},
# #' in which case \code{x} is returned as is.
# #'
# #' @examples
# #' f <- function() sample(letters[1:5], 20L, replace = TRUE)
# #' d <- as.data.frame(replicate(3L, f()))
# #' # rsplit(d, 1:3)
# #'
# #' @keywords internal
# rsplit <- function(x, by = integer(0L), drop = FALSE) {
#   if (length(by) == 0L) {
#     return(x)
#   }
#   l <- split(x, x[[by[1L]]], drop = drop)
#   lapply(l, rsplit, by = by[-1L], drop = drop)
# }


