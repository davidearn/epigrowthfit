# #' Extract `ADREPORT`ed variables
# #'
# #' Constructs from `"sdreport"` objects lists of variables
# #' passed to TMB macro `ADREPORT` in the C++ template.
# #'
# #' @param sdreport An `"sdreport"` object returned by [TMB::sdreport()].
# #'
# #' @details
# #' When reconstructing matrices from returned vectors,
# #' assume column-major order.
# #'
# #' @return
# #' A named list of data frames with variables `estimate` and `se`
# #' giving estimates and standard errors.
# #'
# #' @keywords internal
# split_sdreport <- function(sdreport) {
#   ssdr <- summary(sdreport, select = "report")
#   colnames(ssdr) <- c("estimate", "se")
#   lapply(split(as.data.frame(ssdr), rownames(ssdr)), `row.names<-`, NULL)
# }
#
# #' Construct combined model frame
# #'
# #' Joins in a single data frame all mixed effects model frames
# #' and further variables specified via [egf()] argument `append`.
# #'
# #' @param object An `"egf"` object returned by [egf()].
# #'
# #' @details
# #' If a variable name occurs in multiple mixed effects model frames,
# #' then only one instance is retained. Except in unusual cases,
# #' all instances of a variable name are identical, and no information
# #' is lost.
# #'
# #' Since the data frames being combined each correspond rowwise to
# #' `object$endpoints`, so does the result.
# #'
# #' @return
# #' A data frame combining (in the sense of [cbind()]) all data frames
# #' listed in `object$frame_par` and the data frame `object$frame_append`.
# #'
# #' @keywords internal
# make_combined <- function(object) {
#   stop_if_not(
#     inherits(object, "egf"),
#     m = "`object` must inherit from class \"egf\"."
#   )
#   frames <- c(unname(object$frame_par), list(object$frame_append))
#   combined <- do.call(cbind, frames)
#   combined[!duplicated(names(combined))]
# }
#
# do_wald <- function(estimate, se, level) {
#   q <- qchisq(level, df = 1)
#   n <- length(estimate)
#   lu <- estimate + rep.int(sqrt(q) * c(-1, 1), c(n, n)) * se
#   dim(lu) <- c(n, 2L)
#   colnames(lu) <- c("lower", "upper")
#   lu
# }
#
# apply_inverse_link <- function(x, g) {
#   if (is.data.frame(x)) {
#     f <- function(x, s) {
#       x[] <- lapply(x, match_link(s, inverse = TRUE))
#       x
#     }
#   } else {
#     f <- function(x, s) {
#       match_link(s, inverse = TRUE)(x)
#     }
#   }
#   x_split <- split(x, g, drop = TRUE)
#   fx_split <- Map(f, x = x_split, s = string_extract_link(names(x_split)))
#   unsplit(fx_split, g, drop = TRUE)
# }


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


