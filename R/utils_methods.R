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
