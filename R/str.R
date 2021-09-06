## Redundant upon CRAN release of TMB version 1.7.21
## https://github.com/kaskr/adcomp/issues/352

#' @export
str.egf <- function(object, ...) {
  class(object$sdreport) <- NULL
  NextMethod("str")
}
