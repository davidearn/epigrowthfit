#' Extract fixed effects coefficients
#'
#' Retrieve the coefficients of the fixed effects component
#' of a mixed effects model.
#'
#' @param object An \code{"\link{egf}"} object.
#' @param ... Unused optional arguments.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_fixef"}, with one row per coefficient and variables:
#' \item{name}{
#'   Coefficient name in the full parameter vector \code{object$best}.
#' }
#' \item{par}{
#'   Nonlinear or dispersion model parameter,
#'   from \code{\link{get_par_names}(object, link = TRUE)}.
#' }
#' \item{term}{
#'   Term from fixed effects component of mixed effects \link{formula}
#'   for parameter \code{par}.
#' }
#' \item{colname}{
#'   Corresponding column name in the fixed effects design matrix
#'   \code{object$tmb_args$data$X}.
#' }
#' \item{estimate}{
#'   Coefficient estimate.
#' }
#'
#' @aliases fixef
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.egf <- function(object, ...) {
  beta <- object$best[grep("^beta\\[", names(object$best))]
  d <- data.frame(
    name = names(beta),
    object$tmb_args$data$X_info[c("par", "term", "colname")],
    estimate = beta,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(d) <- c("egf_fixef", "data.frame")
  d
}
