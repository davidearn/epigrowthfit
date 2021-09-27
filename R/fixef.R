#' Extract fixed effect coefficients
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
#' \item{bottom}{
#'   Name of a bottom level mixed effects model parameter
#'   (a fixed effect coefficient).
#' }
#' \item{top}{
#'   Name of the top level nonlinear model parameter whose fitted value
#'   is a function of \code{bottom},
#'   from \code{\link{egf_get_names_top}(object, link = TRUE)}.
#' }
#' \item{term}{
#'   Term from fixed effects component of mixed effects model \link{formula}
#'   for parameter \code{top}.
#' }
#' \item{colname}{
#'   Column name in the fixed effects design matrix
#'   \code{object$tmb_out$env$data$X}.
#' }
#' \item{estimate}{
#'   Coefficient estimate.
#' }
#'
#' @examples
#' example("egf", package = "epigrowthfit", local = TRUE, echo = FALSE)
#' object <- readRDS(system.file("exdata", "egf.rds", package = "epigrowthfit", mustWork = TRUE))
#'
#' zz <- fixef(object)
#' str(zz)
#'
#' @family coefficient extractors
#' @aliases fixef
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.egf <- function(object, ...) {
  res <- data.frame(
    object$info$X[c("bottom", "top", "term", "colname")],
    estimate = object$best[grep("^beta\\[", names(object$best))],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(res) <- c("egf_fixef", "data.frame")
  res
}
