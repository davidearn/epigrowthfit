#' Extract random effects conditional modes
#'
#' Retrieve the modes of the coefficients (unit variance scale)
#' of the random effects component of a mixed effects model,
#' conditional on data and parameter estimates.
#'
#' @inheritParams fixef.egf
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_ranef"}, with one row per coefficient and variables:
#' \item{name}{
#'   Coefficient name in the full parameter vector \code{object$best}.
#' }
#' \item{par}{
#'   Nonlinear or dispersion model parameter,
#'   from \code{\link{get_par_names}(object, link = TRUE)}.
#' }
#' \item{term, group}{
#'   Random effects term from mixed effects model formula for parameter
#'   \code{par}. \code{term} and \code{group} give the left and right
#'   hand sides of the \code{`|`} operator, respectively.
#' }
#' \item{level}{
#'   Level of interaction \code{group}.
#' }
#' \item{colname}{
#'   Corresponding column name in the random effects design matrix
#'   \code{object$tmb_args$data$Z}.
#' }
#' \item{mode}{
#'   Conditional mode (unit variance scale).
#' }
#' \item{sd}{
#'   Estimated standard deviation.
#' }
#' Attribute \code{cov} lists estimated unstructured random effects
#' covariance matrices. There is one matrix per random effects term.
#' The number of rows of a matrix is equal to the number of nonlinear
#' and dispersion model parameters whose formulae
#' contain the corresponding random effects term.
#'
#' @aliases ranef
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.egf <- function(object, ...) {
  stop_if_not(
    has_random(object),
    m = "`object` must fit a random effects model."
  )
  Z_info <- object$tmb_args$data$Z_info
  b <- object$best[grep("^b\\[", names(object$best))]
  r <- object$tmb_out$report(object$best)
  d <- data.frame(
    name = names(b),
    Z_info[c("par", "term", "group", "level", "colname")],
    mode = b,
    sd = unlist(Map(rep.int, x = r$sd_list, times = vapply(r$block_list, ncol, 0L)), FALSE, FALSE),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(d) <- c("egf_ranef", "data.frame")

  cov <- Map(cor2cov, cor = r$cor_list, sd = r$sd_list)
  names(cov) <- levels(Z_info$cor)
  cov <- Map(`rownames<-`, cov, tapply(Z_info$par,   Z_info$cor, unique, simplify = FALSE))
  cov <- Map(`colnames<-`, cov, tapply(Z_info$level, Z_info$cor, unique, simplify = FALSE))
  attr(d, "cov") <- cov
  d
}
