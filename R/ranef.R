#' Extract random effects conditional modes
#'
#' Retrieve the modes of the coefficients (unit variance scale)
#' of the random effects component of a mixed effects model,
#' conditional on data and parameter estimates.
#'
#' @inheritParams fixef.egf
#'
#' @return
#' A data frame inheriting from class `"egf_ranef"`,
#' with one row per coefficient and variables:
#' \item{`name`}{
#'   Coefficient name in the full parameter vector `object$best`.
#' }
#' \item{`par`}{
#'   Nonlinear model parameter, from
#'   `get_par_names(object, link = TRUE)`.
#' }
#' \item{`term`, `group`}{
#'   Term from random effects component of mixed effects formula
#'   for nonlinear model parameter `par`. `term` and `group`
#'   give the left and right hand sides of the ```|``` operator,
#'   respectively.
#' }
#' \item{`level`}{
#'   Level of the interaction specified by `group`.
#' }
#' \item{`colname`}{
#'   Corresponding column name in the random effects design matrix
#'   `object$tmb_args$data$Z`.
#' }
#' \item{`mode`}{
#'   Conditional mode (unit variance scale).
#' }
#' \item{`sd`}{
#'   Estimated standard deviation.
#' }
#'
#' Attribute `cov` lists estimated random effect covariance matrices.
#' There is one matrix per random effects term and, in a given matrix
#' one row per nonlinear model parameter whose mixed effects formula
#' includes the term.
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
  d <- data.frame(
    name = names(b),
    Z_info[c("par", "term", "group", "level", "colname")],
    mode = b,
    sd = unlist(Map(rep.int, x = object$report$sd_list, times = vapply(object$report$block_list, ncol, 0L))),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(d) <- c("egf_ranef", "data.frame")

  cov <- Map(cor2cov, cor = object$report$cor_list, sd = object$report$sd_list)
  names(cov) <- levels(Z_info$cor)
  cov <- Map(`rownames<-`, cov, tapply(Z_info$par,   Z_info$cor, unique, simplify = FALSE))
  cov <- Map(`colnames<-`, cov, tapply(Z_info$level, Z_info$cor, unique, simplify = FALSE))
  attr(d, "cov") <- cov
  d
}
