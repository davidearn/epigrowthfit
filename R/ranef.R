#' Extract random effect conditional modes
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
#' \item{cov}{
#'   Name of a covariance matrix.
#'   This is the interaction of \code{term} and \code{group}.
#' }
#' \item{vec}{
#'   Name of a random vector.
#'   This is the interaction of \code{term}, \code{group}, and \code{level}.
#' }
#' \item{bottom}{
#'   Name of a bottom level mixed effects model parameter
#'   (a random effect coefficient).
#' }
#' \item{top}{
#'   Name of the top level nonlinear model parameter whose fitted value
#'   is a function of \code{bottom}.
#' }
#' \item{term, group}{
#'   Random effects term from mixed effects model \link{formula}
#'   for parameter \code{top}. \code{term} and \code{group} give
#'   the left and right hand sides of the \code{`|`} operator.
#' }
#' \item{level}{
#'   \link[=levels]{Level} of \link{factor} or interaction indicated
#'   by \code{group}.
#' }
#' \item{colname}{
#'   Column name in the random effects design matrix
#'   \code{object$tmb_out$env$Z}.
#' }
#' \item{mode}{
#'   Conditional mode (unit variance scale).
#' }
#' \item{sd}{
#'   Estimated standard deviation.
#' }
#' Attribute \code{Sigma} lists estimated unstructured random effect
#' covariance matrices, corresponding to the \link{levels} of \code{cov}.
#'
#' @examples
#' example("egf", "epigrowthfit")
#' zz <- ranef(object)
#' str(zz)
#'
#' @family coefficient extractors
#' @aliases ranef
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.egf <- function(object, ...) {
  stopifnot(egf_has_random(object))

  info <- object$info$Z
  tt <- table(info$top, info$cov)
  d1 <- object$tmb_out$env$data$block_rows
  d2 <- object$tmb_out$env$data$block_cols

  theta <- split(
    unname(object$best)[grep("^theta\\[", names(object$best))],
    rep.int(seq_along(d1), choose(d1 + 1L, 2L))
  )
  theta2cov <- function(theta) {
    n <- as.integer(0.5 * (-1 + sqrt(1 + 8 * length(theta))))
    R <- diag(n)
    R[upper.tri(R)] <- theta[-seq_len(n)]
    RTR <- t(R) %*% R
    diag_D <- exp(theta[seq_len(n)]) / sqrt(diag(RTR))
    RTR[] <- diag_D * RTR * rep(diag_D, each = n)
    RTR
  }
  Sigma <- lapply(theta, theta2cov)
  names(Sigma) <- levels(info$cov)
  for (i in seq_along(Sigma)) {
    dimnames(Sigma[[i]])[1:2] <- list(rownames(tt)[tt[, i] > 0L])
  }

  res <- data.frame(
    object$info$Z[c("cov", "vec", "bottom", "top", "term", "group", "level", "colname")],
    mode = object$best[grep("^b\\[", names(object$best))],
    sd = unlist(Map(rep.int, lapply(Sigma, function(x) sqrt(diag(x))), d2), FALSE, FALSE),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(res) <- c("egf_ranef", "data.frame")
  attr(res, "Sigma") <- Sigma
  res
}
