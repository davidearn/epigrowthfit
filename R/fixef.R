#' Extract fixed effect coefficients
#'
#' Retrieve the coefficients of the fixed effects component
#' of mixed effects models.
#'
#' @param object An `"egf"` object returned by [egf()].
#' @param ... Unused optional arguments.
#'
#' @return
#' A data frame inheriting from class `"egf_fixef"`,
#' with one row per coefficient and variables:
#' \item{`response`}{
#'   A factor naming the relevant response variable,
#'   one of `get_par_names(object, link = TRUE)`.
#' }
#' \item{`term`}{
#'   A factor naming the relevant term in the mixed effects
#'   model formula for `response`.
#' }
#' \item{`contrast`}{
#'   A character vector describing the coefficients.
#'   Means of within-group means are labeled `"mean(term)"`.
#'   Within-group offsets are labeled `"offset(level)`.
#' }
#' \item{`estimate`}{
#'   A numeric vector of coefficients (zero-sum contrasts).
#' }
#' Row names are taken from the `"beta"` component of `object$par`.
#'
#' @aliases fixef
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.egf <- function(object, ...) {
  beta <- object$par[grep("^beta\\[", names(object$par))]
  pn <- get_par_names(object, link = TRUE)

  f <- function(x, s) {
    c(sprintf("mean(%s)", s),
      sprintf("offset(%s)", levels(x)[-nlevels(x)]))
  }

  with(object$tmb_args$data[c("fnl", "ffr")], {
    d <- data.frame(
      response = rep(factor(pn, levels = pn), fnl),
      term = rep(factor(names(fnl), levels = unique(names(fnl))), fnl),
      contrast = unlist(Map(f, ffr[names(fnl)], names(fnl))),
      estimate = beta,
      row.names = names(beta),
      stringsAsFactors = FALSE
    )
    class(d) <- c("egf_fixef", "data.frame")
    d
  })
}
