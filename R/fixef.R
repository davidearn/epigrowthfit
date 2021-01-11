#' Extract fixed effect coefficients
#'
#' Retrieve coefficients (zero-sum contrasts) of the fixed effects
#' component of each mixed effects model.
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
#'   A numeric vector of coefficients.
#' }
#' Row names are taken from the `"beta"` component of `object$par`.
#'
#' @export
fixef <- nlme::fixef

#' @rdname fixef
#' @export
fixef.egf <- function(object, ...) {
  beta <- object$par[grep("^beta\\[", names(object$par))]
  pn <- get_par_names(object, link = TRUE)

  f <- function(x, s) {
    c(sprintf("mean(%s)", s),
      sprintf("offset(%s)", levels(x)[-nlevels(x)]))
  }

  with(object$tmb_args$data[c("fnl", "ffr")], {
    d <- data.frame(
      response = factor(rep.int(pn, fnl), levels = pn),
      term = factor(rep.int(names(fnl), fnl), levels = unique(names(fnl))),
      contrast = unlist(Map(f, ffr[names(fnl)], names(fnl)), use.names = FALSE),
      estimate = beta,
      row.names = names(beta),
      stringsAsFactors = FALSE
    )
    class(d) <- c("egf_fixef", "data.frame")
    d
  })
}
