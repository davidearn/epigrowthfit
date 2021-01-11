#' Extract random effect conditional modes
#'
#' Retrieve random effect modes conditional on the data and
#' parameter estimates.
#'
#' @inheritParams fixef.egf
#'
#' @return
#' A data frame inheriting from class `"egf_ranef"`,
#' with one row per random effect and variables:
#' \item{`response`}{
#'   A factor naming the relevant response variable,
#'   one of `get_par_names(object, link = TRUE)`.
#' }
#' \item{`term`}{
#'   A factor naming the relevant term in the mixed effects
#'   model formula for `response`.
#' }
#' \item{`level`}{
#'   A character vector naming the relevant level of the
#'   interaction specified by `term`.
#' }
#' \item{`mode`}{
#'   A numeric vector of conditional modes.
#' }
#' \item{`sd`}{
#'   A numeric vector of estimated standard deviations.
#' }
#' Row names are taken from the `"b"` component of `object$par`.
#'
#' Attribute `cov` lists estimated random effect covariance matrices.
#' There is one matrix per random effects term and, for a given term,
#' one row per response whose mixed effects model includes the term.
#'
#' @export
ranef <- nlme::ranef

#' @rdname ranef
#' @export
ranef.egf <- function(object, ...) {
  stop_if_not(
    has_random(object),
    m = "`object` must fit a random effects model."
  )

  log_sd_b <- object$par[grep("^log_sd_b\\[", names(object$par))]
  b        <- object$par[grep("^b\\[",        names(object$par))]
  pn <- get_par_names(object, link = TRUE)

  with(object$tmb_args$data[c("rnl", "rfr", "rid")], {
    rnl_times_rid <- rnl * rid
    d <- data.frame(
      response = factor(rep.int(rep.int(pn, nrow(rid)), t(rnl_times_rid)), levels = pn),
      term = factor(rep.int(rownames(rid), rowSums(rnl_times_rid)), levels = rownames(rid)),
      level = unlist(Map(rep.int, lapply(rfr[rownames(rid)], levels), rowSums(rid)), use.names = FALSE),
      mode = b,
      sd = rep.int(exp(log_sd_b), t(rnl_times_rid)[t(rid) > 0L]),
      row.names = names(b),
      stringsAsFactors = FALSE
    )
    ml <- Map(function(x, y) outer(x, x) * y,
      x = object$report$sd_list,
      y = object$report$cor_list
    )
    names(ml) <- rownames(rid)
    dnl <- lapply(seq_along(ml), function(i) {
      pn <- colnames(rid)[rid[i, ] > 0L]
      rep(list(pn), 2L)
    })
    attr(d, "cov") <- Map("dimnames<-", ml, dnl)
    class(d) <- c("egf_ranef", "data.frame")
    d
  })
}
