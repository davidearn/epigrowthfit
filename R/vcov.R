#' @importFrom stats setNames
vcov.egf <- function(object, term_labels = NULL, standardize = TRUE, ...) {
  a <- rownames(object$madf_data$rid)
  if (is.null(a)) {
    stop("Model must include random effects.")
  }
  if (is.null(term_labels)) {
    term_labels <- a
  } else {
    check(term_labels,
      what = "character",
      len = c(1L, NA),
      opt = a,
      "`term_labels` must be a subset of\n",
      "`rownames(object$madf_data$rid)`."
    )
  }
  check(standardize,
    what = "logical",
    len = 1L,
    no = is.na,
    "`standardize` must be TRUE or FALSE."
  )

  wh <- match(term_labels, a)
  dnl <- lapply(wh, function(i) {
    pn <- colnames(object$madf_data$rid)[object$madf_data$rid[i, ] > 0]
    rep(list(sprintf("log_%s", pn)), 2L)
  })

  r <- object$madf_out$report(object$par)
  sdl <- setNames(r$sd_list[wh], a[wh])
  cml <- setNames(r$cor_list[wh], a[wh])
  cml <- mapply("dimnames<-", cml, dnl, SIMPLIFY = FALSE)

  if (standardize) {
    return(cml)
  }
  mapply(function(x, y) outer(x, x) * y, x = sdl, y = cml, SIMPLIFY = FALSE)
}
