get_par_names <- function(curve = NULL,
                          distr = NULL,
                          include_baseline = NULL) {
  if (is.null(curve) && is.null(distr) && is.null(include_baseline)) {
    return(c("r", "c0", "thalf", "K", "p", "nbdisp", "b"))
  }
  pn <- character(0)
  if (!is.null(curve)) {
    a <- switch(curve,
      exponential = c("r", "c0"),
      logistic    = c("r", "thalf", "K"),
      richards    = c("r", "thalf", "K", "p")
    )
    pn <- c(pn, a)
  }
  if (!is.null(distr) && distr == "nbinom") {
    pn <- c(pn, "nbdisp")
  }
  if (!is.null(include_baseline) && include_baseline) {
    pn <- c(pn, "b")
  }
  pn
}

#' @importFrom stats terms reformulate
get_term_labels <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  tl <- labels(terms(formula))
  if (length(tl) == 0L) { # ~1
    return("(1)")
  }
  has_bar <- grepl("|", tl, fixed = TRUE)
  if (any(has_bar)) {
    tl[has_bar] <- gsub("1 \\| ", "", tl[has_bar])
    tl <- labels(terms(reformulate(tl)))
  }
  tl
}

get_factor <- function(term_label, frame) {
  if (is.null(term_label)) {
    return(NULL)
  }
  if (term_label == "(1)") {
    return(factor(integer(nrow(frame))))
  }
  interaction(frame[unlist(strsplit(term_label, ":"))], drop = TRUE, sep = ":")
}

#' @importFrom Matrix sparseMatrix sparse.model.matrix
factor_to_matrix <- function(x, sparse, intercept) {
  if (is.null(x)) {
    return(NULL)
  }
  f <- if (intercept) ~x else ~-1 + x
  n <- length(x)
  l <- nlevels(x)
  if (sparse) {
    if (l == 1L) {
      sparseMatrix(seq_len(n), rep(1L, n), x = 1)
    } else {
      unname(sparse.model.matrix(f,
        data = list(x = x),
        contrasts.arg = list(x = "contr.sum")
      ))
    }
  } else {
    if (l == 1L) {
      matrix(1L, nrow = n, ncol = 1)
    } else {
      unname(model.matrix(f,
        data = list(x = x),
        contrasts.arg = list(x = "contr.sum")
      ))
    }
  }
}

#' @importFrom stats setNames
make_madf_data <- function(frame, spX, spZ) {
  date_to_integer <- function(x) as.integer(x - min(x))

  a <- attributes(frame)
  date <- frame[, a$date_name]
  cases <- frame[, a$cases_name]
  time <- unlist(lapply(split(date, a$index), date_to_integer), use.names = FALSE)
  w <- structure(a$index, class = NULL, levels = NULL) - 1L

  p <- length(a$fixed_term_labels)
  pn <- names(a$fixed_term_labels)
  pn0 <- get_par_names()

  ftl_incl_dupl <- unlist(a$fixed_term_labels, use.names = FALSE)
  ftl <- unique(ftl_incl_dupl)
  names(ftl) <- ftl
  ffr <- data.frame(lapply(ftl, get_factor, frame = frame), check.names = FALSE)
  fnl <- sapply(ffr, nlevels)[ftl_incl_dupl]
  fmat <- lapply(ffr, factor_to_matrix, sparse = spX, intercept = TRUE)[ftl_incl_dupl]
  X <- do.call(cbind, fmat)

  rtl_incl_dupl <- unlist(a$random_term_labels, use.names = FALSE)
  rtl <- unique(rtl_incl_dupl)
  names(rtl) <- rtl
  rfr <- data.frame(lapply(rtl, get_factor, frame = frame), check.names = FALSE)
  rnl <- vapply(rfr, nlevels, integer(1L))
  rmat <- lapply(rfr, factor_to_matrix, sparse = spZ, intercept = FALSE)
  if (is.null(rtl)) {
    if (spZ) {
      Z <- sparseMatrix(integer(0L), integer(0L), x = numeric(0L), dims = c(nrow(frame), 0L))
    } else {
      Z <- matrix(integer(0L), nrow = nrow(frame))
    }
  } else {
    Z <- do.call(cbind, rmat)
  }

  fid <- matrix(0L, nrow = length(ftl), ncol = p, dimnames = list(ftl, pn))
  for (i in seq_len(nrow(fid))) { # FIXME: outer(FUN = "==")?
    for (j in seq_len(ncol(fid))) {
      fid[i, j] <- 1L * (ftl[[i]] == a$fixed_term_labels[[j]])
    }
  }
  if (is.null(rtl)) {
    rid <- matrix(integer(0L), ncol = p, dimnames = list(NULL, pn))
  } else {
    rid <- matrix(0L, nrow = length(rtl), ncol = p, dimnames = list(rtl, pn))
    for (i in seq_len(nrow(rid))) { # FIXME: outer(FUN = "%in%")?
      for (j in seq_len(ncol(rid))) {
        rid[i, j] <- rtl[[i]] %in% a$random_term_labels[[j]]
      }
    }
  }

  l1 <- list(
    t = time,
    x = cases,
    w = w,
    X = X,
    Z = Z,
    ffr = ffr,
    rfr = rfr,
    fnl = fnl,
    rnl = rnl,
    fid = fid,
    rid = rid,
    spX_flag = 1L * spX,
    spZ_flag = 1L * spZ,
    predict_flag = 0L
  )
  l2 <- setNames(as.list(match(pn0, pn) - 1L), sprintf("i_log_%s", pn0))
  c(l1, l2)
}

#' @importFrom stats setNames
make_madf_parameters(madf_data) {
  pn <- colnames(madf_data$fid)
  pn0 <- get_par_names()
  pid <- setNames(match(pn0, pn), pn0)

  ffr_aug <- data.frame(madf_data[c("t", "x", "ffr")])
  names(ffr_aug)[-(1:2)] <- names(madf_data$ffr) # sigh
  ffr_aug_split <- split(ffr_aug, madf_data$w)

  get_naive_log_r_log_c0 <- function(d) {
    h <- max(2, nrow(d) / 2)
    x <- d$t[seq_len(h)]
    y <- log1p(cumsum(c(0, d$x[-1])[seq_len(h)]))
    ab <- try(silent = TRUE, expr = {
      coef(lm(y ~ x, data = data.frame(x, y), na.action = na.omit))
    })
    if (inherits(ab, "try-error")) {
      return(log(c(1, 0.1)))
    }
    c(ab[[1L]], log(ab[[2L]]))
  }

  log_r_log_c0 <- sapply(ffr_aug_split, get_naive_log_r_log_c0)
  log_r <- log_r_log_c0[2L, ]
  log_c0 <- log_r_log_c0[1L, ]
  log_thalf <- sapply(ffr_aug_split, function(d) log(max(d$t)))
  log_K <- sapply(ffr_aug_split, function(d) log(2) + log(sum(d$x[-1L], na.rm = TRUE)))
  log_p <- log_nbdisp <- log_b <- 0 * log_r

  pfr <- do.call(cbind, as.list(environment())[sprintf("log_%s", pn)])
  pfr <- cbind(pfr,
    do.call(rbind, lapply(ffr_aug_split, "[", 1L, -(1:2), drop = FALSE))
  )

  beta <- unlist(lapply(pn, function(s) {
    tl <- names(madf_data$fnl)[pid[s]]
    ls <- sprintf("log_%s", s)
    unname(sapply(split(pfr, pfr[[tl]]), function(d) mean(d[[ls]])))
  }))
  if (nrow(madf_data$rid) > 0L) {
    b <- rep(0, sum(madf_data$rnl * rowSums(madf_data$rid)))
    sd_b <- rep(1, sum(madf_data$rid))
  } else {
    b <- NA_real_
    sd_b <- NA_real_
  }
  list(beta = beta, b = b, sd_b = sd_b)
}

make_madf_map <- function(madf_data) {
  if (nrow(madf_data$rid) > 0L) {
    return(list())
  }
  list(b = factor(NA), sd_b = factor(NA))
}

make_madf_random <- function(madf_data) {
  if (nrow(madf_data$rid) > 0L) {
    return(NULL)
  }
  "b"
}
