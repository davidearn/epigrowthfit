get_par_names <- function(curve = NULL, distr = NULL, excess = NULL) {
  if (inherits(curve, "egf")) {
    pn <- colnames(curve$madf_args$data$fid)
    return(pn)
  }
  if (is.null(curve) && is.null(distr) && is.null(excess)) {
    pn <- c("r", "c0", "thalf", "K", "p", "nbdisp", "b")
  } else {
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
    if (!is.null(excess) && excess) {
      pn <- c(pn, "b")
    }
  }
  sprintf("log_%s", pn)
}

has_random <- function(madf_data) {
  if (inherits(madf_data, "egf")) {
    madf_data <- madf_data$madf_args$data
  }
  nrow(madf_data$rid) > 0L
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
    return(factor(rep("(1)", nrow(frame))))
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
      matrix(1L, nrow = n, ncol = 1L)
    } else {
      unname(model.matrix(f,
        data = list(x = x),
        contrasts.arg = list(x = "contr.sum")
      ))
    }
  }
}

#' @importFrom stats setNames
make_madf_data <- function(frame, sparse_X, sparse_Z) {
  date_to_integer <- function(x) as.integer(x - min(x))

  a <- attributes(frame)
  date <- frame[, a$date_name]
  cases <- frame[, a$cases_name]
  time <- unlist(lapply(split(date, a$index), date_to_integer), use.names = FALSE)
  w <- as.integer(a$index) - 1L

  p <- length(a$fixed_term_labels)
  pn <- names(a$fixed_term_labels)
  pn0 <- get_par_names()

  ftl_incl_dupl <- unlist(a$fixed_term_labels, use.names = FALSE)
  ftl <- unique(ftl_incl_dupl)
  names(ftl) <- ftl
  ffr <- data.frame(lapply(ftl, get_factor, frame = frame), check.names = FALSE)
  fnl <- vapply(ffr, nlevels, integer(1L))[ftl_incl_dupl]
  fmat <- lapply(ffr, factor_to_matrix, sparse = sparse_X, intercept = TRUE)[ftl_incl_dupl]
  X <- do.call(cbind, fmat)

  rtl_incl_dupl <- unlist(a$random_term_labels, use.names = FALSE)
  rtl <- unique(rtl_incl_dupl)
  names(rtl) <- rtl
  rfr <- data.frame(lapply(rtl, get_factor, frame = frame), check.names = FALSE)
  rnl <- vapply(rfr, nlevels, integer(1L))
  rmat <- lapply(rfr, factor_to_matrix, sparse = sparse_Z, intercept = FALSE)
  if (is.null(rtl)) {
    if (sparse_Z) {
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
    spX_flag = 1L * sparse_X,
    spZ_flag = 1L * sparse_Z,
    predict_flag = 0L
  )
  l2 <- setNames(as.list(match(pn0, pn, 0L) - 1L), sprintf("i_%s", pn0))
  c(l1, l2)
}

#' @importFrom stats setNames
make_madf_parameters <- function(madf_data) {
  ts_split <- split(data.frame(madf_data[c("t", "x")]), madf_data$w)

  get_log_c0_log_r <- function(d) {
    d$x[1L] <- 0L
    h <- max(2, nrow(d) / 2)
    x <- d$t
    y <- log1p(cumsum(d$x))
    ab <- tryCatch(
      expr = coef(lm(y ~ x, data = data.frame(x, y), subset = seq_len(h), na.action = na.omit)),
      error = function(e) c(0, 0.1)
    )
    c(ab[[1L]], log(ab[[2L]]))
  }
  log_c0_log_r <- vapply(ts_split, get_log_c0_log_r, numeric(2L))

  pfr <- data.frame(
    log_r      = log_c0_log_r[2L, ],
    log_c0     = log_c0_log_r[1L, ],
    log_thalf  = vapply(ts_split, function(d) log(max(d$t)), numeric(1L)),
    log_K      = vapply(ts_split, function(d) log(2) + log(sum(d$x[-1L], na.rm = TRUE)), numeric(1L)),
    log_p      = 0,
    log_nbdisp = 0,
    log_b      = 0
  )
  ffr <- do.call(rbind, lapply(split(madf_data$ffr, madf_data$w), "[", 1L, , drop = FALSE))
  pn <- colnames(madf_data$fid)

  beta <- unlist(lapply(pn, function(s) {
    tl <- names(madf_data$fnl)[match(s, pn)]
    m <- vapply(split(pfr, ffr[[tl]]), function(d) mean(d[[s]]), numeric(1L))
    mm <- mean(m)
    unname(c(mm, m[-length(m)] - mm))
  }))
  if (has_random(madf_data)) {
    sd_b <- rep(1, sum(madf_data$rid))
    b <- rep(0, sum(madf_data$rnl * rowSums(madf_data$rid)))
  } else {
    sd_b <- NA_real_
    b <- NA_real_
  }
  list(beta = beta, sd_b = sd_b, b = b)
}

make_madf_map <- function(madf_data) {
  if (has_random(madf_data)) {
    return(list())
  }
  list(sd_b = factor(NA), b = factor(NA))
}

make_madf_random <- function(madf_data) {
  if (has_random(madf_data)) {
    return("b")
  }
  NULL
}

get_par_info <- function(par, madf_data, which = c("beta", "sd_b", "b"),
                         decontrast = FALSE) {
  par <- split(par, names(par))
  d <- data.frame(
    name  = character(0),
    par   = character(0),
    term  = character(0),
    level = character(0),
    stringsAsFactors = TRUE
  )

  with(madf_data[c("ffr", "fnl", "rfr", "rnl", "rid")], {
    pn <- colnames(rid)
    if ("beta" %in% which) {
      if (decontrast) {
        f <- function(x, s) levels(x)
      } else {
        f <- function(x, s) {
          c(sprintf(".mean(%s)", s),
            sprintf(".offset(%s)", levels(x)[-nlevels(x)]))
        }
      }
      d <- rbind(d, data.frame(
        name  = rep("beta", length(par$beta)),
        par   = rep(pn, times = fnl),
        term  = rep(names(fnl), times = fnl),
        level = unlist(Map(f, ffr[names(fnl)], names(fnl))),
        stringsAsFactors = TRUE
      ))
    }
    re <- has_random(madf_data)
    if (re && "sd_b" %in% which) {
      d <- rbind(d, data.frame(
        name  = rep("sd_b", length(par$sd_b)),
        par   = rep(pn, nrow(rid))[t(rid) > 0L],
        term  = rep(rownames(rid), rowSums(rid)),
        level = rep(sprintf(".sd(%s)", rownames(rid)), rowSums(rid)),
        stringsAsFactors = TRUE
      ))
    }
    if (re && "b" %in% which) {
      d <- rbind(d, data.frame(
        name  = rep("b", length(par$b)),
        par   = rep(rep(pn, nrow(rid)), t(rnl * rid)),
        term  = rep(rownames(rid), rowSums(rnl * rid)),
        level = unlist(Map(rep, lapply(rfr[rownames(rid)], levels), rowSums(rid))),
        stringsAsFactors = TRUE
      ))
    }
    row.names(d) <- NULL
    d
  })
}

make_lin_comb <- function(object, parm, decontrast) {
  if (!decontrast) {
    index <- which(object$par_info$name != "b" & object$par_info$par %in% parm)
    return(index)
  }

  re <- has_random(object)
  n <- length(object$inr)
  l1 <- (object$par_info$name == "beta")
  l2 <- (object$par_info$name == "sd_b")

  lcl <- lapply(parm, function(s) {
    lp <- (object$par_info$par == s)
    i1 <- which(l1 & lp)
    i2 <- which(l2 & lp)
    m1 <- length(i1)
    m2 <- length(i2)

    lc <- matrix(0, nrow = m1 + m2, ncol = n)
    lc[seq_len(m1), i1[1L]] <- 1
    if (m1 > 1L) {
      lc[seq_len(m1), i1[-1L]] <- contr.sum(m1)
    }
    if (re) {
      lc[m1 + seq_len(m2), i2] <- diag(rep(1, m2))
    }

    list(lc = lc, i = c(i1, i2))
  })

  lin_comb <- do.call(rbind, lapply(lcl, "[[", "lc"))
  index <- do.call(c, lapply(lcl, "[[", "i"))
  par_info <- get_par_info(
    par = object$par,
    madf_data = object$madf_args$data,
    decontrast = TRUE
  )
  structure(lin_comb, index = index, par_info = par_info)
}
