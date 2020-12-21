get_par_names <- function(curve = NULL, distr = NULL, excess = NULL) {
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

has_random <- function(object) {
  if (inherits(object, "egf")) {
    n <- nrow(object$madf_args$data$rid)
  } else {
    n <- nrow(object$rid)
  }
  n > 0L
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
  s <- unlist(strsplit(term_label, ":"))
  interaction(frame[s], drop = TRUE, sep = ":")
}

#' @importFrom Matrix sparseMatrix sparse.model.matrix
#' @importFrom stats model.matrix
factor_to_matrix <- function(x, sparse, intercept) {
  if (is.null(x)) {
    return(NULL)
  }
  f <- if (intercept) ~x else ~-1 + x
  n <- length(x)
  if (sparse) {
    if (nlevels(x) == 1L) {
      m <- sparseMatrix(seq_len(n), rep(1L, n), x = 1)
    } else {
      m <- unname(sparse.model.matrix(f,
        data = list(x = x),
        contrasts.arg = list(x = "contr.sum")
      ))
    }
  } else {
    if (nlevels(x) == 1L) {
      m <- matrix(1L, nrow = n, ncol = 1L)
    } else {
      m <- unname(model.matrix(f,
        data = list(x = x),
        contrasts.arg = list(x = "contr.sum")
      ))
    }
  }
  m
}

make_madf_data <- function(frame, sparse_X, sparse_Z) {
  date_to_integer <- function(x) as.integer(x - min(x))

  a <- attributes(frame)
  date <- frame[[a$date_name]]
  cases <- frame[[a$cases_name]]
  time <- unlist(lapply(split(date, a$index), date_to_integer), use.names = FALSE)
  w <- as.integer(a$index) - 1L
  wl <- as.vector(table(w))

  p <- length(a$fixed_term_labels)
  pn <- names(a$fixed_term_labels)
  pn0 <- get_par_names()

  ftl_incl_dupl <- unlist(a$fixed_term_labels, use.names = FALSE) # length `p`
  rtl_incl_dupl <- unlist(a$random_term_labels, use.names = FALSE)
  ftl <- unique(ftl_incl_dupl)
  rtl <- unique(rtl_incl_dupl)
  names(ftl) <- ftl
  names(rtl) <- rtl

  ## Term labels can be non-syntactic due to colon delimiters
  ffr <- data.frame(lapply(ftl, get_factor, frame), check.names = FALSE)
  rfr <- data.frame(lapply(rtl, get_factor, frame), check.names = FALSE)

  fnl <- vapply(ffr, nlevels, integer(1L))[ftl_incl_dupl] # length `p`
  rnl <- vapply(rfr, nlevels, integer(1L))

  fmat <- lapply(ffr, factor_to_matrix, sparse = sparse_X, intercept = TRUE)[ftl_incl_dupl] # length `p`
  rmat <- lapply(rfr, factor_to_matrix, sparse = sparse_Z, intercept = FALSE)

  X <- do.call(cbind, fmat)
  if (is.null(rtl_incl_dupl)) {
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
        rid[i, j] <- 1L * (rtl[[i]] %in% a$random_term_labels[[j]])
      }
    }
  }

  l1 <- list(
    t = time,
    x = cases,
    w = w,
    wl = wl,
    X = X,
    Z = Z,
    ffr = ffr,
    rfr = rfr,
    fnl = fnl,
    rnl = rnl,
    fid = fid,
    rid = rid,
    sparseX_flag = 1L * sparse_X,
    sparseZ_flag = 1L * sparse_Z,
    predict_flag = 0L
  )
  l2 <- as.list(match(pn0, pn, 0L) - 1L)
  names(l2) <- sprintf("j_%s", pn0)
  c(l1, l2)
}

make_madf_parameters <- function(madf_data) {
  ts_split  <- split(data.frame(madf_data[c("t", "x")]), madf_data$w)
  ffr_split <- split(madf_data$ffr, madf_data$w)

  get_log_r_log_c0 <- function(d) {
    d$x[1L] <- 0L
    h <- max(2, nrow(d) / 2)
    x <- d$t
    y <- log1p(cumsum(d$x))
    ab <- tryCatch(
      coef(lm(y ~ x, data = data.frame(x, y), subset = seq_len(h), na.action = na.omit)),
      error = function(e) c(0, 0.1)
    )
    c(log(ab[[2L]]), ab[[1L]])
  }
  get_log_thalf <- function(d) log(max(d$t))
  get_log_K     <- function(d) log(2) + log(sum(d$x[-1L], na.rm = TRUE))

  log_r_log_c0 <- vapply(ts_split, get_log_r_log_c0, numeric(2L))
  log_thalf    <- vapply(ts_split, get_log_thalf,    numeric(1L))
  log_K        <- vapply(ts_split, get_log_K,        numeric(1L))

  pn <- colnames(madf_data$fid)
  pl <- list(
    log_r      = log_r_log_c0[1L, ],
    log_c0     = log_r_log_c0[2L, ],
    log_thalf  = log_thalf,
    log_K      = log_K,
    log_p      = 0,
    log_nbdisp = 0,
    log_b      = 0
  )
  pfr <- do.call(data.frame, pl)
  ffr <- do.call(rbind, lapply(ffr_split, "[", 1L, , drop = FALSE))

  beta <- unlist(lapply(pn, function(s) {
    tl <- names(madf_data$fnl)[match(s, pn)]
    m <- vapply(split(pfr[[s]], ffr[[tl]]), mean, numeric(1L))
    mm <- mean(m)
    unname(c(mm, m[-length(m)] - mm))
  }))
  if (has_random(madf_data)) {
    sd <- rep(1, sum(madf_data$rid))
    b <- rep(0, sum(madf_data$rnl * rowSums(madf_data$rid)))
  } else {
    sd <- NA_real_
    b <- NA_real_
  }
  list(beta = beta, sd = sd, b = b)
}

make_madf_map <- function(madf_data) {
  if (has_random(madf_data)) {
    return(list())
  }
  list(sd = factor(NA), b = factor(NA))
}

make_madf_random <- function(madf_data) {
  if (has_random(madf_data)) {
    return("b")
  }
  NULL
}

#' @importFrom stats setNames
rename_par <- function(par) {
  a <- c("beta", "sd", "b")
  which <- a[a %in% names(par)]
  ps <- split(par, names(par))[which]
  f <- function(x, s) setNames(x, sprintf("%s[%d]", s, seq_along(x)))
  do.call(c, Map(f, x = unname(ps), s = names(ps)))
}

get_par_info <- function(object,
                         par = object$par,
                         madf_data = object$madf_args$data,
                         which = c("beta", "sd", "b"),
                         decontrast = FALSE) {
  if (!missing(object) && inherits(object, "egf") && !decontrast) {
    return(object$par_info)
  }

  d <- data.frame(
    name  = character(0),
    par   = character(0),
    term  = character(0),
    level = character(0),
    stringsAsFactors = TRUE
  )

  with(madf_data[c("ffr", "fnl", "rfr", "rnl", "rid")], {
    pn <- colnames(rid)
    ps <- split(par, sub("\\[[0-9]+\\]$", "", names(par)))

    if ("beta" %in% which) {
      if (decontrast) {
        extract_index <- function(s) {
          gsub("^beta\\[([0-9]+)\\]$", "\\1", s)
        }
        f <- function(s) {
          n <- length(s)
          if (n == 1L) {
            return(s)
          }
          ss <- s
          ss[-n] <- sprintf("%s+%s", s[1L], s[-1L])
          ss[n] <- sprintf("%s-sum(beta[%s:%s])", s[1L],  extract_index(s[2L]), extract_index(s[n]))
          ss
        }
        bns <- split(names(ps$beta), rep(seq_along(fnl), fnl))
        dname <- unlist(lapply(unname(bns), f))
        dlevel <- unlist(lapply(ffr[names(fnl)], levels))
      } else {
        f <- function(x, s) {
          c(sprintf("mean(%s)", s),
            sprintf("offset(%s)", levels(x)[-nlevels(x)]))
        }
        dname <- names(ps$beta)
        dlevel <- unlist(Map(f, ffr[names(fnl)], names(fnl)))
      }
      dpar <- rep(pn, times = fnl)
      dterm <- rep(names(fnl), fnl)
      d <- rbind(d, data.frame(
        name  = dname,
        par   = dpar,
        term  = dterm,
        level = dlevel,
        stringsAsFactors = TRUE
      ))
    }
    re <- has_random(madf_data)
    if (re && "sd" %in% which) {
      d <- rbind(d, data.frame(
        name  = names(ps$sd),
        par   = rep(pn, nrow(rid))[t(rid) > 0L],
        term  = rep(rownames(rid), rowSums(rid)),
        level = rep(sprintf("sd(%s)", rownames(rid)), rowSums(rid)),
        stringsAsFactors = TRUE
      ))
    }
    if (re && "b" %in% which) {
      d <- rbind(d, data.frame(
        name  = names(ps$b),
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

split_sdreport <- function(object) {
  ssdr <- summary(object, select = "report")
  f <- factor(rownames(ssdr))
  rownames(ssdr) <- NULL
  colnames(ssdr) <- c("estimate", "se")
  split(as.data.frame(ssdr), f)
}

make_lin_comb <- function(object, parm, decontrast) {
  l1 <- grepl("^beta\\[", names(object$par))
  l2 <- grepl("^sd\\[", names(object$par))
  lp <- (object$par_info$par %in% parm)

  if (!decontrast) {
    index <- which((l1 | l2) & lp)
    return(index)
  }

  re <- has_random(object)
  n <- length(object$nonrandom)

  f <- function(s) {
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
  }

  lcl <- lapply(parm, f)
  lin_comb <- do.call(rbind, lapply(lcl, "[[", "lc"))
  index <- do.call(c, lapply(lcl, "[[", "i"))
  par_info <- get_par_info(object, decontrast = TRUE)
  structure(lin_comb, index = index, par_info = par_info)
}
