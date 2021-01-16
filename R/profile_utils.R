#' @importFrom stats contr.sum
make_lin_comb <- function(object) {
  p1 <- sum(object$tmb_args$data$fnl) # length(beta)
  p2 <- sum(object$tmb_args$data$rid) # length(log_sd_b)
  p <- p1 + p2
  lc <- matrix(0, nrow = p, ncol = p)
  j <- 1L
  for (nl in object$tmb_args$data$fnl) {
    i <- seq.int(from = j, length.out = nl)
    lc[i, j] <- 1
    if (nl > 1L) {
      lc[i, i[-1L]] <- contr.sum(nl)
    }
    j <- j + nl
  }
  if (p2 > 0L) {
    k <- seq_len(p1)
    lc[-k, -k] <- diag(rep(1, p2))
  }
  lc
}

decontrast_beta_names <- function(x) {
  l <- length(x)
  if (l < 2L) {
    return(x)
  }
  y <- x
  y[-l] <- sprintf("%s+%s", x[1L], x[-1L])
  if (l == 2L) {
    y[l] <- sprintf("%s-%s", x[1L], x[2L])
  } else {
    extract_index <- function(s) {
      gsub("^beta\\[([0-9]+)\\]$", "\\1", s)
    }
    y[l] <- sprintf("%s-sum(beta[%s:%s])",
                    x[1L],
                    extract_index(x[2L]),
                    extract_index(x[l]))
  }
  y
}

make_lin_comb_for_parm <- function(object, parm) {
  fnl <- object$tmb_args$data$fnl
  rid <- object$tmb_args$data$rid

  pn <- get_par_names(object, link = TRUE)
  pn_along_par <- c(rep.int(pn, fnl), rep.int(pn, nrow(rid))[t(rid) > 0L])
  lin_comb <- make_lin_comb(object)

  f <- function(s) {
    i <- which(pn_along_par == s)
    k <- seq_len(fnl[match(s, pn)])
    m <- lin_comb[i, , drop = FALSE]
    rn <- names(object$par)[i]
    rn[k] <- decontrast_beta_names(rn[k])
    rownames(m) <- rn
    m
  }
  do.call(rbind, lapply(unname(parm), f))
}

make_index_for_parm <- function(object, parm) {
  parm <- unname(parm)
  pn <- get_par_names(object, link = TRUE)
  pn_along_par <- with(object$tmb_args$data[c("fnl", "rid")], {
    c(rep.int(pn, fnl), rep.int(pn, nrow(rid))[t(rid) > 0L])
  })

  f <- function(s) {
    i <- which(pn_along_par == s)
    names(i) <- names(object$par)[i]
    i
  }
  unlist(lapply(parm, f))
}

check_parallel <- function(parallel, cores, outfile, cl) {
  if (parallel == "multicore" || (parallel == "snow" && is.null(cl))) {
    stop_if_not_positive_integer(cores, n = 2L)
    if (!is.null(outfile)) {
      stop_if_not_character_string(outfile, n = 2L)
    }
  } else if (parallel == "snow" && !is.null(cl)) {
    stop_if_not(
      inherits(cl, "SOCKcluster"),
      m = "`cl` must be a socket cluster or NULL.",
      n = 2L
    )
  }
  invisible(NULL)
}
