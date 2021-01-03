#' @keywords internal
make_lin_comb <- function(fnl, rid) {
  p1 <- sum(fnl) # length(beta)
  p2 <- sum(rid) # length(log_sd_b)
  p <- p1 + p2
  lc <- matrix(0, nrow = p, ncol = p)
  j <- 1L
  for (nl in fnl) {
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

#' @keywords internal
make_index_for_parm <- function(object, parm) {
  i_parm <- which(object$par_info$response %in% parm)
  i_parm_split <- split(i_parm, object$par_info$vector[i_parm])
  index <- unlist(i_parm_split[c("beta", "log_sd_b")], use.names = FALSE)
  names(index) <- object$par_info$name[index]
  index
}

#' @keywords internal
make_lin_comb_for_parm <- function(object, parm, decontrast) {
  i_parm <- which(object$par_info$response %in% parm)
  i_parm_split <- split(i_parm, object$par_info$vector[i_parm])
  index <- unlist(i_parm_split[c("beta", "log_sd_b")], use.names = FALSE)
  index_prime <- unlist(i_parm_split[c(".decontrast(beta)", "log_sd_b")], use.names = FALSE)

  lin_comb <- make_lin_comb(fnl = object$tmb_args$data$fnl,
                            rid = object$tmb_args$data$rid)
  lin_comb <- lin_comb[index, , drop = FALSE]
  rownames(lin_comb) <- object$par_info$name[index_prime]
  lin_comb
}
