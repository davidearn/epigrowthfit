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
