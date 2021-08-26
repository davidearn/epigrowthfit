l <- local({
  latent <- plague_latent_period$relfreq
  m <- length(latent)
  infectious <- plague_infectious_period$relfreq
  n <- length(infectious)
  list(latent = latent, m = m, infectious = infectious, n = n)
})
attach(l, name = "testdata")

test_that("dgi", {
  dgi_ <- function(x) dgi(x = x, latent = latent, infectious = infectious)

  ## Supported on interval [1, m + n)
  expect_equal(dgi_(c(1 - 1e-06, m + n)), c(0, 0))
  ## Constant on intervals [i, i+1)
  expect_equal(dgi_(seq_len(m + n - 1L)), dgi_(seq_len(m + n - 1L) + runif(m + n - 1L, 0, 1)))
  ## Integrates to 1
  expect_equal(sum(dgi_(seq_len(m + n - 1L))), 1)
})

test_that("rgi", {
  set.seed(411422L)
  x <- rgi(1e+8L, latent = latent, infectious = infectious)
  xbin <- .bincode(x, breaks = seq_len(m + n), right = FALSE)
  expect_false(anyNA(xbin))
  freqs <- unname(c(table(xbin))) / length(xbin)
  probs <- dgi(seq_len(m + n - 1L), latent = latent, infectious = infectious)
  expect_equal(freqs, probs, tolerance = 1e-3)
})

detach("testdata")
