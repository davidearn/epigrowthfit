library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

latent <- c(0.026, 0.104, 0.182, 0.246, 0.318, 0.104,
            0.013, 0.004, 0.003)
m <- length(latent)

infectious <- c(0.138, 0.462, 0.256, 0.078, 0.041, 0.007,
                0.004, 0.004, 0.006, 0.004)
n <- length(infectious)


## dgi #################################################################

dgi1 <- function(x) dgi(x = x, latent = latent, infectious = infectious)

stopifnot(exprs = {
	## Supported on interval [1, m + n)
	identical(dgi1(c(1 - 1e-06, m + n)), c(0, 0))
	## Constant on intervals [i, i+1)
	identical(dgi1(seq_len(m + n - 1L)),
	          dgi1(seq_len(m + n - 1L) + runif(m + n - 1L, 0, 1)))
	## Integrates to 1
	all.equal(sum(dgi1(seq_len(m + n - 1L))), 1)
})


## rgi #################################################################

set.seed(411422L)
x <- rgi(1e+06L, latent = latent, infectious = infectious)
x.bin <- .bincode(x, breaks = seq_len(m + n), right = FALSE)
stopifnot(!anyNA(x.bin))

freqs <- tabulate(x.bin, m + n - 1L) / length(x.bin)
probs <- dgi(seq_len(m + n - 1L), latent = latent, infectious = infectious)
stopifnot(all.equal(freqs, probs, tolerance = 1e-02))
