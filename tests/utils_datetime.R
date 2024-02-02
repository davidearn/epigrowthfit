attach(asNamespace("epigrowthfit"))
options(warn = 2L, error = if (interactive()) recover)


## ymd #################################################################

x <- c(0, 0.5, NaN, -Inf, Inf)
Dx <- .Date(x)
names(Dx) <- letters[seq_along(Dx)]

A <- matrix(c(1L, 1L, NA, NA, NA,
              1L, 1L, NA, NA, NA,
              1970L, 1970L, NA, NA, NA),
            ncol = 3L, dimnames = list(names(Dx), c("d", "m", "y")))
R1 <- ymd(Dx, which = "dmy")
R2 <- ymd(Dx, which = "y", drop = TRUE)

stopifnot(exprs = {
	identical(R1, A)
	identical(R2, A[, "y", drop = TRUE])
})


## ceiling.Date, floor.Date ############################################

x <- c(0, 0.5, NaN, -Inf, Inf)
Dx <- .Date(x)
names(Dx) <- letters[seq_along(Dx)]

a1 <- .Date(floor(x))
a2 <- .Date(ceiling(x))
names(a1) <- names(a2) <- names(Dx)
stopifnot(exprs = {
	identical(  floor.Date(Dx, "day"), a1)
	identical(ceiling.Date(Dx, "day"), a2)
})

Dx <- as.Date("1970-07-02")
stopifnot(exprs = {
	identical(  floor.Date(Dx, "month"), as.Date("1970-07-01"))
	identical(ceiling.Date(Dx, "month"), as.Date("1970-08-01"))
	identical(  floor.Date(Dx,  "year"), as.Date("1970-01-01"))
	identical(ceiling.Date(Dx,  "year"), as.Date("1971-01-01"))
})
