library(epigrowthfit)
options(warn = 2L, error = recover)


## ymd ######
x <- c(0, 0.5, NA, NaN, -Inf, Inf)
Dx <- .Date(x)
names(Dx) <- letters[seq_along(Dx)]

A <- matrix(c(1L, 1L, NA, NA, NA, NA,
              1L, 1L, NA, NA, NA, NA,
              1970L, 1970L, NA, NA, NA, NA),
            nrow = 6L, ncol = 3L,
            dimnames = list(names(Dx), c("d", "m", "y")))

R1 <- ymd(Dx, which = "dmy")
identical(R1, A)
R2 <- ymd(Dx, which = "y", drop = TRUE)
identical(R2, A[, "y", drop = TRUE])


## Dround ######
x <- c(0, 0.5, NA, NaN, -Inf, Inf)
Dx <- .Date(x)
names(Dx) <- letters[seq_along(Dx)]

a1 <- .Date(floor(x))
a2 <- .Date(ceiling(x))
names(a1) <- names(a2) <- names(Dx)
identical(Dfloor(  Dx, "day"), a1)
identical(Dceiling(Dx, "day"), a2)

Dx <- as.Date("1970-07-02")
identical(Dfloor(  Dx, "month"), as.Date("1970-07-01"))
identical(Dceiling(Dx, "month"), as.Date("1970-08-01"))
identical(Dfloor(  Dx, "year"),  as.Date("1970-01-01"))
identical(Dceiling(Dx, "year"),  as.Date("1971-01-01"))

