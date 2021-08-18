test_that("ymd", {
  x <- c(0, 0.5, NA, NaN, -Inf, Inf)
  Dx <- .Date(x)
  names(Dx) <- letters[seq_along(Dx)]

  A <- matrix(c(1L, 1L, NA, NA, NA, NA, 1L, 1L, NA, NA, NA, NA, 1970L, 1970L, NA, NA, NA, NA),
              nrow = 6L, ncol = 3L, dimnames = list(names(Dx), c("d", "m", "y")))

  R1 <- ymd(Dx, which = "dmy")
  expect_equal(R1, A)
  R2 <- ymd(Dx, which = "y", drop = TRUE)
  expect_equal(R2, A[, "y", drop = TRUE])
})

test_that("Dround", {
  x <- c(0, 0.5, NA, NaN, -Inf, Inf)
  Dx <- .Date(x)
  names(Dx) <- letters[seq_along(Dx)]

  a1 <- .Date(floor(x))
  a2 <- .Date(ceiling(x))
  names(a1) <- names(a2) <- names(Dx)
  expect_equal(Dfloor(  Dx, "day"), a1)
  expect_equal(Dceiling(Dx, "day"), a2)

  Dx <- as.Date("1970-07-02")
  expect_equal(Dfloor(  Dx, "month"), as.Date("1970-07-01"))
  expect_equal(Dceiling(Dx, "month"), as.Date("1970-08-01"))
  expect_equal(Dfloor(  Dx, "year"),  as.Date("1970-01-01"))
  expect_equal(Dceiling(Dx, "year"),  as.Date("1971-01-01"))
})
