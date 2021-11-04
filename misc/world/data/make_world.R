(path_covid <- Sys.getenv("PATH_COVID"))
(path_rds   <- Sys.getenv("PATH_RDS"))

covid <- readRDS(path_covid)

## Take first differences
world1 <- covid
world1[["cases_new"]] <- NA_integer_
split(world1[["cases_new"]], world1[["country_iso_alpha3"]]) <-
  c(tapply(world1[["cases_total"]], world1[["country_iso_alpha3"]], function(x) c(NA, diff(x)), simplify = FALSE))

## Treat negative numbers as missing
world1[["cases_new"]][world1[["cases_new"]] < 0L] <- NA

## Aggregate weekly
aggr <- function(data, b) {
  stopifnot(is.integer(b), length(b) == 1L, b > 0L)
  n <- nrow(data)
  if (n < b + 1L) {
    return(data[integer(0L), , drop = FALSE])
  }
  res <- data[seq.int(1L, n, b), , drop = FALSE]
  res[["cases_new"]][1L] <- NA
  r <- nrow(res) - 1L
  i <- seq.int(2L, by = 1L, length.out = r * 7L)
  x <- data[["cases_new"]][i]
  f <- gl(r, 7L)
  res[["cases_new"]][-1L] <- tapply(x, f, sum, na.rm = FALSE)
  res
}
world1_split <- split(world1, world1[["country_iso_alpha3"]])
world7_split <- lapply(world1_split, aggr, b = 7L)

## Delete spurious zeros
delz <- function(data, b, tol) {
  f <- function(x, b, tol) {
    z <- !is.na(x) & x == 0
    if (!any(z)) {
      return(z)
    }
    stopifnot(is.integer(b), length(b) == 1L, b > 0L, b %% 2L == 1L,
              is.numeric(tol), length(tol) == 1L, !is.na(tol))
    p <- rep.int(NA_integer_, 0.5 * (b - 1L))
    X <- embed(c(p, x, p), b)
    z[z] <- apply(X[z, , drop = FALSE] > tol, 1L, any, na.rm = TRUE)
    z
  }
  z <- f(x = data[["cases_new"]], b = b, tol = tol)
  data[!z, , drop = FALSE]
}
world1_split <- lapply(world1_split, delz, b = 15L, tol = 15L)
world7_split <- lapply(world7_split, delz, b = 3L,  tol = 90L)
world1 <- do.call(rbind, world1_split)
world7 <- do.call(rbind, world7_split)
row.names(world1) <- NULL
row.names(world7) <- NULL
saveRDS(list(world1 = world1, world7 = world7), file = path_rds)
str(world1)
str(world7)
