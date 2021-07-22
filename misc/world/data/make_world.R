covid19 <- readRDS("rds/covid19.rds")
world <- covid19[c("country_iso_alpha3", "Date", "cases_new")]

## Treat negative numbers as missing
world$cases_new[world$cases_new < 0L] <- NA

## Discard time series with fewer than 1000 cases
totals <- c(tapply(covid19$cases_total, covid19$country_iso_alpha3, max, na.rm = TRUE))
world$country_iso_alpha3 <- factor(world$country_iso_alpha3, levels = names(totals)[totals > 1000])
world <- world[!is.na(world$country_iso_alpha3), , drop = FALSE]

## Aggregate weekly
f <- function(d) {
  n <- nrow(d) - 1L
  n <- n - n %% 7L
  d7 <- d[seq.int(1L, 1L + n, by = 7L), , drop = FALSE]
  d7$cases_new[-1L] <- c(tapply(d$cases_new[-1L][seq_len(n)], gl(n / 7L, 7L), sum, na.rm = FALSE))
  d7
}
world_split <- split(world, world$country_iso_alpha3)
world7_split <- lapply(world_split, f)

## Delete spurious zeros
h_ <- function(x, b, tol) {
  zero <- !is.na(x) & x == 0
  if (any(zero)) {
    if (b < 3 || b %% 2 != 1) {
      stop("`b` must be an odd number greater than 1.")
    }
    p <- rep_len(NA_real_, (b - 1) / 2)
    X <- embed(c(p, x, p), b)
    zero[zero] <- apply(X[zero, , drop = FALSE] > tol, 1L, any, na.rm = TRUE)
  }
  !zero
}
h <- function(d, b, tol) {
  ok <- h_(d$cases_new, b, tol)
  d[ok, , drop = FALSE]
}
world_split  <- lapply(world_split,  h, b = 15, tol = 15)
world7_split <- lapply(world7_split, h, b = 3,  tol = 90)
world  <- do.call(rbind, world_split)
world7 <- do.call(rbind, world7_split)
row.names(world)  <- NULL
row.names(world7) <- NULL

saveRDS(list(world = world, world7 = world7), file = "rds/world.rds")
