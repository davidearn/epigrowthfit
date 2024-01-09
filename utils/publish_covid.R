canadacovid_raw <- read.csv("canadacovid_raw.csv")

## Extract cumulative incidence by province
canadacovid <- data.frame(
  province = factor(canadacovid_raw$Province),
  Date = as.Date(canadacovid_raw$Date),
  cases_tot = canadacovid_raw$confirmed_positive
)

## Replace non-postal abbreviation
m <- match("PEI", levels(canadacovid$province), 0L)
levels(canadacovid$province)[m] <- "PE"

## Order
o <- order(canadacovid$province, canadacovid$Date, -canadacovid$cases_tot)
canadacovid <- canadacovid[o, , drop = FALSE]

## Process by province
f <- function(d) {
  ## Ignore leading missing values and duplicated dates
  k <- seq.int(which.max(!is.na(d$cases_tot)), nrow(d))
  k <- k[!duplicated(d$Date[k])]
  d <- d[k, , drop = FALSE]

  ## Compute interval incidence
  d$cases_new <- c(NA, diff(d$cases_tot))
  d
}
canadacovid <- do.call(rbind, by(canadacovid, canadacovid$province, f, simplify = FALSE))
row.names(canadacovid) <- NULL

## Save
save(canadacovid, file = "../data/canadacovid.RData", compress = "xz")
