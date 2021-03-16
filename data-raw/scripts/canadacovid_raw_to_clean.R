load("../canadacovid_raw.RData")
d <- canadacovid_raw

## Restrict to cumulative incidence by province
d <- data.frame(
  date = as.Date(d$Date),
  time = lubridate::decimal_date(as.Date(d$Date)),
  province = factor(d$Province),
  tot_confirmed = d$confirmed_positive
)

## Use postal abbreviations
m <- match("PEI", levels(d$province), 0L)
levels(d$province)[m] <- "PE"

## Split data frame and clean by province
d_split <- split(d, d$province)
d_split <- lapply(d_split, function(x) {

  ## Order by `date` (increasing) and by `tot_confirmed` (decreasing)
  ord <- order(x$date, -x$tot_confirmed)
  x <- x[ord, , drop = FALSE]

  ## Ignore leading missing values
  x <- x[which.max(!is.na(x$tot_confirmed)):nrow(x), , drop = FALSE]

  ## Remove duplicated dates:
  ## maximum of `tot_confirmed` is retained due to row order
  x <- x[!duplicated(x$date), , drop = FALSE]

  ## Compute interval incidence
  x$new_confirmed <- c(NA, diff(x$tot_confirmed))
  x

})

## Merge and clean
d <- do.call(rbind, d_split)
row.names(d) <- NULL


canadacovid <- d
save(canadacovid, file = "../../data/canadacovid.RData")
