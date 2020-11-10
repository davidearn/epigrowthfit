load("../canadacovid_raw.RData")
df <- canadacovid_raw


## Restrict to cumulative incidence by province
df <- data.frame(
  date = as.Date(df$Date),
  time = lubridate::decimal_date(as.Date(df$Date)),
  province = factor(df$Province),
  tot_confirmed = df$confirmed_positive
)

## Use postal abbreviations
levels(df$province) <- sub("PEI", "PE", levels(df$province))

## Split data frame and clean by province
df_split <- split(df, df$province)
df_split <- lapply(df_split, function(x) {

  ## Order by `date` (increasing) and by `tot_confirmed` (decreasing)
  ord <- order(x$date, x$tot_confirmed,
    decreasing = c(FALSE, TRUE),
    na.last = TRUE
  )
  x <- x[ord, ]

  ## Ignore leading missing values
  x <- x[min(which(!is.na(x$tot_confirmed))):nrow(x), ]

  ## Remove duplicated dates:
  ## maximum of `tot_confirmed` is retained due to row order
  x <- x[!duplicated(x$date), ]

  ## Compute interval incidence
  x$new_confirmed <- c(NA, diff(x$tot_confirmed))
  x

})

## Merge and clean
df <- do.call(rbind, df_split)
row.names(df) <- NULL


canadacovid <- df
save(canadacovid, file = "../../data/canadacovid.RData")
