## Load
d_full <- read.csv(file = "2021-03-26/IPHIS_REPORT_2years.CSV")
d <- d_full[c("CASE_REPORTED_DATE", "COUNTY_NAME", "age_grp")]
names(d) <- c("date", "region", "age")

## Clean
d$date <- as.Date(d$date, "%d%b%Y:00:00:00")
d$region <- factor(d$region, exclude = "")
d$age <- factor(d$age, exclude = ".")
levels(d$age) <- sub("^age", "", levels(d$age))
levels(d$age) <- sub("PLUS$", "-Inf", levels(d$age))
d <- d[complete.cases(d), , drop = FALSE]

## Aggregate by region
tt <- table(interaction(d[c("date", "region")], drop = TRUE, sep = "::", lex.order = TRUE))
on_by_region <- data.frame(
  region = factor(substr(names(tt), 13L, nchar(names(tt)))), # strip YYYY-MM-DD::
  date = as.Date(substr(names(tt), 1L, 10L)), # extract YYYY-MM-DD
  cases = as.integer(tt)
)

## Clean and save
o <- do.call(order, on_by_region[c("region", "date")])
on_by_region <- on_by_region[o, , drop = FALSE]
row.names(on_by_region) <- NULL
save(on_by_region, file = "on_by_region.RData")

## Aggregate by age group
tt <- table(interaction(d[c("date", "age")], drop = TRUE, sep = "::", lex.order = TRUE))
on_by_age <- data.frame(
  age = factor(substr(names(tt), 13L, nchar(names(tt)))),
  date = as.Date(substr(names(tt), 1L, 10L)),
  cases = as.vector(tt)
)

## Clean and save
al <- levels(on_by_age$age)
o_al <- order(as.numeric(sub("^([0-9]{1,3})-.*$", "\\1", al)))
on_by_age$age <- ordered(on_by_age$age, levels = al[o_al])
o <- do.call(order, on_by_age[c("age", "date")])
on_by_age <- on_by_age[o, , drop = FALSE]
row.names(on_by_age) <- NULL
save(on_by_age, file = "on_by_age.RData")
