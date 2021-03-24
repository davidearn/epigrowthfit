d0 <- read.csv(file = "2021-03-15/IPHIS_REPORT_2years.CSV")

d <- d0[c("CASE_REPORTED_DATE", "COUNTY_NAME", "age_grp")]
names(d) <- c("date", "region", "age")
d$date <- factor(as.Date(d$date, "%d%b%Y:00:00:00"))
d$region <- factor(d$region, exclude = "")
d$age <- factor(d$age, exclude = ".")
levels(d$age) <- sub("^age", "", levels(d$age))
levels(d$age) <- sub("PLUS$", "-Inf", levels(d$age))
d <- d[complete.cases(d), , drop = FALSE]

tt <- table(interaction(d[c("date", "region")], drop = TRUE, sep = "::"))
on_by_region <- data.frame(
  region = factor(substr(names(tt), 13L, nchar(names(tt)))),
  date = as.Date(substr(names(tt), 1L, 10L)),
  cases = as.vector(tt)
)
on_by_region <- on_by_region[order(on_by_region$region, on_by_region$date), ]
row.names(on_by_region) <- NULL
save(on_by_region, file = "on_by_region.RData")

tt <- table(interaction(d[c("date", "age")], drop = TRUE, sep = "::"))
on_by_age <- data.frame(
  age = factor(substr(names(tt), 13L, nchar(names(tt)))),
  date = as.Date(substr(names(tt), 1L, 10L)),
  cases = as.vector(tt)
)
on_by_age$age <- ordered(on_by_age$age, levels = levels(on_by_age$age)[order(as.numeric(sub("^([0-9]{1,3})-.*$", "\\1", levels(on_by_age$age))))])
on_by_age <- on_by_age[order(on_by_age$age, on_by_age$date), ]
row.names(on_by_age) <- NULL
save(on_by_age, file = "on_by_age.RData")
