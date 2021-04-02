l <- list(
  TORONTO = list(
    c("2020-03-05", "2020-04-23"),
    c("2020-08-07", "2020-09-05"),
    c("2020-09-06", "2020-10-05"),
    c("2020-10-11", "2021-01-13")
  ),
  PEEL = list(
    c("2020-03-14", "2020-04-23"),
    c("2020-08-02", "2020-10-10"),
    c("2020-10-16", "2020-11-19"),
    c("2020-12-23", "2021-01-08")
  ),
  YORK = list(
    c("2020-03-09", "2020-04-18"),
    c("2020-08-12", "2020-10-26"),
    c("2020-12-12", "2021-01-02")
  ),
  OTTAWA = list(
    c("2020-03-09", "2020-04-29"),
    c("2020-07-01", "2020-07-31"),
    c("2020-08-31", "2020-10-09"),
    c("2020-12-19", "2021-01-17")
  ),
  DURHAM = list(
    c("2020-03-17", "2020-04-25"),
    c("2020-08-16", "2020-09-29"),
    c("2020-10-02", "2021-01-13")
  ),
  HALTON = list(
    c("2020-03-19", "2020-04-10"),
    c("2020-08-30", "2020-10-19"),
    c("2020-10-25", "2020-11-13"),
    c("2020-11-24", "2021-01-07")
  ),
  HAMILTON = list(
    c("2020-03-10", "2020-04-15"),
    c("2020-09-09", "2020-10-13"),
    c("2020-10-16", "2021-01-06")
  ),
  WATERLOO = list(
    c("2020-03-03", "2020-04-25"),
    c("2020-09-04", "2020-09-23"),
    c("2020-10-24", "2020-11-27"),
    c("2020-12-20", "2021-01-13")
  ),
  SIMCOE = list(
    c("2020-03-11", "2020-04-23"),
    c("2020-08-20", "2020-10-12"),
    c("2020-10-16", "2021-01-19")
  ),
  MIDDLESEX = list(
    c("2020-03-16", "2020-04-17"),
    c("2020-11-15", "2021-01-10")
  ),
  NIAGARA = list(
    c("2020-03-12", "2020-04-16"),
    c("2020-09-01", "2020-10-03"),
    c("2020-10-15", "2020-11-10"),
    c("2020-11-26", "2021-01-14")
  ),
  ESSEX = list(
    c("2020-03-15", "2020-04-08"),
    c("2020-10-20", "2021-01-02")
  )
)
l <- l[sort(names(l))]
se <- matrix(unlist(l), ncol = 2L, byrow = TRUE, dimnames = list(NULL, c("start", "end")))
se <- as.data.frame(se, stringsAsFactors = FALSE)
se[] <- lapply(se, as.Date)
region <- rep.int(gl(length(l), 1L, labels = names(l)), lengths(l))
f <- function(x) x[do.call(order, x), , drop = FALSE]
se <- unsplit(by(se, region, f, simplify = FALSE), region)
wave <- factor(unsplit(lapply(l, seq_along), region))
load("populations.RData")
population <- populations[match(region, names(populations), 0L)]
endpoints <- data.frame(se, region, wave, population, row.names = NULL)
save(endpoints, file = "endpoints.RData")
