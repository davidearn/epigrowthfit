library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


## object ######
o <- egf_cache("egf-1.rds")
po <- egf_cache("profile-egf-1.rds")
co <- coef(o, full = TRUE)

is.list(po)
identical(oldClass(po), c("egf_profile", "data.frame"))
length(po) == 6L
identical(names(po), c("top", "ts", "window", "linear_combination", "value", "deviance"))

is.integer(po$linear_combination)
identical(oldClass(po$linear_combination), "factor")
identical(levels(po$linear_combination), c("1", "2"))

identical(po$top, factor(po$linear_combination, labels = c("log(r)", "log(c0)")))
identical(po$ts, rep.int(factor("A", levels = LETTERS[1:10]), nrow(po)))
identical(po$window, rep.int(factor("window_01", levels = sprintf("window_%02d", seq_len(20L))), nrow(po)))

is.double(po$value)
is.double(po$deviance)

f <- function(d) d$value[[which.min(d$deviance)]]
identical(c(by(po, po$linear_combination, f)), co[1:2],
                 ignore_attr = "names")


## confint ######
po <- egf_cache("profile-egf-1.rds")
cpo <- confint(po, level = 0.95)
n <- nlevels(po$linear_combination)
cpo_expected <- structure(po[!duplicated(po$linear_combination), c("top", "ts", "window"), drop = FALSE],
                          A = attr(po, "A"),
                          x = attr(po, "x"),
                          level = 0.95,
                          row.names = seq_len(n),
                          class = "data.frame")
cpo_expected$linear_combination <- seq_len(n)
cpo_expected[c("estimate", "lower", "upper")] <- list(double(n))
all.equal(cpo, cpo_expected, tolerance = Inf)
all(cpo$lower < cpo$estimate)
all(cpo$estimate < cpo$upper)



## parallel ######
skip_on_cran()
o <- egf_cache("egf-1.rds")
po <- egf_cache("profile-egf-1.rds")

po_multicore <-
    profile(o,
            subset = (country == "A" & wave == 1),
            parallel = egf_parallel(method = "multicore", cores = 2L))
all.equal(po_multicore, po)

## skip_if_not(is.null(pkgload::dev_meta("epigrowthfit")))
po_snow <-
    profile(o,
            subset = (country == "A" & wave == 1),
            parallel = egf_parallel(method = "snow", cores = 2L))
all.equal(po_snow, po)


## plot ######
skip_on_cran()
po <- egf_cache("profile-egf-1.rds")

f <- function() {
    plot(po, type = "o", bty = "u", las = 1, main = "")
}
## vdiffr::expect_doppelganger("plot-egf_profile-1", fig = f)

TRUE # otherwise test is considered skipped

