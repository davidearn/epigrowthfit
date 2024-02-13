library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

o.1  <- egf_cache(       "egf-1.rds")
o.1f <- egf_cache("fitted-egf-1.rds")


## object ##############################################################

o.1f.e <- data.frame(top = gl(2L, 20L, labels = c("log(r)", "log(c0)")),
                     ts = gl(10L, 2L, 40L, labels = LETTERS[1:10]),
                     window = gl(20L, 1L, 40L, labels = sprintf("window_%02d", 1:20)),
                     estimate = as.double(o.1[["sdreport"]][["value"]]),
                     se = as.double(o.1[["sdreport"]][["sd"]]))
attr(o.1f.e, "se") <- TRUE
class(o.1f.e) <- c("fitted.egf", "data.frame")
stopifnot(identical(o.1f, o.1f.e))


## confint #############################################################

o.1fc <- confint(o.1f, level = 0.95)
o.1fc.e <- o.1f
o.1fc.e[c("lower", "upper")] <-
	epigrowthfit:::wald(o.1f[["estimate"]], o.1f[["se"]], level = 0.95)
o.1fc.e[["se"]] <- NULL
class(o.1fc.e) <- "data.frame"
attr(o.1fc.e, "level") <- 0.95
attr(o.1fc.e, "se") <- NULL
stopifnot(identical(o.1fc, o.1fc.e))
