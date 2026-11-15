library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


## time origin #########################################################

## The model measures time from the start of the fitting window, so the
## origin of the time variable should not matter.  In particular, Date
## times (days since 1970-01-01) should give the same fit as times
## measured from the start of the fitting window.

data(covid19.ontario, package = "epigrowthfit")
x <- c(NA, diff(covid19.ontario[["cases"]]))
D <- covid19.ontario[["date"]]
D0 <- as.Date("2021-12-14")
D1 <- as.Date("2022-01-11")

f <-
function(time, start, end, ...)
	egf(egf_model(...),
	    formula_ts = cbind(time, x) ~ 1,
	    formula_windows = cbind(start, end) ~ 1,
	    data_ts = data.frame(time = time, x = x),
	    data_windows = data.frame(start = start, end = end),
	    na_action_ts = "pass", se = TRUE)

for (args in list(list(curve = "logistic"),
                  list(curve = "logistic", excess = TRUE),
                  list(curve = "exponential"))) {
	m0 <- do.call(f, c(list(time = as.double(D - D0), start = 0,
	                        end = as.double(D1 - D0)), args))
	m1 <- do.call(f, c(list(time = D, start = D0, end = D1), args))
	stopifnot(exprs = {
		egf_has_converged(m0)
		egf_has_converged(m1)
		all.equal(coef(m0), coef(m1), tolerance = 1e-06)
	})
}
