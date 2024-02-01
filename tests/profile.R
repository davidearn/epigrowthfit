library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

o.1  <- egf_cache(        "egf-1.rds")
o.1p <- egf_cache("profile-egf-1.rds")


## object ##############################################################

o.1c <- coef(o.1, full = TRUE)

stopifnot(exprs = {
	is.list(o.1p)
	identical(oldClass(o.1p), c("egf_profile", "data.frame"))
	length(o.1p) == 6L
	identical(names(o.1p), c("top", "ts", "window", "linear_combination", "value", "deviance"))

	typeof(o.1p[["top"]]) == "integer"
	identical(oldClass(o.1p[["top"]]), "factor")
	identical(levels(o.1p[["top"]]), c("log(r)", "log(c0)"))

	identical(o.1p[["ts"]],
	          rep.int(factor("A", levels = LETTERS[1:10]), nrow(o.1p)))
	identical(o.1p[["window"]],
	          rep.int(factor("window_01", levels = sprintf("window_%02d", seq_len(20L))), nrow(o.1p)))
	identical(o.1p[["linear_combination"]],
	          factor(o.1p[["top"]], levels = "log(r)", labels = "1"))

	is.double(o.1p[["value"]])
	is.double(o.1p[["deviance"]])

	identical(unname(c(by(o.1p, o.1p[["linear_combination"]],
	                      function(d) d[["value"]][[which.min(d[["deviance"]])]]))),
	          unname(o.1c[1L]))
})


## confint #############################################################

o.1pc <- confint(o.1p, level = 0.95)
n <- nlevels(o.1p[["linear_combination"]])

o.1pc.e <- structure(o.1p[!duplicated(o.1p[["linear_combination"]]), c("top", "ts", "window"), drop = FALSE],
                     A = attr(o.1p, "A"),
                     x = attr(o.1p, "x"),
                     level = 0.95,
                     row.names = seq_len(n),
                     class = "data.frame")
o.1pc.e[["linear_combination"]] <- seq_len(n)
o.1pc.e[c("estimate", "lower", "upper")] <- list(double(n))

stopifnot(exprs = {
	all.equal(o.1pc, o.1pc.e, tolerance = Inf)
	all(o.1pc[["lower"]] < o.1pc[["estimate"]])
	all(o.1pc[["estimate"]] < o.1pc[["upper"]])
})


## parallel ############################################################

f <-
function(method, cores)
	profile(o.1, top = "log(r)",
	        subset = quote(country == "A" & wave == 1),
	        parallel = egf_parallel(method = method, cores = cores))

windows <- .Platform[["OS.type"]] == "windows"
stopifnot(exprs = {
	all.equal(o.1p, f("multicore", if (windows) 1L else 2L))
	all.equal(o.1p, f("snow", 2L))
})


## plot ################################################################

plot(o.1p, type = "o", bty = "u", las = 1, main = "")
