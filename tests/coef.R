library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

o.2  <- egf_cache("egf-2.rds")
nms <- c("beta", "theta", "b")


## object ##############################################################

o.2c <- coef(o.2)
o.2c.e <- structure(as.double(o.2c),
                    names = rep.int(nms, c(2L, 2L, 40L)),
                    lengths = c(beta = 2L, theta = 2L, b = 40L),
                    map = list(beta = NULL, theta = c(1L, 2L, NA), b = NULL),
                    full = FALSE,
                    class = "coef.egf")
stopifnot(identical(o.2c, o.2c.e))


## print ###############################################################

vv <- withVisible(print(o.2c))
stopifnot(exprs = {
	identical(vv[["value"]], o.2c)
	identical(vv[["visible"]], FALSE)
})


## as.list #############################################################

o.2cl <- as.list(o.2c)
o.2cl.e <- split(as.double(o.2c), factor(names(o.2c), levels = nms))
for (s in nms)
	attr(o.2cl.e[[s]], "map") <- attr(o.2c, "map")[[s]]
attr(o.2cl.e, "full") <- attr(o.2cl, "full")
stopifnot(identical(o.2cl, o.2cl.e))
