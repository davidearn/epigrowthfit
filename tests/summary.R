library(epigrowthfit)
options(warn = 2L, error = recover)


## basic ######
o <- egf_cache("egf-1.rds")
so <- egf_cache("summary-egf-1.rds")
fo <- egf_cache("fitted-egf-1.rds")

is.list(so)
identical(oldClass(so), "egf_summary")
length(so) == 5L
identical(names(so), c("fitted", "convergence", "value", "gradient", "hessian"))
identical(so$convergence, o$optimizer_out$convergence)
identical(so$value, o$value)
identical(so$gradient, o$gradient)
identical(so$hessian, o$hessian)

is.double(so$fitted)
identical(dim(so$fitted), c(6L, nlevels(fo$top)))
identical(dimnames(so$fitted),
                 list(names(summary(0)), levels(fo$top)))
for (s in colnames(so$fitted)) {
    eval(bquote(all.equal(so$fitted[, .(s)],
                             c(summary(fo$estimate[fo$top == .(s)])))))
}


## print ######
so <- egf_cache("summary-egf-1.rds")
capture.output({
    expect_condition(print(so), regexp = NA)
    identical(print(so), so)
    expect_invisible(print(so))
})

