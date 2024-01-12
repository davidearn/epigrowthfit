library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


## egf_model ######
x <- egf_model()
is.list(x)
identical(oldClass(x), "egf_model")
length(x) == 4L
names_ <- c("curve", "excess", "family", "day_of_week")
identical(names(x), names_)

is.character(x$curve) && length(x$curve) == 1L
x$curve %in% eval(formals(egf_model)$curve)

is.logical(x$excess) && length(x$excess) == 1L
!is.na(x$excess)

is.character(x$family) && length(x$family) == 1L
x$family %in% eval(formals(egf_model)$family)

is.integer(x$day_of_week) && length(x$day_of_week) == 1L
x$day_of_week %in% 0:7


## egf_control ######
x <- egf_control()
is.list(x)
identical(oldClass(x), "egf_control")
length(x) == 6L
names_ <- c("optimizer", "inner_optimizer", "trace",
            "profile", "sparse_X", "omp_num_threads")
identical(names(x), names_)

is.list(x$optimizer)
identical(oldClass(x$optimizer), "egf_optimizer")

is.list(x$inner_optimizer)
length(x$inner_optimizer) > 0L
for (i in seq_along(x$inner_optimizer)) {
    is.list(x$inner_optimizer[[i]])
    identical(oldClass(x$inner_optimizer[[i]]), "egf_inner_optimizer")
}

is.integer(x$trace) && length(x$trace) == 1L
x$trace %in% 0:2

is.logical(x$profile) && length(x$profile) == 1L
!is.na(x$profile)

is.logical(x$sparse_X) && length(x$sparse_X) == 1L
!is.na(x$sparse_X)

is.integer(x$omp_num_threads) && length(x$omp_num_threads) == 1L
x$omp_num_threads > 0L


## egf_optimizer ######
x <- egf_optimizer()
is.list(x)
identical(oldClass(x), "egf_optimizer")
length(x) == 3L
identical(names(x), c("f", "args", "control"))

is.function(x$f) && !is.primitive(x$f)
names_ <- c("par", "fn", "gr", "control", "...")
identical(names(formals(x$f)), names_)

is.list(x$args)
!any(names(args) %in% names_)

is.list(x$control)


## egf_inner_optimizer ######
x <- egf_inner_optimizer()
identical(oldClass(x), "egf_inner_optimizer")
is.list(x)
length(x) == 2L
identical(names(x), c("method", "control"))

is.character(x$method) && length(x$method) == 1L
x$method %in% c("newton", eval(formals(optim)$method))

is.list(x$control)
if (x$method == "newton") {
    reserved <- c("par", "fn", "gr", "he", "env", "...")
    !any(names(x$control) %in% reserved)
}


## egf_parallel ######
x <- egf_parallel()
is.list(x)
identical(oldClass(x), "egf_parallel")
length(x) == 5L
names_ <- c("method", "outfile", "cores", "args", "cl")
identical(names(x), names_)

is.character(x$method) && length(x$method) == 1L
x$method %in% eval(formals(egf_parallel)$method)

is.character(x$outfile) && length(x$outfile) == 1L
!is.na(x$outfile)

is.integer(x$cores) && length(x$cores) == 1L
x$cores > 0L

is.list(x$args)

if (!is.null(x$cl)) {
    is.list(x$cl)
    identical(oldClass(x$cl), c("SOCKcluster", "cluster"))
}


## egf_plot_control ######
x <- egf_plot_control()
is.list(x)
identical(oldClass(x), "egf_plot_control")
length(x) == 9L
names_ <- c("window", "data", "predict", "asymptote",
            "box", "axis", "title", "tdoubling", "heat")
identical(names(x), names_)

nest <- list(data = c("main", "short", "long"),
             predict = c("estimate", "ci"),
             axis = c("x", "y"),
             title = c("main", "sub", "xlab", "ylab", "plab"),
             tdoubling = c("legend", "estimate", "ci"),
             heat = c("pal", "bg", "ul"))
for (s in names(x)) {
    eval(bquote(is.list(x[[.(s)]])))
}
for (s in names(nest)) {
    eval(bquote(length(x[[.(s)]]) == .(length(nest[[s]]))))
    eval(bquote(identical(names(x[[.(s)]]), .(nest[[s]]))))
}

