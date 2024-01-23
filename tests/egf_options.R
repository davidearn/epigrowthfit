library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


## egf_model ###########################################################

x <- egf_model()
nms <- c("curve", "excess", "family", "day_of_week")

stopifnot(exprs = {
	is.list(x)
	identical(oldClass(x), "egf_model")
	length(x) == 4L
	identical(names(x), nms)

	is.character(x[["curve"]])
	length(x[["curve"]]) == 1L
	any(x[["curve"]] == eval(formals(egf_model)[["curve"]]))

	is.logical(x[["excess"]])
	length(x[["excess"]]) == 1L
	!is.na(x[["excess"]])

	is.character(x[["family"]])
	length(x[["family"]]) == 1L
	any(x[["family"]] == eval(formals(egf_model)[["family"]]))

	is.integer(x[["day_of_week"]])
	length(x[["day_of_week"]]) == 1L
	any(x[["day_of_week"]] == 0:7)
})


## egf_control #########################################################

x <- egf_control()
nms <- c("optimizer", "inner_optimizer", "trace",
            "profile", "sparse_X", "omp_num_threads")

stopifnot(exprs = {
	is.list(x)
	identical(oldClass(x), "egf_control")
	length(x) == 6L
	identical(names(x), nms)

	is.list(x[["optimizer"]])
	identical(oldClass(x[["optimizer"]]), "egf_optimizer")

	is.list(x[["inner_optimizer"]])
	length(x[["inner_optimizer"]]) > 0L
	vapply(x[["inner_optimizer"]],
	       function(y)
	       	is.list(y) && identical(oldClass(y), "egf_inner_optimizer"),
	       FALSE)

	is.integer(x[["trace"]])
	length(x[["trace"]]) == 1L
	any(x[["trace"]] == 0:2)

	is.logical(x[["profile"]])
	length(x[["profile"]]) == 1L
	!is.na(x[["profile"]])

	is.logical(x[["sparse_X"]])
	length(x[["sparse_X"]]) == 1L
	!is.na(x[["sparse_X"]])

	is.integer(x[["omp_num_threads"]])
	length(x[["omp_num_threads"]]) == 1L
	x[["omp_num_threads"]] > 0L
})


## egf_optimizer #######################################################

x <- egf_optimizer()
nms <- c("f", "args", "control")

stopifnot(exprs = {
	is.list(x)
	identical(oldClass(x), "egf_optimizer")
	length(x) == 3L
	identical(names(x), nms)

	is.function(x[["f"]]) && !is.primitive(x[["f"]])
	identical(names(formals(x[["f"]])), c("par", "fn", "gr", "control", "..."))

	is.list(x[["args"]])
	match(names(x[["args"]]), c("par", "fn", "gr", "control", "..."), 0L) == 0L

	is.list(x[["control"]])
})


## egf_inner_optimizer #################################################

x <- egf_inner_optimizer()
nms <- c("method", "control")

stopifnot(exprs = {
	is.list(x)
	identical(oldClass(x), "egf_inner_optimizer")
	length(x) == 2L
	identical(names(x), nms)

	is.character(x[["method"]])
	length(x[["method"]]) == 1L
	any(x[["method"]] == c("newton", eval(formals(optim)[["method"]])))

	is.list(x[["control"]])
	x[["method"]] != "newton" ||
		match(names(x[["control"]]), c("par", "fn", "gr", "he", "env", "..."), 0L) == 0L
})


## egf_parallel ########################################################

x <- egf_parallel()
nms <- c("method", "outfile", "cores", "args", "cl")

stopifnot(exprs = {
	is.list(x)
	identical(oldClass(x), "egf_parallel")
	length(x) == 5L
	identical(names(x), nms)

	is.character(x[["method"]])
	length(x[["method"]]) == 1L
	any(x[["method"]] == eval(formals(egf_parallel)[["method"]]))

	is.character(x[["outfile"]])
	length(x[["outfile"]]) == 1L
	!is.na(x[["outfile"]])

	is.integer(x[["cores"]])
	length(x[["cores"]]) == 1L
	x[["cores"]] > 0L

	is.list(x[["args"]])

	is.null(x[["cl"]]) ||
		(is.list(x[["cl"]]) && identical(oldClass(x[["cl"]]), c("SOCKcluster", "cluster")))
})


## egf_plot_control ####################################################

reference <- list(window = NULL,
                  data = list(main = NULL, short = NULL, long = NULL),
                  predict = list(estimate = NULL, ci = NULL),
                  asymptote = NULL,
                  box = NULL,
                  axis = list(x = NULL, y = NULL),
                  title = list(main = NULL, sub = NULL,
                               xlab = NULL, ylab = NULL, plab = NULL),
                  tdoubling = c(legend = NULL, estimate = NULL, ci = NULL),
                  heat = c(pal = NULL, bg = NULL, ul = NULL))

recurseOK <-
function(x, reference) {
	if (is.null(reference))
		stopifnot(is.list(x))
	else
		stopifnot(exprs = {
			is.list(x)
			length(x) == length(reference)
			identical(names(x), names(reference))
			mapply(recurseOK, x, reference)
		})
	TRUE
}

x <- egf_plot_control()

stopifnot(exprs = {
	is.list(x)
	identical(oldClass(x), "egf_plot_control")
	recurseOK(x, reference)
})
