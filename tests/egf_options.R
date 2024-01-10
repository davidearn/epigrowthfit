test_that("egf_model", {
    x <- egf_model()
    expect_type(x, "list")
    expect_s3_class(x, "egf_model", exact = TRUE)
    expect_length(x, 4L)
    names_ <- c("curve", "excess", "family", "day_of_week")
    expect_named(x, names_, ignore.order = TRUE)

    expect_vector(x$curve, character(0L), 1L)
    expect_true(x$curve %in% eval(formals(egf_model)$curve))

    expect_vector(x$excess, logical(0L), 1L)
    expect_false(is.na(x$excess))

    expect_vector(x$family, character(0L), 1L)
    expect_true(x$family %in% eval(formals(egf_model)$family))

    expect_vector(x$day_of_week, integer(0L), 1L)
    expect_true(x$day_of_week %in% 0:7)
})

test_that("egf_control", {
    x <- egf_control()
    expect_type(x, "list")
    expect_s3_class(x, "egf_control", exact = TRUE)
    expect_length(x, 6L)
    names_ <- c("optimizer", "inner_optimizer", "trace",
                "profile", "sparse_X", "omp_num_threads")
    expect_named(x, names_, ignore.order = TRUE)

    expect_type(x$optimizer, "list")
    expect_s3_class(x$optimizer, "egf_optimizer")

    expect_type(x$inner_optimizer, "list")
    expect_gt(length(x$inner_optimizer), 0L)
    for (i in seq_along(x$inner_optimizer)) {
        expect_type(x$inner_optimizer[[i]], "list")
        expect_s3_class(x$inner_optimizer[[i]], "egf_inner_optimizer")
    }

    expect_vector(x$trace, integer(0L), 1L)
    expect_true(x$trace %in% 0:2)

    expect_vector(x$profile, logical(0L), 1L)
    expect_false(is.na(x$profile))

    expect_vector(x$sparse_X, logical(0L), 1L)
    expect_false(is.na(x$sparse_X))

    expect_vector(x$omp_num_threads, integer(0L), 1L)
    expect_gt(x$omp_num_threads, 0L)
})

test_that("egf_optimizer", {
    x <- egf_optimizer()
    expect_type(x, "list")
    expect_s3_class(x, "egf_optimizer", exact = TRUE)
    expect_length(x, 3L)
    expect_named(x, c("f", "args", "control"), ignore.order = TRUE)

    expect_type(x$f, "closure")
    names_ <- c("par", "fn", "gr", "control", "...")
    expect_named(formals(x$f), names_)

    expect_type(x$args, "list")
    expect_false(any(names(args) %in% names_))

    expect_type(x$control, "list")
})

test_that("egf_inner_optimizer", {
    x <- egf_inner_optimizer()
    expect_s3_class(x, "egf_inner_optimizer", exact = TRUE)
    expect_type(x, "list")
    expect_length(x, 2L)
    expect_named(x, c("method", "control"), ignore.order = TRUE)

    expect_vector(x$method, character(0L), 1L)
    expect_true(x$method %in% c("newton", eval(formals(optim)$method)))

    expect_type(x$control, "list")
    if (x$method == "newton") {
        reserved <- c("par", "fn", "gr", "he", "env", "...")
        expect_false(any(names(x$control) %in% reserved))
    }
})

test_that("egf_parallel", {
    x <- egf_parallel()
    expect_type(x, "list")
    expect_s3_class(x, "egf_parallel", exact = TRUE)
    expect_length(x, 5L)
    names_ <- c("method", "outfile", "cores", "args", "cl")
    expect_named(x, names_, ignore.order = TRUE)

    expect_vector(x$method, character(0L), 1L)
    expect_true(x$method %in% eval(formals(egf_parallel)$method))

    expect_vector(x$outfile, character(0L), 1L)
    expect_false(is.na(x$outfile))

    expect_vector(x$cores, integer(0L), 1L)
    expect_gt(x$cores, 0L)

    expect_type(x$args, "list")

    if (!is.null(x$cl)) {
        expect_type(x$cl, "list")
        expect_s3_class(x$cl, "SOCKcluster")
    }
})

test_that("egf_plot_control", {
    x <- egf_plot_control()
    expect_type(x, "list")
    expect_s3_class(x, "egf_plot_control", exact = TRUE)
    expect_length(x, 9L)
    names_ <- c("window", "data", "predict", "asymptote",
                "box", "axis", "title", "tdoubling", "heat")
    expect_named(x, names_, ignore.order = TRUE)

    nest <- list(data = c("main", "short", "long"),
                 predict = c("estimate", "ci"),
                 axis = c("x", "y"),
                 title = c("main", "sub", "xlab", "ylab", "plab"),
                 tdoubling = c("legend", "estimate", "ci"),
                 heat = c("pal", "bg", "ul"))
    for (s in names(x)) {
        eval(bquote(expect_type(x[[.(s)]], "list")))
    }
    for (s in names(nest)) {
        eval(bquote(expect_length(x[[.(s)]], .(length(nest[[s]])))))
        eval(bquote(expect_named(x[[.(s)]], .(nest[[s]]), ignore.order = TRUE)))
    }
})
