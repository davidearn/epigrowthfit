library(epigrowthfit)
options(warn = 2L, error = recover)


data <- data.frame(x = 1:10,
                   y = letters[11:20],
                   z = .Date(0:9),
                   row.names = letters[1:10])

test_that("egf_eval_subset", {
    expect_identical(egf_eval_subset(quote(x > 5L), data), 6:10)
    expect_identical(egf_eval_subset(NULL, data), seq_len(10L))
    expect_identical(egf_eval_subset(1:3, data), 1:3)
    expect_identical(egf_eval_subset(letters[4:6], data = data), 4:6)
})

test_that("egf_eval_order", {
    expect_identical(egf_eval_order(quote(order(x, decreasing = TRUE)), data),
                     10:1)
    expect_identical(egf_eval_order(NULL, data), seq_len(10L))
    expect_identical(egf_eval_order(c(1:5, 10:6), data), c(1:5, 10:6))
    expect_error(egf_eval_order(1:5, data))
})

test_that("egf_eval_append", {
    expect_identical(egf_eval_append(quote(c(x, y)), data), 1:2)
    expect_identical(egf_eval_append(quote(-x), data), 2:3)
    expect_identical(egf_eval_append(NULL, data), integer(0L))
    expect_identical(egf_eval_append(1:2, data), 1:2)
    expect_identical(egf_eval_append(c("y", "z"), data), 2:3)
})

test_that("egf_eval_label", {
    expect_identical(egf_eval_label(quote(paste("Date:", z)), data),
                     paste("Date:", data$z))
    expect_null(egf_eval_label(NULL, data))
    expect_identical(egf_eval_label("1", data), rep.int("1", 10L))
    expect_identical(egf_eval_label(1, data), rep.int("1", 10L))
})
