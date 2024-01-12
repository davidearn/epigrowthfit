library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


data <- data.frame(x = 1:10,
               y = letters[11:20],
               z = .Date(0:9),
               row.names = letters[1:10])

## egf_eval_subset ######
identical(egf_eval_subset(quote(x > 5L), data), 6:10)
identical(egf_eval_subset(NULL, data), seq_len(10L))
identical(egf_eval_subset(1:3, data), 1:3)
identical(egf_eval_subset(letters[4:6], data = data), 4:6)


## egf_eval_order ######
identical(egf_eval_order(quote(order(x, decreasing = TRUE)), data),
                 10:1)
identical(egf_eval_order(NULL, data), seq_len(10L))
identical(egf_eval_order(c(1:5, 10:6), data), c(1:5, 10:6))
assertError(egf_eval_order(1:5, data))


## egf_eval_append ######
identical(egf_eval_append(quote(c(x, y)), data), 1:2)
identical(egf_eval_append(quote(-x), data), 2:3)
identical(egf_eval_append(NULL, data), integer(0L))
identical(egf_eval_append(1:2, data), 1:2)
identical(egf_eval_append(c("y", "z"), data), 2:3)


## egf_eval_label ######
identical(egf_eval_label(quote(paste("Date:", z)), data),
                 paste("Date:", data$z))
is.null(egf_eval_label(NULL, data))
identical(egf_eval_label("1", data), rep.int("1", 10L))
identical(egf_eval_label(1, data), rep.int("1", 10L))

