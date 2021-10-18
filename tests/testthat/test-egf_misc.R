test_that("egf_get_names_top", {
  names_top_all <- egf_get_names_top(NULL, link = FALSE)
  expect_type(names_top_all, "character")
  expect_gt(length(names_top_all), 0L)
  expect_false(anyNA(names_top_all))
  expect_named(names_top_all, NULL)

  model <- egf_model()
  names_top <- egf_get_names_top(model, link = FALSE)
  expect_type(names_top, "character")
  expect_gt(length(names_top), 0L)
  expect_true(all(names_top %in% names_top_all))

  object <- list(model = model)
  class(object) <- "egf"
  names_top_again <- egf_get_names_top(object, link = FALSE)
  expect_identical(names_top_again, names_top)
})

test_that("egf_has_random", {
  e <- new.env()
  e$data <- list()
  object <- list(tmb_out = list(env = e))
  class(object) <- "egf"

  e$data$Z <- matrix(numeric(9L), 3L, 3L)
  expect_true(egf_has_random(object))

  e$data$Z <- matrix(numeric(0L), 3L, 0L)
  expect_false(egf_has_random(object))
})
