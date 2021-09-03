test_that("egf_get_names_top", {
  names_top_all <- egf_get_names_top(NULL, link = FALSE)
  expect_type(names_top_all, "character")
  expect_gt(length(names_top_all), 0L)
  expect_false(anyNA(names_top_all))
  expect_named(names_top_all, NULL)

  model <- egf_model(
    curve = "richards",
    excess = TRUE,
    family = "nbinom",
    day_of_week = TRUE
  )
  names_top <- egf_get_names_top(egf_model(), link = FALSE)
  expect_type(names_top, "character")
  expect_gt(length(names_top), 0L)
  expect_true(all(names_top %in% names_top_all))
})

test_that("egf_has_random", {
  f <- function(Z) {
    e <- new.env()
    e$data <- list(Z = Z)
    object <- list(tmb_out = list(env = e))
    class(object) <- "egf"
    object
  }
  expect_false(egf_has_random(f(matrix(numeric(0L), 6L, 0L))))
  expect_true(egf_has_random(f(matrix(numeric(6L), 6L, 1L))))
})

test_that("egf_combine_frames", {
  object <- list(
    frame_parameters = list(a = data.frame(x = 1:6), b = data.frame(y = letters[1:6])),
    frame_append = data.frame(x = rnorm(1:6), z = TRUE)
  )
  class(object) <- "egf"
  expect_identical(egf_combine_frames(object), data.frame(x = 1:6, y = letters[1:6], z = TRUE))
})
