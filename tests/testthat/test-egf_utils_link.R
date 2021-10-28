test_that("egf_link_(get|add|remove|extract)", {
  names_top_link0 <- egf_get_names_top(NULL, link = FALSE)
  names_link <- egf_link_get(names_top_link0)

  expect_identical(names_link, ifelse(names_top_link0 == "p", "logit", "log"))
  expect_identical(egf_link_get("invalid name"), NA_character_)
  expect_identical(egf_link_get("log(r)"), NA_character_)

  names_top_link1 <- egf_link_add(names_top_link0)
  expect_identical(egf_link_add(names_top_link0), sprintf("%s(%s)", names_link, names_top_link0))
  expect_identical(egf_link_add("invalid name"), NA_character_)
  expect_identical(egf_link_add("log(r)"), NA_character_)

  expect_identical(egf_link_remove(names_top_link1), names_top_link0)
  expect_identical(egf_link_remove("invalid name"), NA_character_)
  expect_identical(egf_link_remove("r"), NA_character_)

  expect_identical(egf_link_extract(names_top_link1), names_link)
  expect_identical(egf_link_extract("invalid name"), NA_character_)
  expect_identical(egf_link_extract("r"), NA_character_)

  expect_identical(names_top_link1, egf_get_names_top(NULL, link = TRUE))
})

test_that("egf_link_match", {
  expect_identical(egf_link_match("identity"), identity)
  expect_identical(egf_link_match("log"), log)
  expect_identical(egf_link_match("logit"), function(p) qlogis(p), ignore_function_env = TRUE)
  expect_error(egf_link_match("invalid name"))

  expect_identical(egf_link_match("identity", inverse = TRUE), identity)
  expect_identical(egf_link_match("log", inverse = TRUE), exp)
  expect_identical(egf_link_match("logit", inverse = TRUE), function(q) plogis(q), ignore_function_env = TRUE)
  expect_error(egf_link_match("invalid name", inverse = TRUE))
})
