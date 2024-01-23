attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## egf_link_(get|add|remove|extract) ###################################

x0 <- egf_get_names_top(NULL, link = FALSE)
x1 <- egf_link_add(x0)
link <- egf_link_get(x0)

identical(link, ifelse(x0 == "p", "logit", "log"))
identical(egf_link_get("invalid name"), NA_character_)
identical(egf_link_get("log(r)"), NA_character_)

identical(x1, sprintf("%s(%s)", link, x0))
identical(egf_link_add("invalid name"), NA_character_)
identical(egf_link_add("log(r)"), NA_character_)

identical(egf_link_remove(x1), x0)
identical(egf_link_remove("invalid name"), NA_character_)
identical(egf_link_remove("r"), NA_character_)

identical(egf_link_extract(x1), link)
identical(egf_link_extract("invalid name"), NA_character_)
identical(egf_link_extract("r"), NA_character_)


## egf_link_match ######################################################

identical(egf_link_match("identity"), identity)
identical(egf_link_match("log"), log)
identical(egf_link_match("logit"),
                 function(p) qlogis(p),
                 ignore_function_env = TRUE)
assertError(egf_link_match("invalid name"))

identical(egf_link_match("identity", inverse = TRUE), identity)
identical(egf_link_match("log", inverse = TRUE), exp)
identical(egf_link_match("logit", inverse = TRUE),
                 function(q) plogis(q),
                 ignore_function_env = TRUE)
assertError(egf_link_match("invalid name", inverse = TRUE))
