## MJ: priority for refactoring

egf_link_get <-
function(s) {
	ok <- s %in% egf_top(NULL, link = FALSE)
	s[ok] <- replace(rep.int("log", sum(ok)), s[ok] == "p", "logit")
	s[!ok] <- NA
	s
}

egf_link_add <-
function(s) {
	ok <- s %in% egf_top(NULL, link = FALSE)
	s[ok] <- sprintf("%s(%s)", egf_link_get(s[ok]), s[ok])
	s[!ok] <- NA
	s
}

egf_link_remove <-
function(fs) {
	ok <- fs %in% egf_top(NULL)
	fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\2", fs[ok])
	fs[!ok] <- NA
	fs
}

egf_link_extract <-
function(fs) {
	ok <- fs %in% egf_top(NULL)
	fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\1", fs[ok])
	fs[!ok] <- NA
	fs
}

egf_link_match <-
function(f, inverse = FALSE) {
	if (inverse)
		switch(f,
		       identity = identity,
		       log = exp,
		       logit = plogis,
		       stop("link not implemented"))
	else
		switch(f,
		       identity = identity,
		       log = log,
		       logit = qlogis,
		       stop("link not implemented"))
}

egf_top <-
function(object, ...)
	UseMethod("egf_top", object)

egf_top.default <-
function(object, link = TRUE, ...) {
	stopifnot(is.null(object))
	top <- c("r", "alpha", "c0", "tinfl", "K",
	         "p", "a", "b", "disp", paste0("w", 1:6))
	if (link) egf_link_add(top) else top
}

egf_top.egf_model <-
function(object, link = TRUE, ...) {
	top <- switch(object[["curve"]],
	              exponential = c("r", "c0"),
	              subexponential = c("alpha", "c0", "p"),
	              gompertz = c("alpha", "tinfl", "K"),
	              logistic = c("r", "tinfl", "K"),
	              richards = c("r", "tinfl", "K", "a"))
	if (object[["excess"]])
		top <- c(top, "b")
	if (object[["family"]] == "nbinom")
		top <- c(top, "disp")
	if (object[["day_of_week"]] > 0L)
		top <- c(top, paste0("w", 1:6))
	if (link) egf_link_add(top) else top
}

egf_top.egf <-
function(object, link = TRUE, ...)
	egf_top(object[["model"]], link = link, ...)

egf_top.egf_no_fit <- egf_top.egf
