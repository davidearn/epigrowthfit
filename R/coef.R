coef.egf <-
function(object, random = TRUE, full = FALSE, ...) {
	stopifnot(is_true_or_false(random), is_true_or_false(full))
	if (full) {
		ans <- egf_expand_par(object[["tmb_out"]], par = object[["best"]])
		len <- attr(ans, "len")
	}
	else {
		ans <- object[["best"]]
		len <- lengths(object[["tmb_out"]][["env"]][["parameters"]])
	}
	map <- lapply(object[["tmb_out"]][["env"]][["parameters"]], attr, "map")
	if (!egf_has_random(object))
		map[names(map) != "beta"] <- list(NULL)
	if (!random) {
		ans <- ans[names(ans) != "b"]
		len <- len[names(len) != "b"]
		map <- map[names(map) != "b"]
	}
	for (i in seq_along(map))
		if (!is.null(m <- map[[i]])) {
			m <- as.integer(m) + 1L
			m[m == 0L] <- NA
			map[[i]] <- m
		}
	attr(ans, "len") <- len
	attr(ans, "map") <- map
	class(ans) <- "coef.egf"
	ans
}

coef.egf_no_fit <-
function(object, ...) {
	object[["best"]] <- object[["init"]]
	coef.egf(object, ...)
}

print.coef.egf <-
function(x, ...) {
	y <- x
	x <- as.double(x)
	names(x) <- names(y)
	NextMethod("print")
	invisible(y)
}

as.list.coef.egf <-
function(x, ...) {
	len <- attr(x, "len")
	map <- attr(x, "map")
	ans <- split(as.double(x),
	             rep.int(gl(length(len), 1L, labels = names(len)), len))
	for (s in names(ans))
		attr(ans[[s]], "map") <- map[[s]]
	ans
}
