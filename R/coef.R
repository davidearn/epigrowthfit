coef.egf <-
function(object, full = FALSE, ...) {
	stopifnot(is_true_or_false(full))
	if (full)
		res <- egf_expand_par(object$tmb_out, par = object$best)
	else {
		res <- object$best
		attr(res, "lengths") <- lengths(object$tmb_out$env$parameters)
	}
	map <- lapply(object$tmb_out$env$parameters, attr, "map")
	if (!egf_has_random(object))
		map[names(map) != "beta"] <- list(NULL)
	for (i in seq_along(map))
		if (!is.null(map[[i]])) {
			map[[i]] <- as.integer(map[[i]]) + 1L
			map[[i]][map[[i]] == 0L] <- NA
		}
	attr(res, "map") <- map
	attr(res, "full") <- full
	class(res) <- "coef.egf"
	res
}

coef.egf_no_fit <-
function(object, full = FALSE, ...) {
	object$best <- object$init
	coef.egf(object, full = full, ...)
}

print.coef.egf <-
function(x, ...) {
	y <- x
	attributes(x)[c("full", "lengths", "map", "class")] <- NULL
	NextMethod("print")
	invisible(y)
}

as.list.coef.egf <-
function(x, ...) {
	names(x) <- NULL
	len <- attr(x, "lengths")
	f <- rep.int(gl(length(len), 1L, labels = names(len)), len)
	res <- split(x, f)
	map <- attr(x, "map")
	for (s in names(res))
		attr(res[[s]], "map") <- map[[s]]
	attr(res, "full") <- attr(x, "full")
	res
}

as.data.frame.coef.egf <- as.data.frame.vector
