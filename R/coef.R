coef.egf <-
function(object, random = FALSE, full = FALSE, ...) {
	stopifnot(isTrueFalse(random), isTrueFalse(full))
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
		ans <- ans[rep.int(names(len) != "b", len)]
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
	attributes(x) <- NULL
	names(x) <- labels(y)
	NextMethod("print")
	invisible(y)
}

as.list.coef.egf <-
function(x, ...) {
	len <- attr(x, "len")
	map <- attr(x, "map")
	ans <- split(x, rep.int(gl(length(len), 1L, labels = names(len)), len))
	for (s in names(ans))
		attr(ans[[s]], "map") <- map[[s]]
	ans
}

labels.coef.egf <-
function(object, disambiguate = FALSE, ...) {
	len <- attr(object, "len")
	nms <- rep.int(names(len), len)
	if (!disambiguate)
		return(nms)
	map <- attr(object, "map")
	f <-
	function(len, map) {
		if (is.null(map) || length(map) == len)
			seq_len(len)
		else if (length(map) > len)
			## 'map' is a sample with replacement from c(seq_len(len), NA)
			match(unique(if (anyNA(map)) map[!is.na(map)] else map), map)
		else stop("should never happen")
	}
	## FIXME: want fixed width for integer part
	sprintf("%s[%d]", nms, unlist1(.mapply(f, list(len, map), NULL)))
}
