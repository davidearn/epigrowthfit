egf_has_random <-
function(object, check = TRUE) {
	if (check) stopifnot(inherits(object, c("egf", "egf_no_fit")))
	ncol(object[["tmb_out"]][["env"]][["data"]][["Z"]]) > 0L
}

egf_has_converged <-
function(object, check = TRUE, tol = 1) {
	if (check) stopifnot(inherits(object, "egf"))
	object[["optimizer_out"]][["convergence"]] == 0L &&
		is.finite(object[["value"]]) &&
		all(is.finite(range(object[["gradient"]]))) &&
		max(abs(range(object[["gradient"]]))) < tol &&
		object[["hessian"]]
}

egf_par_expand <-
function(obj, par) {
	names(par) <- NULL
	l <- obj[["env"]][["parList"]](par[obj[["env"]][["lfixed"]]()], par)
	if (ncol(obj[["env"]][["data"]][["Z"]]) == 0L)
		l[names(l) != "beta"] <- list(double(0L))
	ans <- unlist1(l)
	attr(ans, "len") <- lengths(l)
	ans
}

egf_par_condense <-
function(obj, par) {
	names(par) <- NULL
	parameters <- obj[["env"]][["parameters"]]
	f <-
	function(x) {
		if (is.null(map <- attr(x, "map"))) {
			ans <- seq_along(x)
			attr(ans, "n") <- length(x)
		}
		else {
			ans <- match(seq_len(attr(x, "nlevels")) - 1L, map)
			attr(ans, "n") <- length(map)
		}
		ans
	}
	index <- lapply(parameters, f)
	len <- vapply(index, attr, 0L, "n")
	l <- split(par, rep.int(gl(length(len), 1L, labels = names(len)), len))
	l <- Map(`[`, l, index)
	ans <- unlist1(l)
	attr(ans, "len") <- lengths(l)
	ans
}

egf_report <-
function(object, check = TRUE) {
	if (check) stopifnot(inherits(object, "egf"))
	e <- object[["tmb_out"]][["env"]]
	if (!exists(".__egf__.", where = e, mode = "environment", inherits = FALSE))
		e[[".__egf__."]] <- new.env(parent = emptyenv())
	ans <- e[[".__egf__."]][["report"]]
	if (is.null(ans)) {
		call <- quote(object[["tmb_out"]][["report"]](object[["best"]]))
		ans <- e[[".__egf__."]][["report"]] <-
			tryCatch(eval(call), error = function(.) `[[<-`(., "call", call))
	}
	if (inherits(ans, "error"))
		stop(ans)
	ans
}

egf_adreport <-
function(object, check = TRUE) {
	if (check) stopifnot(inherits(object, "egf"))
	e <- object[["tmb_out"]][["env"]]
	if (!exists(".__egf__.", where = e, mode = "environment", inherits = FALSE))
		e[[".__egf__."]] <- new.env(parent = emptyenv())
	ans <- e[[".__egf__."]][["adreport"]]
	if (is.null(ans)) {
		if (egf_has_random(object, check = check))
			message("computing a Hessian matrix ...")
		call <- quote(sdreport(object[["tmb_out"]],
		                       par.fixed = object[["best"]][!object[["random"]]],
		                       getReportCovariance = FALSE))
		ans <- e[[".__egf__."]][["adreport"]] <-
			tryCatch(eval(call), error = function(.) `[[<-`(., "call", call))
	}
	if (inherits(ans, "error"))
		stop(ans)
	ans
}

egf_cache <-
function(file, object, topic = NULL, clear = FALSE, clearAll = FALSE, ...) {
	pkgdir <- R_user_dir("epigrowthfit", "cache")
	verdir <- file.path(pkgdir, packageVersion("epigrowthfit"))
	if (clearAll)
		return(unlink(pkgdir, recursive = TRUE))
	stopifnot(identical(file, basename(file)), grepl("^[A-Za-z]", file))
	path <- file.path(verdir, file)
	if (clear)
		return(unlink(path))
	if (file.exists(path))
		return(readRDS(path))
	if (missing(object)) {
		if (is.null(topic)) {
			name <- sub("-\\d+\\.rds$", "", file)
			pattern <- paste0("^", gsub("-", "\\\\W", name), "$")
			hs <- help.search(pattern,
			                  ignore.case = FALSE,
			                  package = "epigrowthfit",
			                  fields = "name",
			                  types = "help",
			                  verbose = FALSE)
			topic <- hs[["matches"]][["Topic"]]
			if (length(topic) != 1L)
				stop(gettextf("resource not created because '%s' matches %d topics",
				              "file", length(topic)),
				     domain = NA)
		}
		example(topic, character.only = TRUE, package = "epigrowthfit",
		        local = TRUE, echo = FALSE)
		if (file.exists(path))
			return(readRDS(path))
		else stop(gettextf("examples for topic \"%s\" were sourced but resource was not created",
		                   topic),
		          domain = NA)
	}
	if (!dir.exists(verdir))
		dir.create(verdir, recursive = TRUE)
	saveRDS(object, file = path, ...)
	object
}
