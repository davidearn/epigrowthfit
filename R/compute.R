compute_final_size <-
function(R0, S0, I0) {
	suggest("emdbook", "lambertW")
	if (missing(S0)) {
		if (missing(I0)) {
			stop("At least one of 'S0' and 'I0' must be supplied.")
		} else {
			stopifnot(is.numeric(R0), is.numeric(I0))
			S0 <- 1 - I0
		}
	} else {
		if (missing(I0)) {
			stopifnot(is.numeric(R0), is.numeric(S0))
			I0 <- 1 - S0
		} else {
			stopifnot(is.numeric(R0), is.numeric(S0), is.numeric(I0))
		}
	}

	len <- c(length(R0), length(S0), length(I0))
	if (min(len) == 0L) {
		return(double(0L))
	}
	n <- max(len)
	if (length(R0) < n) {
		R0 <- rep_len(R0, n)
	}
	if (length(S0) < n) {
		S0 <- rep_len(S0, n)
	}
	if (length(I0) < n) {
		I0 <- rep_len(I0, n)
	}
	Z <- rep.int(NA_real_, n)

	ok <- !(is.na(R0) | is.na(S0) | is.na(I0) |
	        R0 < 0 | S0 < 0 | I0 < 0 | S0 + I0 > 1)

	if (any(l1 <- ok & R0 == 0)) {
		Z[l1] <- 0
	}
	if (any(l2 <- ok & R0 == Inf)) {
		Z[l2] <- S0[l2]
	}
	if (any(l3 <- ok & !(l1 | l2))) {
		Z[l3] <- S0[l3] + emdbook::lambertW(-R0[l3] * S0[l3] * exp(-R0[l3] * (S0[l3] + I0[l3]))) / R0[l3]
	}
	if (!all(ok)) {
		warning("NA returned for invalid (R0, S0, I0) triples.")
	}
	Z
}

compute_R0 <-
function(r, breaks, probs) {
	stopifnot(is.numeric(r),
	          is.numeric(breaks),
	          length(breaks) >= 2L,
	          is.finite(breaks),
	          diff(breaks) > 0,
	          is.numeric(probs),
	          length(probs) == length(breaks) - 1L,
	          is.finite(probs),
	          probs >= 0,
	          any(probs > 0))
	if (length(r) == 0L) {
		return(double(0L))
	}
	R0 <- r
	R0[] <- NA_real_
	R0[r == 0] <- 1
	R0[r == Inf] <- Inf
	if (any(ok <- is.finite(r) & r > 0)) {
		n <- length(breaks)
		e1 <- exp(tcrossprod(breaks[-n], -r[ok]))
		e2 <- exp(tcrossprod(breaks[-1L], -r[ok]))
		R0[ok] <- r[ok] / colSums((probs / sum(probs)) * (e1 - e2) / (breaks[-1L] - breaks[-n]))
	}
	if (any(r < 0, na.rm = TRUE)) {
		warning("NA returned for negative elements of 'r'.")
	}
	R0
}

compute_tdoubling <-
function(r, per = NULL) {
	stopifnot(is.numeric(r), is.null(per) || is_number(per, "positive"))
	if (any(r < 0, na.rm = TRUE)) {
		r[r < 0] <- NA_real_
		warning("NA returned for negative 'r'.")
	}
	res <- log(2) / r
	class(res) <- c("tdoubling", oldClass(r))
	attr(res, "per") <- per
	res
}

print.tdoubling <-
function(x, ...) {
	y <- x
	if (!is.null(per <- attr(x, "per"))) {
		units <- switch(as.character(per),
		                `1` = "days",
		                `7` = "weeks",
		                `365` = "years (1 year = 365 days)",
		                sprintf("units dt = %g days", per))
		cat("doubling times in ", units, ":\n\n", sep = "")
	}
	class(x) <- NULL
	attr(x, "per") <- NULL
	NextMethod("print")
	invisible(y)
}

as.data.frame.tdoubling <-
function(x, row.names = NULL, optional = FALSE, ...) {
	class(x) <- setdiff(oldClass(x), "tdoubling")
	as.data.frame(x)
}
