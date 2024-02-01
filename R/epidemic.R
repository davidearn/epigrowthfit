finalsize <-
function(R0, S0, I0) {
	stopifnot(requireNamespace("emdbook", quietly = TRUE))
	if (missing(S0))
		S0 <- 1 - I0
	if (missing(I0))
		I0 <- 1 - S0
	stopifnot(is.double(R0), is.double(S0), is.double(I0))
	n <- max(length(R0), length(S0), length(I0))
	if (length(R0) < n)
		R0 <- rep_len(R0, n)
	if (length(S0) < n)
		S0 <- rep_len(S0, n)
	if (length(I0) < n)
		I0 <- rep_len(I0, n)
	if (any(nok <- is.na(R0) | is.na(S0) | is.na(I0) | R0 < 0 | S0 < 0 | I0 < 0 | S0 + I0 > 1)) {
		warning(gettextf("NaNs produced: invalid (%s[i], %s[i], %s[i])",
		                 "R0", "S0", "I0"),
		        domain = NA)
		R0[nok] <- NaN
	}
	ans <- R0
	## ans[R0 == 0] <- 0
	ans[R0 == Inf] <- S0[R0 == Inf]
	if (any(ok <- !nok & is.finite(R0) & R0 > 0)) {
		R0 <- R0[ok]
		S0 <- S0[ok]
		I0 <- I0[ok]
		ans[ok] <- S0 + emdbook::lambertW(-R0 * S0 * exp(-R0 * (S0 + I0))) / R0
	}
	ans
}

R0 <-
function(r, breaks, probs) {
	stopifnot(exprs = {
		is.double(r)
		is.double(breaks)
		length(breaks) >= 2L
		all(is.finite(range(breaks)))
		min(diff(breaks)) > 0
		is.double(probs)
		length(probs) == length(breaks) - 1L
		all(is.finite(mm <- range(probs))) && mm[1L] >= 0 && mm[2L] > 0
	})
	if (min(0, r, na.rm = TRUE) < 0) {
		warning(gettextf("NaNs produced: negative %s[i]", "r"),
		        domain = NA)
		r[r < 0] <- NaN
	}
	ans <- r
	ans[r == 0] <- 1
	ans[r == Inf] <- Inf
	if (any(ok <- is.finite(r) & r > 0)) {
		r <- r[ok]
		n <- length(breaks)
		e1 <- exp(tcrossprod(breaks[-n], -r))
		e2 <- exp(tcrossprod(breaks[-1L], -r))
		ans[ok] <- r / colSums((probs / sum(probs)) * (e1 - e2) / (breaks[-1L] - breaks[-n]))
	}
	ans
}

timescale <-
function(r, units) {
	stopifnot(is.double(r))
	if (min(0, r, na.rm = TRUE) < 0) {
		warning(gettextf("NaNs produced: negative %s[i]", "r"),
		        domain = NA)
		r[r < 0] <- NaN
	}
	ans <- 1 / r
	if (missing(units))
		ans
	else .difftime(ans, units)
}
