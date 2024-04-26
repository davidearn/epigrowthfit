## emdbook::lambertW(z = [-exp(-1), 0], b = 0L), simplified and adapted
Wp <-
function(z, iter.max = 10L, eps = 100 * .Machine[["double.eps"]]) {
	stopifnot(is.double(z), all(is.finite(r <- range(0, z))),
	          r[1L] >= -exp(-1), r[2L] <= 0,
	          is.integer(iter.max), length(iter.max) == 1L, !is.na(iter.max),
	          is.double(eps), length(eps) == 1L, is.finite(eps), eps >= 0)
	w <- sqrt(exp(1) * z + 1) - 1 # a number in [-1, 0]
	iter <- 0L
	done <- FALSE
	while (iter < iter.max) {
		p <- exp(w)
		t <- w * p - z
		f <- w != -1
		w <- w - f * t / (p * (w + f) - 0.5 * (w + 2) * t / (w + f))
		aok <- all(ok <- is.finite(t) & is.finite(w)) # ever not TRUE?
		if (all(abs(if (aok) t else t[ok]) <= eps * (1 + abs(if (aok) w else w[ok])))) {
			done <- TRUE
			break
		}
		iter <- iter - 1L
	}
	if (!done)
		warning(gettextf("iteration limit (%d) reached in Lambert W approximation",
		                 iter.max),
		        domain = NA)
	w
}

finalsize <-
function(R0, S0, I0) {
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
		ans[ok] <- S0 + Wp(-R0 * S0 * exp(-R0 * (S0 + I0))) / R0
	}
	ans
}

R0 <-
function(r, breaks, probs) {
	stopifnot(is.double(r),
	          is.double(breaks),
	          length(breaks) >= 2L,
	          all(is.finite(range(breaks))),
	          min(diff(breaks)) > 0,
	          is.double(probs),
	          length(probs) == length(breaks) - 1L,
	          all(is.finite(range(probs))),
	          min(probs) >= 0,
	          max(probs) >  0)
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
