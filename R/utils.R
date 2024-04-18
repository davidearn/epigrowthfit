## A frequent case of 'unlist' usage
unlist1 <-
function(x)
	unlist(x, recursive = FALSE, use.names = FALSE)

## Format probabilities as percentages {a copy of stats:::.format_perc}
formatp <-
function(probs, digits)
	paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
	      "%")

## Creates and prints headings;
## used to section 'print' method output
heading <-
function(text, width = 0.9 * getOption("width"), symbol = ".") {
	line <- strrep(symbol, max(0L, width - 1L - nchar(text)))
	cat(text, " ", line, "\n", sep = "")
}

## Conditionally pluralizes most singular English nouns;
## used to make sure that 'print' and 'plot' method output is grammatical
pluralize <-
function(word, n) {
	plural <- n > 1L
	word[plural] <- paste0(word[plural], "s")
	word
}

## Formats tables of text;
## used by 'print' methods instead of hacking 'print.default'
align <-
function(..., justify = "right", gap = 1L) {
	dots <- list(...) # list of column vectors
	dots <- .mapply(format, list(x = dots, justify = justify), NULL)
	dots[["sep"]] <- strrep(" ", gap)
	do.call(paste, dots)
}

## Disambiguates duplicated strings (or names)
disambiguate <-
function(x, nms = FALSE, fmt = "%s[%d]") {
	if (nms) {
		if (!is.null(nx <- names(x)))
			names(x) <- disambiguate(nx, nms = FALSE, fmt = fmt)
		return(x)
	}
	x <- as.character(x)
	f <- factor(x)
	n <- tabulate(f, nlevels(f))
	i <- unsplit(lapply(n, seq_len), f)
	sprintf(fmt, x, i)
}

if (FALSE) {
## A drop-in replacement for 'rle' that regards NA as equal to previous NA
## and (in double vectors) NaN as equal to previous NaN
rle1 <-
function(x) {
	n <- length(x)
	if (n == 0L)
		return(list(lengths = integer(0L), values = x))
	l <- x[-n] != x[-1L]
	if (any(argna <- is.na(x))) {
		l[is.na(l)] <- FALSE
		if (is.double(x)) {
			argnan <- is.nan(x)
			argna <- argna & !argnan
			l <- l | (argnan[-n] & !argnan[-1L]) | (!argnan[-n] & argnan[-1L])
		}
		l <- l | (argna[-n] & !argna[-1L]) | (!argna[-n] & argna[-1L])
	}
	i <- c(which(l), n)
	list(lengths = diff(c(0L, i)), values = x[i])
}

## Last observation carried forward, with optional replacement of leading NA
locf <-
function(x, x0 = NULL) {
	if (!anyNA(x))
		return(x)
	rle.x <- rle1(x)
	y <- rle.x[["values"]]
	if (is.na(y[1L]) && !is.null(x0))
		y[1L] <- x0
	if (anyNA(y[-1L])) {
		argna.y <- which(c(FALSE, is.na(y[-1L])))
		y[argna.y] <- y[argna.y - 1L]
	}
	rle.x[["values"]] <- y
	inverse.rle(rle.x)
}
}

## Computes confidence intervals from point estimates, standard errors
wald <-
function(value, se, level) {
	h <- 0.5 * (1 - level)
	p <- c(h, 1 - h)
	q <- qnorm(p)
	ans <- value + rep(q, each = length(se)) * se
	dim(ans) <- c(length(se), 2L)
	ans
}

## Real symmetric positive definite matrix   to packed representation
cov2theta <-
function(Sigma) {
	n <- dim(Sigma)[1L]
	R <- chol(Sigma)
	R1 <- R * rep(1 / diag(R, names = FALSE), each = n)
	c(0.5 * log(diag(Sigma, names = FALSE)), R1[upper.tri(R1)])
}

## Real symmetric positive definite matrix from packed representation
theta2cov <-
function(theta) {
	n <- as.integer(round(0.5 * (-1 + sqrt(1 + 8 * length(theta)))))
	h <- seq_len(n)
	R1 <- diag(n)
	R1[upper.tri(R1)] <- theta[-h]
	S <- crossprod(R1)
	scale <- exp(theta[h] - 0.5 * log(diag(S, names = FALSE)))
	scale * S * rep(scale, each = n)
}
