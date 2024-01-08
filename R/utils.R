### A frequent case of 'unlist' usage
unlist1 <- function(x) {
    unlist(x, recursive = FALSE, use.names = FALSE)
}

### For backwards compatibility
if (getRversion() < "4.0.0") {
    deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
}

### Coerces a sequence of R objects to character, concatenates them with
### no separator, and formats the resulting string so that, when printed,
### the user sees nicely wrapped paragraphs. Used inside of 'stop', etc.
### so that condition messages print nicely.
wrap <- function(..., width = 0.9 * getOption("width")) {
    dots <- list(...)
    x <- vapply(dots, paste0, "", collapse = "")
    x1 <- paste0(x, collapse = "")
    y <- unlist(strsplit(x1, "\n[ \t\n]*\n"), FALSE, FALSE)
    y1 <- paste0(y, collapse = "\n\n")
    z <- strwrap(y1, width = width)
    paste0(z, collapse = "\n")
}

### Tests that a package
### (usually one listed in the 'Suggests' field of DESCRIPTION)
### is installed and throws error if not.
suggest <- function(package, name = NULL) {
    if (!requireNamespace(package, quietly = TRUE)) {
        name <- if (is.null(name)) package else paste0(package, "::", name)
        stop1(sQuote(name), " is needed, but ", sQuote(package),
              " is not installed. Run ",
              sQuote(paste0("install.packages(", dQuote(package), ")")),
              " to install ", sQuote(package), ", then try again.")
    }
    invisible(NULL)
}

### Creates and prints headings. Used to section 'print' method output.
heading <- function(text, width = 0.9 * getOption("width"), symbol = ".") {
    line <- strrep(symbol, max(0L, width - 1L - nchar(text)))
    cat(text, " ", line, "\n", sep = "")
}

### Conditionally pluralizes most singular English nouns.
### Used to make sure 'print' and 'plot' method output is grammatical.
pluralize <- function(word, n) {
    plural <- (n > 1L)
    word[plural] <- paste0(word[plural], "s")
    word
}

### Formats tables of text.
### Used by 'print' methods instead of hacking 'print.default'.
align <- function(..., justify = "right", gap = 1L) {
    dots <- list(...) # list of column vectors
    dots <- Map(format, x = dots, justify = justify, USE.NAMES = FALSE)
    dots$sep <- strrep(" ", gap)
    do.call(paste, dots)
}

### Disambiguates duplicated strings (or names)
disambiguate <- function(x, nms = FALSE, fmt = "%s[%d]") {
    if (nms) {
        nx <- names(x)
        if (is.null(nx)) {
            return(x)
        }
        names(x) <- disambiguate(nx, nms = FALSE, fmt = fmt)
        return(x)
    }
    x <- as.character(x)
    f <- factor(x)
    n <- tabulate(f, nlevels(f))
    i <- unsplit(lapply(n, seq_len), f)
    sprintf(fmt, x, i)
}

### A drop-in replacement for 'rle' that regards NA as equal to previous NA
### and (in double vectors) NaN as equal to previous NaN
literal_rle <- function(x) {
    n <- length(x)
    if (n == 0L) {
        return(list(lengths = integer(0L), values = x))
    }
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

### Last observation carried forward, with optional replacement of leading NA
locf <- function(x, x0 = NULL) {
    if (!anyNA(x)) {
        return(x)
    }
    rle_x <- literal_rle(x)
    y <- rle_x$values
    if (is.na(y[1L]) && !is.null(x0)) {
        y[1L] <- x0
    }
    if (anyNA(y[-1L])) {
        argna_y <- which(c(FALSE, is.na(y[-1L])))
        y[argna_y] <- y[argna_y - 1L]
    }
    rle_x$values <- y
    inverse.rle(rle_x)
}

### Computes Wald confidence intervals from estimates and standard errors
wald <- function(estimate, se, level) {
    q <- qchisq(level, df = 1)
    n <- length(estimate)
    lu <- estimate + rep.int(sqrt(q) * c(-1, 1), c(n, n)) * se
    dim(lu) <- c(n, 2L)
    colnames(lu) <- c("lower", "upper")
    lu
}

cov2theta <- function(S) {
    n <- dim(S)[1L]
    log_sd <- 0.5 * log(diag(S, names = FALSE))
    R <- chol(S)
    R <- R * rep.int(1 / diag(R), rep.int(n, n))
    c(log_sd, R[upper.tri(R)])
}

theta2cov <- function(theta) {
    n <- as.integer(round(0.5 * (-1 + sqrt(1 + 8 * length(theta)))))
    i <- seq_len(n)
    R <- diag(n)
    R[upper.tri(R)] <- theta[-i]
    S <- crossprod(R)
    scale <- exp(theta[i] - 0.5 * log(diag(S, names = FALSE)))
    scale * S * rep.int(scale, rep.int(n, n))
}

##' Apply length-preserving functions to ragged vectors
##'
##' Modifies ragged vectors in place by replacing each group of elements
##' with the result of applying a length-preserving function to that group
##' of elements.
##'
##' @param x
##'   A vector or data frame.
##' @param index
##'   A factor (insofar as \code{as.factor(index)} is a factor)
##'   defining a grouping of the elements or rows of \code{x}.
##' @param f
##'   A function or list of one or more functions to be applied to the subsets
##'   of \code{x} defined by \code{index}.
##'   Functions are recycled to the number of levels of \code{index}.
##'   Each function must accept a vector argument matching \code{typeof(x)}
##'   (or \code{typeof(x[[j]])} for all \code{j} if \code{x} is a data frame)
##'   and return a vector of the same length.
##'
##' @details
##' Let \code{f} be a list of \code{nlevels(index)} functions,
##' and let \code{k = split(seq_along(index), index)}.
##' For vectors \code{x}, function \code{f[[i]]} is applied to
##' \code{x[k[[i]]]}.
##' For data frames \code{x}, function \code{f[[i]]} is applied to
##' \code{x[[j]][k[[i]]]} for all \code{j}.
##'
##' @return
##' If \code{x} is a vector, then a vector of the same length.
##' If \code{x} is a data frame, then a data frame with the same dimensions.

in_place_ragged_apply <- function(x, index, f) {
    if (!is.list(f)) {
        f <- list(f)
    }
    if (is.data.frame(x)) {
        do_call <- function(g, y) { y[] <- lapply(y, g); y }
    } else {
        do_call <- function(g, y) { g(y) }
    }
    split(x, index) <- Map(do_call, y = split(x, index), g = f)
    x
}
