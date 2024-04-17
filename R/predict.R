predict.egf <-
function(object, newdata = NULL, class = FALSE, se = FALSE, ...) {
	if (is.null(newdata)) {
		ans <- fitted(object, class = class, se = se)
		oldClass(ans)[oldClass(ans) == "fitted.egf"] <- "predict.egf"
		return(ans)
	}
	.NotYetImplemented() # TODO
}

confint.predict.egf <- confint.fitted.egf
