#' a class representing a trajectory model used to fit exponential growth rates
#'
#' @slot expr an expression specifying the model
#' @slot grad gradient with respect to the parameters
#' @slot hessian Hessian of the model
#' @slot inits initial values
#' @slot input input variable name
#' @slot name short name of the model
#' @slot output output variable name
#' @slot par parameter names
#' @importFrom Deriv Deriv
#' @export
Model <- setClass(
  "Model",
  slots = c(
    name = "character",
    expr = "expression",
    output = "character",
    input = "character",
    par = "character",
    inits = "numeric",
    grad = "list",
    hessian = "list"
  )
)

tr2str <- function(x) {
    paste0(x@name,"=",deparse(x@expr[[1]]))
}

tmpf <- function(x,s) {
    if (length(slot(x,s))==0) "(none)" else paste(slot(x,s),collapse=", ")
}

#' @export
setMethod("show",
          "Model",
          definition = function(object) {
    cat(paste("model",sQuote(object@name)),":\n")
    cat("input:\n",tmpf(object,"input"),"\n")
    cat("parameters:\n",tmpf(object,"par"),"\n")
    if ("transforms" %in% slotNames(object)) {
        cat("transforms:\n",
            paste(sapply(object@transforms,tr2str),collapse=", "),"\n")
    }
})

#' second derivative of lbeta
#' use Taylor expansion of digamma(a+b) for a>>b
#' discontinuity in second derivative, but ... probably OK
#' @importFrom Deriv drule
dfun <- function(x,y,mag=1e8) {
    return(ifelse(x/y>mag,
                  -y*trigamma(x),
                  digamma(x)-digamma(x+y)))
}
dfun2 <- function(x,y,mag=1e8,focal="x") {
    return(switch(focal,
             x=ifelse(x/y>mag,
                  -y*psigamma(x,2),
                  trigamma(x)-trigamma(x+y)),
             y=ifelse(x/y>mag,
                  -trigamma(x),
                  -trigamma(x+y))))
}

w_lbeta <- function(a,b) {
    ## when we have an effectively-Poisson casee
    ## lbeta gives "underflow occurred in 'lgammacor'" frequently ...
    ## suppressWarnings() causes an obscure error ?
    ## using w_lbeta rather than lbeta causes obscure errors from Deriv()
    op <- options(warn=-1)
    on.exit(options(op))
    return(lbeta(a,b))
}

## drule[["NBconst"]] <- alist(x=-dfun(x,k)-1/x, ## blows up if X=0, but not used!
##                             k=-dfun(k,x))

## drule[["lbeta"]] <- drule[["w_lbeta"]] <- alist(a=dfun(a,b),
##                                                 b=dfun(b,a))

## drule[["dfun"]] <- alist(x=dfun2(x,y),
##                          y=dfun2(y,x))


my_deriv <- function(expr,vars) {
    d <- lapply(vars, function(p){Deriv(expr, p)})
    names(d) <- vars
    return(d)
}

#' the initializer for Model
#' @param model the formula specifying the model
#' @param input the input variable name
#' @param inits initial values
#' @param par the parameter names
#' @param keep_hessian maintain the hessian as part of the model?
#' @docType methods
#' @export Model
setMethod(
  "initialize",
  "Model",
  definition = function(.Object, name, model, input, par=NULL,
                        inits = rep(NA_real_,length(par)),
                        keep_hessian=FALSE) {
    ## using Deriv::Deriv, hessians are better, but still big
    ## and take too long to compute!
    .Object@name <- name
    if (!is(model, "formula"))
      stop("model must be a formula")
    f <- as.list(model)
    .Object@output <- as.character(f[[2]])
    .Object@expr <- as.expression(f[[3]])
    .Object@input <- input
    .Object@par <- par <- as.character(par)
    if (length(inits) <- length(par)) {
        oinits <- inits
        inits <- setNames(rep(NA_real_,length(par)),par)
        inits[names(oinits)] <- oinits
    }
    .Object@inits <- inits
    # compute the gradient
    vars <- c(input, par)
    .Object@grad <- my_deriv(.Object@expr,vars)
    if (keep_hessian) {
        .Object@hessian <- lapply(.Object@grad, my_deriv, vars=vars)
        names(.Object@hessian) <- vars
    } else .Object@hessian <- list() ## ugh, must be a list
    .Object
  }
)

#' Evaluate a model
#'
#' @docType methods
#' @param obj a model object
#' @rdname Eval
#' @export
setGeneric(
  "Eval",
  def = function(obj, ...) {
    standardGeneric("Eval")
  }
)

#' Evaluate a model
#'
#' @param data a dataframe object holding input values, if NULL, take from ...
#' @param par a named vector (or list) of parameter values, if NULL, take from ...
#' @param ... the input values and parameter values
#' @return numeric
#' @docType methods
#' @export
setMethod(
  "Eval",
  "Model",
  definition = function(obj, data=NULL, par = NULL, ...) {
      if (is.null(par)) browser()
      if (!is.null(par) && is.null(names(par))) {
          ## HACK to restore names when missing, e.g. within nloptr
          names(par) <- obj@par
      }
      frame <- c(as.list(par), as.list(data), list(...))
      eval(obj@expr, frame)
  }
)

#' @docType methods
#' @export
setGeneric(
  "grad",
  def = function(obj, ...) {
    standardGeneric("grad")
  }
)

#' compute gradients
#'
#' compute the gradients of the model negative log-likelihood as a function of the parameters
#' @param data a dataframe object holding inut values, if NULL, take from ...
#' @param par a named vector (or list) of parameters to compute the derivatives
#' @param ... the input values and parameter values
#' @return a dataframe with each column as a partial derivative values
#' @docType methods
#' @export
setMethod(
  "grad",
  "Model",
  definition <- function(obj, data=NULL, par, ...) {
      if (!is.null(par) && is.null(names(par))) {
          ## HACK to restore names when missing, e.g. within nloptr
          names(par) <- obj@par
      }
      par <- as.list(par)
      frame <- c(par, as.list(data), list(...))
      var <- names(par)
      l <- lapply(obj@grad[var], function(deriv) { eval(deriv, frame) })
      as.data.frame(l)
  }
)

#' @docType methods
#' @export
setGeneric(
  "hessian",
  def = function(obj, ...) {
    standardGeneric("hessian")
  }
)

#' the gradients w.r.t. to parameters
#' compute the gradient of the model as a function of the parameters
#' @param data a dataframe object holding inut values, if NULL, take from ...
#' @param par a named vector (or list) of parameters to compute the derivatives
#' @param ... the input values and parameter values
#' @return a dataframe with each column as a partial derivative values
#' @docType methods
#' @export
setMethod(
  "hessian",
  "Model",
  definition <- function(obj, data=NULL, par, ...) {
    par <- as.list(par)
    frame <- c(par, as.list(data), list(...))
    var <- names(par)
    l <- lapply(obj@hessian[var], function(dd) {
      sapply(dd[var], function(e) eval(e, frame))})
    n <- length(l)
    mat <- array(dim=c(dim(l[[1]]), n), dimnames=list(NULL, var, var))
    for (i in 1:n)
      mat[,,i] <- l[[i]]
    mat
  }
)

#' a model with transformed parameters
#'
#' @slot orig_expr original expression
#' @slot transforms the parameter transforms
#' @slot inverses the inverse transforms
#' @export
TransformedModel <- setClass(
  "TransformedModel",
  contains = "Model",
  slots = c(
    transforms = "list",
    inverses = "list",
    orig_expr = "expression"
  )
)

#' initializer for TransformedModel
#' @param model the model to be transformed
#' @param transforms the parameter transforms, a list of formula
#' @param inverses the inverse transforms, a list of formula
#' @importFrom methods callNextMethod is
#' @importFrom stats D
#' @export TransformedModel
#' @docType methods
setMethod(
  "initialize",
  "TransformedModel",
  definition = function(.Object, model, transforms=NULL, inverses=NULL) {
    # create models from transforms
    trans <- function(formulae, allvars) {
      # extract vars from an expression
      l <- list()
      for (f in formulae) {
        if (!is(f, "formula"))
          stop("transforms must be formula: ", as.character(f))
        var <- as.character(f[[2]]) # LHS
        if (!var %in% allvars) next
        input <- all.vars(f[[3]]) # extract vars from RHS
        l[[var]] <- new("Model", var, f, input = input)
      }
      l
    }
    # substitute the transformation expressions
    subst <- function(e) {
      if (is.numeric(e)) return(e)
      if (is.name(e)) {
        v <- as.character(e)
        expr <- transforms[[v]]
        if (is.null(expr)) return(e)
        return (expr@expr[[1]])
      }
      if (!is.call(e)) stop("unknown class: ", class(e))
      l <- list(e[[1]])
      for (i in 2:length(e))
        l <- c(l, subst(e[[i]]))
      as.call(l)
    }
    # if no transform, return model
    if (length(transforms) == 0)
      return(model)
    allvars <- c(model@input, model@par)
    .Object@orig_expr <- model@expr
    .Object@transforms <- transforms <- trans(transforms, allvars)
    .Object@inverses <- inverses <- trans(inverses, allvars)
    # recover the formula of the model
    f <- c(as.symbol("~"), as.symbol(model@output), subst(model@expr[[1]]))
    f <- as.formula(as.call(f))
    callNextMethod(.Object, model@name, f, model@input, model@par, model@inits)
  }
)

#' transform parameter values
#'
#' Map parameter values of a transformed model to/from the parameter values of the original model.
#' @param model a list of Model objects defining parameter transformations
#' @return a named list of transformed parameter values
#' @export
setGeneric(
  "transformPar",
  def = function(model, ...) {
    standardGeneric("transformPar")
  })

#' transform parameter values
#'
#' return the same parameter values as input, but in a named list
#'
#' @param model a Model object
#' @param par a named vector (or list) of parameters to be transformed
#' @param inverse logical, if TRUE, map from the TransformedModel parameters to the original model parameters; if FALSE, map from the original model parameters to the TransformedModel parameters
#' @return a named list of transformed parameter values
#' @export
setMethod(
  "transformPar",
  "Model",
  definition = function(model, par, inverse=FALSE) {
    as.list(par)
  }
)

#' transform parameter values
#'
#' Map parameter values of a transformed model to/from the parameter values of the original model.
#' @param trans a list of Model objects defining parameter transformations
#' @param par a named vector (or list) of parameters to be transformed
#' @param inverse logical, if TRUE, map from the TransformedModel parameters to the original model parameters; if FALSE, map from the original model parameters to the TransformedModel parameters
#' @return a named list of transformed parameter values
#' @export
setMethod(
  "transformPar",
  "TransformedModel",
  definition = function(model, par, inverse=FALSE) {
    trans <- if (inverse) model@inverses else model@transforms
    l <- lapply(names(par), function(v) {
        tr <- trans[[v]]
        if (is.null(tr)) par[[v]] else Eval(tr, par=par)
    })
    names(l) <- names(par)
    return(l)
}
)

## naming follows family()
link_defs <- list(log=list(linkfun=~log(x),
                           linkinv=~exp(x)),
                  atan=list(linkfun=~tan((x*2-1)*pi/2),
                            linkinv=~(atan(x)*2/pi+1)/2),
                  ## use explicit expressions rather than plogis/qlogis;
                  ##  not sure if we have Derivs package for extended
                  ##  derivatives table ...
                  logit=list(linkfun=~log(x/(1-x)),
                             linkinv=~1/(1+exp(-x))),
                  identity=list(linkfun=~x,
                                linkinv=~x))

## replace one symbol with another in a formula
## recursion seems to be a necessary evil here
rf <- function(x,from="x",to="r",debug=FALSE) {
    if (is.name(x) && x==as.symbol(from)) {
        if (debug) browser()
        return(as.symbol(to))
    }
    if (length(x)>1) {
        for (i in seq_along(x)) {
            x[[i]] <- rf(x[[i]],from,to)
        }
    }
    return(x)
}

## return a link/inverse link formula for a variable
makelink <- function(link="log",v="r",which="linkfun") {
    res <- rf(link_defs[[link]][[which]],"x",v)
    res[[3]] <- res[[2]]
    res[[2]] <- as.symbol(v)
    return(res)
}

## FIXME: implement proper naming (r -> log.r, log.r -> r)
maketrans <- function(model,links, trans=list(), itrans=list()) {
    for (i in seq_along(links)) {
        L <- links[[i]]
        v <- names(links)[i]
        trans  <- append(trans,makelink(links[[i]],names(links)[i],"linkinv"))
        itrans <- append(itrans,makelink(links[[i]],names(links)[i],"linkfun"))
    }
    return(TransformedModel(model,
                     transforms=trans,
                     inverses=itrans))
}


#' add a constant incidence term (b*t) to RHS of formula
#' @param cform base formula
#' @keywords internal
add_const_incidence <- function(cform) {
    cform[[3]] <- substitute(x+y,list(x=quote(b*t),
                                      y=cform[[3]]))
    return(cform)
}


#' return a model object (deterministic/mean trajectory) for fitting
#' @param model base model name ("exp","logistic","richards")
#' @param add_const (logical) add a constant incidence term?
#' @param s_lower (numeric) for \code{model=="richards"} only
#' (otherwise ignored): add lower bound above zero for shape (s)
#' parameter?
#' @param link_vals link functions
#' @export
get_model <- function(model,
                      add_const=FALSE,
                      s_lower=0.1,
                      link_vals=NULL) {

    default_link_vals <- c(r="log",
                           x0=if (model=="exp") "log" else "atan",
                           b="log", K="log", s="log")

    ## FIXME: modify model name/make it easier to determine which
    ##  variants (constant incidence term, s bounds) were added?

    ## "link" here means the transformation between the natural
    ## parameter and the (unconstrained) scale on which we fit the
    ## model.  Unconstrained scales are better because we can't fit
    ## past an impossible boundary; also, this usually means the
    ## parameters are scaled similarly.
    if (is.null(link_vals)) link_vals <- default_link_vals
    for (i in setdiff(names(default_link_vals), names(link_vals))) {
        link_vals[[i]] <- default_link_vals[[i]]
    }

    ## build model piecewise
    parinfo <- list(r=list(link=link_vals[["r"]], init=NA_real_))
    if (model == "exp") {
        form <- X ~ x0*exp(r * t)
        parinfo <- append(parinfo, list(x0=list(link=link_vals[["x0"]],
                                                init=NA_real_)))
    } else {
        ## logistic or Richards: add carrying capacity
        ##  x0 is scaled to carrying capacity (!!)
        ##  e.g. see machinery in epigrowthfit initialize method
        ##  (currently l. 392 in epigrowthfit.R); scales estimated init(x0)
        ##  to init(x0)/init(K)
        ##     hence we fit it on (0,1) scale:
        parinfo <- append(parinfo, list(x0=list(link=link_vals[["x0"]],init=NA_real_),
                                        K=list(link=link_vals[["K"]],init=NA_real_)))
        if (model =="logistic") {
            form <- X ~ K / (1 + (1/x0-1) * exp(-r*t))
        } else {
            ## Richards: add shape parameter
            parinfo <- append(parinfo, list(s=list(link=link_vals[["s"]],
                                                              init=1.001)))
            form <- X ~ K / (1 + (1/x0-1) * exp(-r*s*t))^(1/s)
        }
    }
    ## add constant incidence if requested
    if (add_const) {
        parinfo <- append(parinfo, list(b=list(link=link_vals[["b"]],init=0.001)))
        ## form <- update(form,.~.+b*t) ## doesn't work because of simplification
        form <- add_const_incidence(form)
    }
    links <- sapply(parinfo,"[[","link")
    mod <- new("Model", model,
               form, input="t",
               inits=sapply(parinfo,"[[","init"),
               par=names(links))
    trans <- itrans <- list()
    ## if we impose a lower bound on s then we can't just use a log link:
    if (model=="richards" && s_lower != 0) {
        ## add fancy transforms for s; subtract from list of links
        trans <- as.formula(substitute(s ~ exp(s)+s_lower,list("s_lower"=s_lower)))
        itrans <- as.formula(substitute(s ~ log(s-s_lower),list("s_lower"=s_lower)))
        links <- links[names(links) != "s"]
    }

    mod <- maketrans(mod,as.list(links), trans, itrans)

    return(mod)
}

#' @export
NBconst <- function(k,x) {
    return(ifelse(x==0,0,-lbeta(k,x)-log(x)))
}

#' @export
get_loglik <- function(loglik) {
    switch(loglik,
           poisson = {
        new("Model", "poisson",
            LL ~ -lgamma(X+1) + X*log(lambda) - lambda,
            input = c("lambda", "X"),
            par = c())
    },
    nbinom = {
        nb <- new ("Model", "nbinom",
                   LL ~ NBconst(ll.k,X) +
                       ll.k*(-log1p(mu/ll.k)) +
                       X * log(mu) - X * log(ll.k + mu),
                   input = c("mu", "X"),
                   par = "ll.k",
                   inits = c(ll.k=1))

        TransformedModel(
            nb,
            transforms = list(ll.k ~ exp(ll.k)),
            inverses = list(ll.k ~ log(ll.k))
        )
    },
    nbinom1 = {
        ## negative binomial '1' likelihood
        ## var proportional to mean
        ## v=mu*(1+mu/k), k>0
        ## v=mu*(1+phi), phi>0
        ## i.e. mu/k=phi -> k=mu/phi
        nb1 <- new("Model", "nbinom1",
                   LL ~ lgamma(mu/ll.phi+X) - lgamma(mu/ll.phi) - lgamma(X+1) +
                       mu/ll.phi*log(mu/ll.phi) - mu/ll.phi*log(mu/ll.phi+mu) + X*log(mu) -
                       X*log(mu/ll.phi+mu),
                   input = c("mu", "X"),
                   par = "ll.phi",
                   inits = c(ll.phi=1))
        TransformedModel(
            nb1,
            transforms = list(ll.phi ~ exp(ll.phi)),
            inverses = list(ll.phi ~ log(ll.phi))
        )
    },
    stop("unknown model type ",loglik,
         "(poisson, nbinom, nbinom1 are current choices)")
    ) ## end switch
}


#' set parameter values for a trajectory model
#'
#' @description
#' transformed parameter values are returned
#'
#' @param model an object of class \code{Model}
#' @param par a named parameter vector (or \code{NULL} for default values)
#' @export
setModelPars <- function(model,par=NULL) {
  r <- K <- x0 <- s <- NULL ## for R CMD check. Don't do this after assigning!!
  ## default parameter values:
  x0.default <- 1e-2
  r.default <- 1
  K.default <- 1
  s.default <- 1+1e-10 # need s>1 to get good fits
  ## set pars that were listed in arg list:
  for (p in names(par)) assign(p,as.numeric(par[p]))
  ## for any pars not passed as args, set to default value:
  still.to.set <- setdiff(model@par, names(par))
  ## FIXME: ugh ... can we do this without assign/get?
  for (p in still.to.set) assign(p,get(paste0(p,".default")))
  par <- switch(model@name,
                exp = c(r=r, x0=log(x0)), # meaning of x0 different for exp!
                logistic = c(r=r, x0=x0, K=K),
                richards = c(r=r, x0=x0, K=K, s=s)
                )
  return(transformPar(model, par, inverse=TRUE))
}

#' plot a model used for fitting exponential growth rates
#' @param x a \code{Model} object
#' @param par a named vector or list of parameter values associated
#'        with the model.  Default values will be chosen for any that
#'        are not specified.
#' @param times vector of times at which to plot the model (currently,
#'        only \code{range(times)} affects the resulting plot)
#' @param add logical: if \code{TRUE} then add to existing plot
#' @param cumulative logical: if \code{TRUE} show cumulative curve
#' @param ... further arguments to be passed to \code{\link{plot}}
#'        and \code{\link{curve}}
#' @docType methods
#' @export
setMethod("plot",
          signature(x="Model",y="missing"),
          definition=function(x,par=NULL,times=NULL,add=FALSE,
                      cumulative=TRUE, ...) {
  ## set any unset model pars to default values and return transformed values:
  par <- setModelPars(x,par=par)
  ## spew transformed and original parameter values:
  cat("Transformed parameter values:\n")
  print(unlist(par))
  cat("Original parameter values:\n")
  print(unlist(transformPar(x, par)))
  ## set default times if not specified in arg
  if (is.null(times)) times <- seq(0,15,length=10)
  ## grad returns a dataframe, whereas Eval returns a vector
  if (cumulative) {
    f <- function(t) {Eval(x, par=par, t=t)}
  } else {
    f <- function(t) {grad(x, par=list(t=t), data = as.list(par))$t}
  }
  tt <- times
  y <- f(tt)
  if (!add) plot(tt,y,type="n",ylim=c(0,max(y)),xlab="Time",ylab="Model",...)
  curve(f(x),add=TRUE,...)
})
