#' @importFrom TMB MakeADFun sdreport
egf <- function(formula,
                fixed = NULL,
                random = NULL,
                data = parent.frame(),
                index = NULL,
                curve = c("richards", "logistic", "exponential"),
                distr = c("nbinom", "pois"),
                excess = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "CG"),
                na_action = c("exclude", "fail"),
                sparse_X = FALSE,
                sparse_Z = TRUE,
                date_format = "%Y-%m-%d",
                ...) {
  curve <- match.arg(curve)
  distr <- match.arg(distr)
  stop_if_not(
    inherits(excess, "logical"),
    length(excess) == 1L,
    !is.na(excess),
    m = "`excess` must be TRUE or FALSE."
  )
  na_action <- match.arg(na_action)
  method <- match.arg(method)

  par_names <- get_par_names(curve, distr, excess)
  formula <- check_formula(formula)
  fixed <- check_fixed(fixed, par_names)
  random <- check_random(random, par_names)
  frame <- check_data(formula, fixed, random, data, index,
                      na_action, date_format)

  madf_data <- make_madf_data(frame, sparse_X, sparse_Z)
  madf_args <- list(
    data = madf_data,
    parameters = make_madf_parameters(madf_data),
    map = make_madf_map(madf_data),
    random = make_madf_random(madf_data),
    DLL = "epigrowthfit",
    silent = TRUE
  )
  madf_out <- do.call(MakeADFun, madf_args)

  if (method == "nlminb") {
    optim_out <- nlminb(
      start = madf_out$par,
      objective = madf_out$fn,
      gradient = madf_out$gr,
      ...
    )
    nll <- optim_out$objective
  } else if (method == "nlm") {
    optim_out <- nlm(
      f = structure(madf_out$fn, gradient = madf_out$gr),
      p = madf_out$par,
      ...
    )
    nll <- optim_out$minimum
  } else {
    optim_out <- optim(
      par = madf_out$par,
      fn = madf_out$fn,
      gr = madf_out$gr,
      method = method,
      ...
    )
    nll <- optim_out$value
  }

  par_init <- madf_out$env$par
  par <- madf_out$env$last.par.best
  par_info <- get_par_info(par, madf_data, decontr = FALSE)
  inr <- which(par_info$name != "b")

  nll_func <- evalq(
    expr = function(x = par) as.numeric(madf_out$fn(x[inr])),
    envir = list(par = par, inr = inr, fn = madf_out$fn),
    enclos = baseenv()
  )
  nll_grad <- evalq(
    expr = function(x = par) as.vector(madf_out$gr(x[inr])),
    envir = list(par = par, inr = inr, gr = madf_out$gr),
    enclos = baseenv()
  )

  out <- list(
    frame = frame,
    index = attr(frame, "index"),
    curve = curve,
    distr = distr,
    excess = excess,
    madf_args = madf_args,
    madf_out = madf_out,
    madf_report = madf_out$report(par),
    madf_sdreport = sdreport(madf_out),
    optim_out = optim_out,
    par_init = par_init,
    par = par,
    par_info = par_info,
    inr = inr,
    nll = nll,
    nll_func = nll_func,
    nll_grad = nll_grad,
    call = match.call()
  )
  class(out) <- c("egf", "list")
  out
}
