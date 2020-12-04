#' @importFrom TMB MakeADFun
egf <- function(formula,
                data = parent.frame(),
                index = NULL,
                fixed = NULL,
                random = NULL,
                spX = FALSE,
                spZ = FALSE,
                curve = c("richards", "logistic", "exponential"),
                distr = c("nbinom", "pois"),
                include_baseline = FALSE,
                na_action = c("exclude", "fail"),
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "CG"),
                dfmt = "%Y-%m-%d",
                ...) {
  curve <- match.arg(curve)
  distr <- match.arg(distr)
  check(include_baseline,
    what = "logical",
    len = 1L,
    no = is.na,
    "`include_baseline` must be TRUE or FALSE."
  )
  na_action <- match.arg(na_action)
  method <- match.arg(method)

  par_names <- get_par_names(curve, distr, include_baseline)
  formula <- check_formula(formula)
  fixed <- check_fixed(fixed, par_names)
  random <- check_random(random, par_names)
  frame <- check_data(formula, data, index, fixed, random, dfmt, na_action)

  madf_data <- make_madf_data(frame, spX, spZ)
  madf_parameters <- make_madf_parameters(madf_data)
  madf_map <- make_madf_map(madf_data)
  madf_random <- make_madf_random(madf_data)

  madf_out <- MakeADFun(
    data = madf_data,
    parameters = madf_parameters,
    map = madf_map,
    random = madf_random,
    DLL = "epigrowthfit",
    silent = TRUE
  )

  if (method == "nlminb") {
    optim_out <- nlminb(
      start = madf_out$par,
      objective = madf_out$fn,
      gradient = madf_out$gr,
      ...
    )
    par <- optim_out$par
    nll <- optim_out$objective
  } else if (method == "nlm") {
    optim_out <- nlm(
      f = structure(madf_out$fn, gradient = madf_out$gr),
      p = madf_out$par,
      ...
    )
    par <- optim_out$estimate
    nll <- optim_out$minimum
  } else {
    optim_out <- optim(
      par = madf_out$par,
      fn = madf_out$fn,
      gr = madf_out$gr,
      method = method,
      ...
    )
    par <- optim_out$par
    nll <- optim_out$value
  }

  out <- list(
    frame = frame,
    index = attr(frame, "index"),
    curve = curve,
    distr = distr,
    include_baseline = include_baseline,
    par0 = madf_out$par,
    madf_data = madf_data,
    madf_out = madf_out,
    optim_out = optim_out,
    par = par,
    nll = nll,
    nll_func = function(par = par) madf_out$fn(par),
    nll_grad = function(par = par) as.vector(madf_out$gr(par)),
    call = match.call()
  )
  structure(out, class = c("egf", "list"))
}
