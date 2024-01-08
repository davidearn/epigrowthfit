\docType{package}
\name{epigrowthfit-package}
\alias{epigrowthfit-package}
\title{epigrowthfit: Nonlinear Mixed Effects Models of Epidemic Growth}
\description{
Maximum likelihood estimation of nonlinear mixed effects models
of epidemic growth using Template Model Builder (\pkg{TMB}).
}
\details{
At the top level, \pkg{epigrowthfit} fits phenomenological nonlinear
models of population growth to segments of one or more epidemic time
series. At a lower level, it jointly estimates generalized linear
mixed effects models of all nonlinear and dispersion model parameters.
This framework supports straightforward estimation of initial epidemic
growth rates, while enabling broader characterization of the effects
of covariates of interest.
}
\references{
Kristensen K, Nielsen A, Berg CW, Skaug H, Bell BM.
TMB: Automatic differentiation and Laplace approximation.
J Stat Softw. 2016;70:1--21.

Ma J, Dushoff J, Bolker BM, Earn DJD.
Estimating initial epidemic growth rates.
Bull Math Biol. 2014;76:246--60.

Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM.
Acceleration of plague outbreaks in the second pandemic.
Proc Natl Acad Sci U S A. 2020;117:27703--11.
}
\keyword{internal}
