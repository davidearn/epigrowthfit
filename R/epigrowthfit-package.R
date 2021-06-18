#' epigrowthfit: Nonlinear Mixed Effects Models of Epidemic Growth
#'
#' A description of the \pkg{epigrowthfit} package.
#'
#' @details
#' At the top level, \pkg{epigrowthfit} fits phenomenological nonlinear
#' models of population growth to segments of one or more epidemic time
#' series. At a lower level, it jointly estimates generalized linear
#' mixed effects models of all nonlinear and dispersion model parameters.
#' This framework supports straightforward estimation of initial epidemic
#' growth rates, while enabling broader characterization of the effects
#' of covariates of interest.
#'
#' Models are fit by maximum likelihood estimation using Template Model
#' Builder (\pkg{TMB}).
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @name epigrowthfit-package
#' @docType package
#' @keywords internal
#' @importFrom Rdpack reprompt
NULL

.epigrowthfit <- new.env(parent = emptyenv())
