#include <TMB.hpp>
#include <iostream>
#include <string>
#include "curve_enum.h"
#include "distr_enum.h"

// eventual goal: flexible epigrowthfit-style fitting
// possibility of including multiple time series with random effects parameters
// compare cumulative vs incidence fitting?
// regularization? priors?

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// https://kaskr.github.io/adcomp/_book/Tutorial.html
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(t);      // times (includes step [-1])
  DATA_VECTOR(x);      // counts
  DATA_INTEGER(debug); // debug flag (1=print debugging statements)
  DATA_INTEGER(curve_flag);  // cumulative curve model
  DATA_INTEGER(distr_flag);  // response distribution model
  // DATA_INTEGER(baseline_flag);  // 0 or 1 depending on baseline or not
  
  // Parameters  
  PARAMETER(log_thalf);    // log half-max time
  PARAMETER(log_K);        // log carrying capacity
  PARAMETER(log_r);        // log growth rate
  PARAMETER(log_p);        // log Richards shape parameter
  PARAMETER(log_x0);       // log initial value
  // PARAMETER(log_b);     // log baseline value
  PARAMETER(log_nbdisp);   // *inverse* NB dispersion parameter

  // FIXME: add background-count model, parameters?

  vector<Type> cumcurve(t.size()); // includes starting time (1 longer than x)
  vector<Type> inccurve(x.size());
  vector<Type> log_inccurve(x.size()); // redundant but useful for plots
  // FIXME: implement dpois with log_lambda specified?
  
  // Objective function
  Type jnll = 0;

  // inverse-link transform parameters
  Type r = exp(log_r);
  Type K = exp(log_K);         // not for exponential
  Type thalf = exp(log_thalf); // not for exponential
  Type x0 = exp(log_x0);       // only for exponential
  Type p = exp(log_p);         // only for richards

  Type log_nbexcessvar; // NB excess variance; only used for NB model
  
  // parameter est

  // std::cout << "thalf " << thalf << " K " << K << " r " << r << " p " << p << " nbdisp " << nb_disp << std::endl;

  
  // first calculate cumulative curve at observed points
  for (int i = 0; i < t.size(); i++) {
	  switch(curve_flag) {
	  case exponential:
		  cumcurve(i) = x0*exp(r*(t(i)));
		  break;
	  case logistic:
		  cumcurve(i) = K/(1+exp(-r*(t(i)-thalf)));
		  break;
	  case richards:
		  cumcurve(i) = K/pow(1+exp(-r*p*(t(i)-thalf)),1/p);
		  // std::cout << cumcurve(i) << std::endl;
		  break;
	  default:
		  std::cout << curve_flag << std::endl;
		  error("curve not implemented");
	  }
  }

  // difference cumulative curve to get expected incidence
  for (int i = 0; i < x.size(); i++) {
	  if (!isNA(x(i))) {
		  switch(distr_flag) {
		  case nbinom2:
			  log_inccurve(i) = log(cumcurve(i+1)-cumcurve(i));  // compute on *LOG SCALE*
			  log_nbexcessvar = log_inccurve(i)-log_nbdisp; // compute NB *excess* variance from mean, dispersion param
			  jnll -= dnbinom_robust(x(i),log_inccurve(i),log_nbexcessvar,1); // 1 = log
			  break;
		  case poisson:
			  inccurve(i) = cumcurve(i+1)-cumcurve(i);
			  log_inccurve(i) = inccurve(i); // for reporting only
			  jnll -= dpois(x(i),inccurve(i),1);
			  break;
		  default:
			  error("distribution not implemented");
		  }
	    if (debug==1) {
		    std::cout << i << " " << x(i) << " " << cumcurve(i) << " " <<
			    log_inccurve(i)  << " " << log_nbexcessvar << " " << jnll << std::endl;
	    }
	  }
  } // loop over observations

  // 
  REPORT(log_inccurve);
  ADREPORT(log_inccurve);
  return jnll;
}
