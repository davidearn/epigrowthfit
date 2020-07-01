#include <TMB.hpp>
#include <iostream>
#include <string>

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
  
  // Parameters  
  PARAMETER(log_thalf);    // log half-max time
  PARAMETER(log_K);        // log carrying capacity
  PARAMETER(log_r);        // log growth rate
  PARAMETER(log_p);        // log Richards shape parameter
  PARAMETER(log_nb_disp);  // log NB dispersion parameter

  vector<Type> cumcurve(t.size()); // includes starting time
  vector<Type> inccurve(x.size());
  
  // Objective function
  Type jnll = 0;

  // inverse-link transform parameters
  Type thalf = exp(log_thalf);
  Type K = exp(log_K);
  Type r = exp(log_r);
  Type p = exp(log_p);
  Type nb_disp = exp(log_nb_disp); 

  Type nb_var; // NB variance
  
  // parameter est

  // std::cout << "thalf " << thalf << " K " << K << " r " << r << " p " << p << " nbdisp " << nb_disp << std::endl;

  // first calculate cumulative curve at observed points
  for (int i = 0; i < t.size(); i++) {
    cumcurve(i) = K/pow(1+exp(-r*p*(t(i)-thalf)),1/p);  // Richards
  }

  // difference cumulative curve to get expected incidence
  for (int i = 0; i < x.size()-1; i++) {
    inccurve(i) = cumcurve(i+1)-cumcurve(i);
    if (!isNA(x(i))) {
      nb_var = inccurve(i)*(1+inccurve(i)/nb_disp); // compute NB variance from mean, dispersion param
      jnll -= dnbinom2(x(i),inccurve(i),nb_var,1);
      if (debug==1) {
	std::cout << i << " " << x(i) << " " << cumcurve(i) << " " <<
	  inccurve(i)  << " " << nb_var << " " << jnll << std::endl;
      }
    }
  }
  
  REPORT(inccurve);
  ADREPORT(inccurve);
  return jnll;
}
