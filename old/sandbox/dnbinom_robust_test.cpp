#include <TMB.hpp>
#include <iostream>
#include <string>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(x);      // counts
  
  // Parameters
  PARAMETER(log_mu);       // log mean
  PARAMETER(log_nbdisp);  // log NB dispersion parameter

  // Objective function
  Type jnll = 0;
  Type log_nb_excess_var = 2*log_mu-log_nbdisp;
  
  jnll -= sum(dnbinom_robust(x,log_mu,log_nb_excess_var,1)); // 1 = log

  return jnll;
}
