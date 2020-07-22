#include <TMB.hpp>
#include <iostream>
#include <string>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(x);      // counts
  
  // Parameters
  PARAMETER(log_mu);      // log mean
  PARAMETER(inv_nbdisp);  // inv NB dispersion parameter

  // Objective function
  Type jnll = 0;
  Type mu = exp(log_mu);
  Type nb_var = mu*(1+mu*inv_nbdisp);

  jnll -= sum(dnbinom2(x,mu,nb_var,1)); // 1 = log

  return jnll;
}
