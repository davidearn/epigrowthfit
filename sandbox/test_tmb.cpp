#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(log_sigma);
  Type nll = -sum(dnorm(y, a + b * x, exp(log_sigma), true));

  // Test code below

  

  return nll;
}
