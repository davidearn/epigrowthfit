#include <TMB.hpp>

template<class Type>
Type add_terms(const vector<Type> &y,
	       const vector<Type> &x,
	       Type a,
	       Type b,
	       Type log_sigma)
{
    Type res = Type(0.0);
    Type sigma = exp(log_sigma);
    for (int i = 0; i < y.size(); ++i)
    {
        res -= dnorm(y(i), a + b * x(i), sigma, true);
    }
    return res;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
#ifdef _OPENMP
    this->max_parallel_regions = omp_get_max_threads();
#endif
    
    DATA_VECTOR(y);
    DATA_VECTOR(x);
    PARAMETER(a);
    PARAMETER(b);
    PARAMETER(log_sigma);
    
    Type nll = Type(0.0);
    PARALLEL_REGION
    {
        nll += add_terms(y, x, a, b, log_sigma);
    }
    return nll;
}
