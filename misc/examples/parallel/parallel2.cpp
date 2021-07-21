#include <TMB.hpp>

template<class Type>
void add_terms(Type &nll,
	       objective_function<Type> *obj,
	       const vector<Type> &y,
	       const vector<Type> &x,
	       Type a,
	       Type b,
	       Type log_sigma)
{
    Type sigma = exp(log_sigma);
    for (int i = 0; i < y.size(); ++i)
    {
        if (obj->parallel_region())
	{
	    nll -= dnorm(y(i), a + b * x(i), sigma, true);
	}
    }
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
    add_terms(nll, this, y, x, a, b, log_sigma);    
    return nll;
}
