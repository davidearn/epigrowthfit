namespace egf
{

/* Similar to R function `is.na` */
template<class Type>
bool is_na(Type x)
{
    return R_IsNA(asDouble(x));
}

/* Similar to R function `is.finite` */
template<class Type>
bool is_finite(Type x)
{
    return R_finite(asDouble(x));
}

/* Poisson density with robust parametrization */
template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

/* Negative binomial sampling with robust parametrization */
template<class Type>
Type rnbinom_robust(Type log_mu, Type log_size)
{
    Type log_prob = log_size - logspace_add(log_mu, log_size);
    return rnbinom(exp(log_size), exp(log_prob));
    /* usage: rnbinom(size, prob) */
}

/* Compute `log(diff(x))` given `log(x)` */
template<class Type>
vector<Type> logspace_diff(vector<Type> log_x)
{
    int n = log_x.size() - 1;
    for (int i = 0; i < n; ++i)
    {
        log_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    log_x.conservativeResize(n);
    return log_x;
}

template<class Type>
void logspace_diff(vector<Type> *log_x)
{
    *log_x = logspace_diff(*log_x);
}

/* Compute `log(cumsum(x))` given `log(x)` */
template<class Type>
vector<Type> logspace_cumsum(vector<Type> log_x)
{
    int n = log_x.size();
    for (int i = 1; i < n; ++i)
    {
        log_x(i) = logspace_add(log_x(i-1), log_x(i));
    }
    return log_x;
}

template<class Type>
void logspace_cumsum(vector<Type> *log_x)
{
    *log_x = logspace_diff(*log_x);
}

} // namespace egf
