template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

template<class Type>
Type rnbinom_robust(Type log_mu, Type log_disp)
{
    Type log_a = log_disp - logspace_add(log_mu, log_disp);
    return rnbinom(exp(log_disp), exp(log_a));
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA_real_(Type x)
{
    return R_IsNA(asDouble(x));
}
