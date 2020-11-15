#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "curve_enum.h"
#include "distr_enum.h"

/* TODO:
- fit multiple time series with mixed effects (latent variables)
- compare fits to cumulative/interval incidence
- regularization? priors? 
*/

template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA(Type x)
{
    return R_IsNA(asDouble(x));
}

// https://kaskr.github.io/adcomp/Introduction.html
template<class Type>
Type objective_function<Type>::operator() ()
{
    // SET UP ==================================================================

    // Data
    DATA_VECTOR(t);              // time,  length N+1
    DATA_VECTOR(x);              // cases, length N
    DATA_INTEGER(curve_flag);    // cum. inc. model    (curve_enum.h)
    DATA_INTEGER(distr_flag);    // observation model  (distr_enum.h)
    DATA_INTEGER(baseline_flag); // baseline growth    (1=include)

    // Parameters
    PARAMETER(log_r);      // log initial growth rate
    PARAMETER(log_c0);     // log initial cum. inc.     (exponential)
    PARAMETER(log_K);      // log final size            (logistic, richards)
    PARAMETER(log_thalf);  // log time half final size  (logistic, richards)
    PARAMETER(log_p);      // log Richards shape        (richards)
    PARAMETER(log_nbdisp); // log nb dispersion         (nbinom)
    PARAMETER(log_b);      // log baseline growth rate  (baseline_flag)
    
    // Inverse link-transformed parameters
    Type r = exp(log_r);
    Type thalf = exp(log_thalf);
    Type p = exp(log_p);

    // Other stuff
    Type log_nbxsvar; // log nb excess variance (nbinom)
    
    // Fitted curves
    vector<Type> log_cum_inc(t.size()); // log cumulative incidence, length N+1
    vector<Type> log_int_inc(x.size()); // log interval incidence,   length N

    // Objective function
    Type nll = Type(0);

    
    // COMPUTE NEGATIVE LOG LIKELIHOOD =========================================

    // Evaluate log cumulative incidence at observation times
    for (int i = 0; i < t.size(); i++)
    {
        switch (curve_flag)
        {
        case exponential:
	    log_cum_inc(i) = log_c0 * r * t(i);
	    break;
        case logistic:
	    log_cum_inc(i) = log_K - logspace_add(Type(0), -r * (t(i) - thalf));
	    break;
        case richards:
	    log_cum_inc(i) = log_K - logspace_add(Type(0), logspace_sub(p * log(2), Type(0)) - r * p * (t(i) - thalf)) / p;
	    break;
        }

        if (baseline_flag == 1)
        {
	    log_cum_inc(i) = logspace_add(log_b + log(t(i)), log_cum_inc(i));
        }
    }

    // Compute log interval incidence and increment negative log likelihood
    for (int i = 0; i < x.size(); i++)
    {
        if (!isNA(x(i)))
        {
	    log_int_inc(i) = logspace_sub(log_cum_inc(i+1), log_cum_inc(i));

            switch (distr_flag)
            {
            case pois:
	        nll -= dpois_robust(x(i), log_int_inc(i), true);
		break;
            case nbinom:
	        log_nbxsvar = Type(2) * log_int_inc(i) - log_nbdisp;
		nll -= dnbinom_robust(x(i), log_int_inc(i), log_nbxsvar, true);
		break;
	    }
        }
    }

    ADREPORT(log_cum_inc)
    ADREPORT(log_int_inc)
    return nll;   
}
