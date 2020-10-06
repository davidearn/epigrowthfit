#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "curve_enum.h"
#include "distr_enum.h"

/* TODO:
- fit multiple time series with random effects (latent variables)
- compare fits to cumulative/interval incidence
- regularization? priors? 
- better alignment of printed text with printf
- compute everything on log scale and use logspace_sub for differencing
  - no so great when baseline mortality is included
*/

/* Wrap body of second loop in `if (!isNA(x(i))) { }` 
   if you decide to remove check for NA on R side
   https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA(Type x)
{
    return R_IsNA(asDouble(x));
}
*/

/* Use if you decide to compute log cumulative incidence
template<class Type>
Type dpois_robust(Type x, Type log_lambda)
{
    return x * log_lambda - exp(log_lambda) - lfactorial(x);
}
*/

// https://kaskr.github.io/adcomp/Introduction.html
template<class Type>
Type objective_function<Type>::operator() ()
{
    /*SET UP =================================================================*/

    // Data
    DATA_VECTOR(t);              // time,  length N+1
    DATA_VECTOR(x);              // cases, length N
    DATA_INTEGER(curve_flag);    // cum. inc. model     (curve_enum.h)
    DATA_INTEGER(baseline_flag); // baseline growth     (1=include)
    DATA_INTEGER(distr_flag);    // observation model   (distr_enum.h)

    // Parameters
    PARAMETER(log_r);      // log initial growth rate
    PARAMETER(log_x0);     // log initial cases         (exponential)
    PARAMETER(log_K);      // log final size            (logistic, richards)
    PARAMETER(log_thalf);  // log time half final size  (logistic, richards)
    PARAMETER(log_p);      // log Richards shape        (richards)
    PARAMETER(log_b);      // log baseline growth rate  (baseline_flag)
    PARAMETER(log_nbdisp); // log nb dispersion         (nbinom)

    // Inverse-link transform parameters
    Type r = exp(log_r);
    Type x0 = exp(log_x0);
    Type K = exp(log_K);
    Type thalf = exp(log_thalf);
    Type p = exp(log_p);
    Type b = exp(log_b);
    Type nbdisp = exp(log_nbdisp);
    Type log_nbxsvar;
    
    // Objective function
    Type nll = Type(0);
    // Fitted curves
    vector<Type> cum_inc(t.size());     // cumulative incidence,    length N+1
    vector<Type> int_inc(x.size());     // interval incidence,      length N
    vector<Type> log_int_inc(x.size()); /* log interval incidence,  length N
					   - redundant but useful for plots
					   - though plots are confusing if 
					     time points are unequally spaced */

    
    /*COMPUTE NEGATIVE LOG LIKELIHOOD ========================================*/

    // Evaluate cumulative incidence at observation times
    for (int i = 0; i < t.size(); i++)
    {
        switch (curve_flag)
        {
        case exponential:
            cum_inc(i) = x0 * exp(r * t(i));
            break;
        case logistic:
            cum_inc(i) = K / (Type(1) + exp(-r * (t(i) - thalf)));
            break;
        case richards:
            cum_inc(i) = K / pow(Type(1) + (pow(Type(2), p) - Type(1)) *
                exp(-r * p * (t(i) - thalf)), 1 / p);
            break;
        }

        if (baseline_flag == 1)
        {
            cum_inc(i) += b * t(i);
        }
    }

    // Compute interval incidence and increment negative log likelihood
    for (int i = 0; i < x.size(); i++)
    {
        int_inc(i) = cum_inc(i + 1) - cum_inc(i);
        log_int_inc(i) = log(int_inc(i));

        switch (distr_flag)
        {
        case pois:
            nll -= dpois(x(i), int_inc(i), true);
            break;
        case nbinom:
            log_nbxsvar = Type(2) * log_int_inc(i) - log_nbdisp;
            nll -= dnbinom_robust(x(i), log_int_inc(i), log_nbxsvar, true);
            break;
	}
    }

    REPORT(log_int_inc);
    ADREPORT(log_int_inc);
    return nll;   
}
