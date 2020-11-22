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
    DATA_VECTOR(t); // time,  length N+1
    DATA_VECTOR(x); // cases, length N

    // Flags
    DATA_INTEGER(curve_flag);    // cum. inc. model    (curve_enum.h)
    DATA_INTEGER(distr_flag);    // observation model  (distr_enum.h)
    DATA_INTEGER(baseline_flag); // baseline growth    (1=include)
    DATA_INTEGER(predict_flag);  // predict or fit     (1=predict)

    // Parameters
    PARAMETER(log_r);      // log initial growth rate
    PARAMETER(log_c0);     // log initial cum. inc.     (exponential)
    PARAMETER(log_thalf);  // log time half final size  (logistic, richards)
    PARAMETER(log_K);      // log final size            (logistic, richards)
    PARAMETER(log_p);      // log Richards shape        (richards)
    PARAMETER(log_nbdisp); // log nb dispersion         (nbinom)
    PARAMETER(log_b);      // log baseline growth rate  (baseline_flag)

    // probably makes sense to parameterize this as
    //  log_r + sd_log_r*eps_log_r(i)
    // where eps_log_r(i) ~ N(0,1)
    // *not* log_r + eps_log_r(i);  eps_log_r(i) ~ N(0, sd_log_r)
    // (if we end up going to tmbstan etc.)

    // as general as possible, for every grouping variable (continent/country)
    // we will have a *matrix* which is n_grp X n_pars
    // each row will be MVN(0,C)  where C is the correlation matrix
    // then we can still say log_r + sd_log_r*epsmat(country(i),0)
    //  {where the column index is an integer corresponding to a parameter}
    //  might get a little tricky dealing with models with different numbers
    // of parameters i.e. {r=0, K=1, nbdisp=2} or {r=0, K=1, p=2, nbdisp=3} ?

    // if we have multiple grouping variables then we might need a list??
    // does TMB even do lists??

    // in order to know the mapping between parameters and columns of the
    // matrix I might pass a vector of integers, or even a bunch of separate
    // integer values r_col, K_col, thalf_col, nbdisp_col ...
    // log_r + sd_log_r* epsmat(country(i),r_col)
    
    // when you get to the end and you want to add the Likelihood(parameters|MVN structure) you would add dmvnorm([row i], mean 0, Correlation matrix)
    
    // Inverse link-transformed parameters
    Type r = exp(log_r);
    Type thalf = exp(log_thalf);
    Type p = exp(log_p);

    // Other stuff
    Type log_nbxsvar; // log nb excess variance (nbinom)
    
    // Incidence curves
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
	    log_cum_inc(i) = log_c0 + r * t(i);
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
        log_int_inc(i) = logspace_sub(log_cum_inc(i+1), log_cum_inc(i));
	
        if (!isNA(x(i)))
        {
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


    // PREDICT AT NEW TIMES ====================================================

    // FIXME: DRY!
    if (predict_flag == 1)
    {
        // Data
        DATA_VECTOR(t_new); // time,  length M+1
	DATA_VECTOR(x_new); // cases, length M

	vector<Type> log_cum_inc_new(t_new.size()); // log cumulative incidence, length M+1
	vector<Type> log_int_inc_new(x_new.size()); // log interval incidence,   length M

	// Evaluate log cumulative incidence at new times
	for (int i = 0; i < t_new.size(); i++)
	{
	    switch (curve_flag)
	    {
	    case exponential:
	        log_cum_inc_new(i) = log_c0 + r * t_new(i);
		break;
	    case logistic:
	        log_cum_inc_new(i) = log_K - logspace_add(Type(0), -r * (t_new(i) - thalf));
		break;
	    case richards:
	        log_cum_inc_new(i) = log_K - logspace_add(Type(0), logspace_sub(p * log(2), Type(0)) - r * p * (t_new(i) - thalf)) / p;
		break;
	    }

	    if (baseline_flag == 1)
	    {
	        log_cum_inc_new(i) = logspace_add(log_b + log(t_new(i)), log_cum_inc_new(i));
	    }
	}

	// Compute log interval incidence
	for (int i = 0; i < x_new.size(); i++)
	{
	    log_int_inc_new(i) = logspace_sub(log_cum_inc_new(i+1), log_cum_inc_new(i));
        }

        ADREPORT(log_cum_inc_new)
	ADREPORT(log_int_inc_new)
    }

    
    return nll;
}
