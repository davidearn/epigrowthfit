#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "enum.h"
#include "curve.h"
#include "utils.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Flags
    DATA_INTEGER(curve_flag);    // curve     (enum curve)
    DATA_INTEGER(distr_flag);    // distr     (enum distr)
    DATA_INTEGER(excess_flag);   // excess    (1=yes,      0=no)
    DATA_INTEGER(weekday_flag);  // weekday   (1=yes,      0=no)
    DATA_INTEGER(sparse_X_flag); // X format  (1=sparse,   0=dense)
    DATA_INTEGER(predict_flag);  // predict   (1=yes,      0=no)
    bool excess   = (excess_flag   == 1);
    bool weekday  = (weekday_flag  == 1);
    bool sparse_X = (sparse_X_flag == 1);
    bool predict  = (predict_flag  == 1);
    
    // Data
    // time series
    DATA_VECTOR(t); // length=N
    DATA_VECTOR(x); // length=N
    // window lengths
    DATA_IVECTOR(wl); // length=w
    // day-of-week of earliest date: 0 -> Sunday, etc.
    DATA_IVECTOR(dow0); // length=w
    // number of coefficients for each nonlinear model parameter
    DATA_IVECTOR(fnc); // length=p
    DATA_IVECTOR(rnc); // length=p
    // design matrices
    DATA_MATRIX(Xd);
    DATA_SPARSE_MATRIX(Xs); // nrow=N,  ncol=sum(fnc)
    DATA_SPARSE_MATRIX(Z);  // nrow=N,  ncol=sum(rnc)
    // parameter indices
    DATA_INTEGER(j_log_r);
    DATA_INTEGER(j_log_alpha);
    DATA_INTEGER(j_log_c0);
    DATA_INTEGER(j_log_tinfl);
    DATA_INTEGER(j_log_K);
    DATA_INTEGER(j_logit_p);
    DATA_INTEGER(j_log_a);
    DATA_INTEGER(j_log_b);
    DATA_INTEGER(j_log_nbdisp);
    DATA_INTEGER(j_log_w1);
    DATA_INTEGER(j_log_w2);
    DATA_INTEGER(j_log_w3);
    DATA_INTEGER(j_log_w4);
    DATA_INTEGER(j_log_w5);
    DATA_INTEGER(j_log_w6);
    // misc.
    int w = wl.size();
    int p = fnc.size();
    bool anyRE = (rnc.sum() > 0);
    
    // Parameters
    // fixed effects
    PARAMETER_VECTOR(beta); // length=sum(fnc)
    // log sd random effects
    PARAMETER_VECTOR(log_sd_b); // length=with(Zinfo, nlevels(interaction(par, term, group, drop = TRUE)))
    // random effects
    PARAMETER_VECTOR(b); // length=sum(rnc) 

    
    // Compute parameter values at every time point ============================
    // NOTE: parameters are constant within but not across fitting windows

    // (N x np) matrix of parameter values
    matrix<Type> Y(t.size(), p);
    
    // Fixed effects component:
    // matrix multiply
    // X * beta2
    // = (N x sum(fnc)) * (sum(fnc) x p)
    // = (N x p)
    // need to construct block matrix `beta2` from vector `beta`

    Eigen::SparseMatrix<Type> beta2(fnc.sum(), p);
    beta2.reserve(fnc);
    
    for (int j = 0, i = 0; j < p; j++) // loop over parameters
    {
	for (int l = 0; l < fnc(j); l++) // loop over coefficients
	{
	    beta2.insert(i + l, j) = beta(i + l);
	}
	i += fnc(j);
    }
    
    if (sparse_X)
    {
        Y = Xs * beta2;
    }
    else
    {
        Y = Xd * beta2;
    }
    matrix<Type> Y_sim = Y; // preserve fixed effects component for SIMULATE
    
    // Random effects component:
    // matrix multiply
    // Z * b2
    // = (N x sum(rnc)) * (sum(rnc) x p)
    // = (N x p)
    // need to construct block matrix `b2` from vector `b`
    // will require some setting up...
    
    matrix<Type> b2(rnc.sum(), p);
    
    // A list of matrices gathering related elements of vector `b`:
    // rows of each matrix follow a zero-mean, unit-variance multivariate
    // normal distribution with a common covariance matrix
    vector< matrix<Type> > re_list(rid.rows());

    // A list of s.d. vectors corresponding elementwise to `re_list`
    vector< vector<Type> > sd_list(rid.rows());

    // A list of correlation matrices corresponding elementwise to `re_list`
    vector< matrix<Type> > cor_list(rid.rows());
	
    if (anyRE)
    {

        // Form random effect matrices -----------------------------------------

	for (int r = 0, k = 0; r < re_list.size(); r++) // loop over RE
	{
	    matrix<Type> re_list_el(rnl(r), rnp(r));
	    for (int j = 0; j < rnp(r); j++) // loop over parameters
	    {
		re_list_el.col(j) = b.segment(k, rnl(r));
		k += rnl(r);
	    }
	    re_list(r) = re_list_el;
	}
	REPORT(re_list);


	// Form scale vectors --------------------------------------------------

	for (int r = 0, k = 0; r < sd_list.size(); r++)
	{
	    vector<Type> sd_list_el = exp(log_sd_b.segment(k, rnp(r)));
	    sd_list(r) = sd_list_el;
	    k += rnp(r);
	}
	REPORT(sd_list);

	
	// Form correlation matrices -------------------------------------------

	// Correlation matrix in univariate case
	matrix<Type> cor_list_el(1, 1);
	cor_list_el(0, 0) = Type(1);

	for (int r = 0; r < cor_list.size(); r++) // loop over RE
	{
	    // UNSTRUCTURED_CORR() does not tolerate `ltri.size() == 0`
	    if (rnp(r) == 1)
	    {
		cor_list(r) = cor_list_el;
		continue;
	    }
	    else
	    {
	        vector<Type> ltri(rnp(r) * (rnp(r) - 1) / 2);
		cor_list(r) = density::UNSTRUCTURED_CORR(ltri).cov();
	    }
	}
	REPORT(cor_list);


	// Construct `b2` ------------------------------------------------------

        b2.fill(Type(0));
	for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	{
	    for (int j = 0, k = 0; j < np; j++) // loop over parameters
	    {
		if (rid(r, j) == 1)
		{
		    b2.block(i, j, rnl(r), 1) = sd_list(r)(k) * re_list(r).col(k);
		    k++;
		}
	    }
	    i += rnl(r);
	}
	Y += Z * b2;
    }


    // Report parameter values in each fitting window ==========================
    // NOTE: parameters are constant within but not across fitting windows

    matrix<Type> Y_short = prune_dupl_rows(Y, wl);
    vector<Type> Y_short_as_vector = Y_short.vec();
    ADREPORT(Y_short_as_vector);
    

    // Compute likelihood ======================================================

    // Log cumulative incidence
    vector<Type> log_cum_inc = eval_log_cum_inc(t, Y, curve_flag, excess,
						j_log_r, j_log_alpha, j_log_c0, j_log_tinfl, j_log_K, j_logit_p, j_log_a, j_log_b);
    // Log interval incidence
    vector<Type> log_int_inc = eval_log_int_inc(log_cum_inc, wl);

    // Negative log likelihood
    Type nll = Type(0);
    Type log_var_minus_mu;
    
    for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++) // loop over time series
    {
        for (int j = 0; j < wl(k) - 1; j++) // loop over time points
	{
	    if (!isNA_real_(x(i2+j)))
	    {
	        switch (distr_flag)
		{
		case pois:
		    nll -= dpois_robust(x(i2+j), log_int_inc(i1+j), true);
		    // usage: dpois_robust(x, log_lambda, give_log)
		    break;
		case nbinom:
		    log_var_minus_mu = Type(2) * log_int_inc(i1+j) - Y(i2+j, j_log_nbdisp);
		    nll -= dnbinom_robust(x(i2+j), log_int_inc(i1+j), log_var_minus_mu, true);
		    // usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)
		    break;
		}
	    }
	}
	i1 += wl(k) - 1;
	i2 += wl(k);
    }
    if (anyRE)
    {
    	for (int r = 0; r < re_list.size(); r++) // loop over RE
    	{
    	    for (int l = 0; l < rnl(r); l++) // loop over levels
    	    {
    	        nll += density::MVNORM(cor_list(r))(re_list(r).row(l));
    	    }
    	}
    }


    // Simulate incidence ======================================================

    SIMULATE
    {
        if (anyRE)
	{
	    matrix<Type> b2_sim(rnl.sum(), np);
	    b2_sim.fill(Type(0));
	    
	    for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	    {
		vector<Type> z(sd_list(r).size());
		for (int l = 0; l < rnl(r); l++) // loop over levels
		{
		    z = sd_list(r) * density::MVNORM(cor_list(r)).simulate();
		    for (int j = 0, k = 0; j < np; j++) // loop over parameters
		    {
			if (rid(r, j) == 1)
			{
			    b2_sim(i + l, j) = z(k);
			    k++;
			}
		    }
		}
		i += rnl(r);
	    }
	    Y_sim += Z * b2_sim;
	}
	// matrix<Type> Y_sim_short = prune_dupl_rows(Y_sim, wl);
	// REPORT(Y_sim_short);

        log_cum_inc = eval_log_cum_inc(t, Y_sim, curve_flag, excess,
				       j_log_r, j_log_alpha, j_log_c0, j_log_tinfl, j_log_K, j_logit_p, j_log_a, j_log_b);
	log_int_inc = eval_log_int_inc(log_cum_inc, wl);

	for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++) // loop over time series
	{
	    for (int j = 0; j < wl(k) - 1; j++) // loop over time points
	    {
	        switch(distr_flag)
		{
		case pois:
		    x(i2+j) = rpois(exp(log_int_inc(i1+j)));
		    // usage: rpois(mu)
		    break;
		case nbinom:
		    x(i2+j) = rnbinom_robust(log_int_inc(i1+j), Y_sim(i2+j, j_log_nbdisp));
		    // usage: rnbinom_robust(log_mu, log_disp)
		    break;
		}
	    }
	    i1 += wl(k) - 1;
	    i2 += wl(k);
	}
	REPORT(x);
    }
    

    // Predict incidence =======================================================
    // NOTE: enforcing on R side that there is only one time series here
    
    if (predict)
    {
        DATA_VECTOR(t_new); // length=N_new
	DATA_SPARSE_MATRIX(Xs_new); // nrow=N_new,  ncol=sum(nlevels(FE i))
    	DATA_MATRIX(Xd_new); 
	DATA_SPARSE_MATRIX(Z_new); // nrow=N_new,  ncol=sum(nlevels(RE i))
	DATA_IVECTOR(predict_lci_lii_lrt_flag); // predict (log_cum_inc, log_int_inc, log_rt)  (1=do,  0=don't)
	DATA_INTEGER(se_flag); // report  (1=ADREPORT,  0=REPORT)

	matrix<Type> Y_new(t_new.size(), np);
	if (sparse_X)
	{
	    Y_new = Xs_new * beta2;
	}
	else
	{
	    Y_new = Xd_new * beta2;
	}
        if (anyRE)
    	{
	    Y_new += Z_new * b2;
    	}

	vector<Type> log_cum_inc_new = eval_log_cum_inc(t_new, Y_new, curve_flag, excess,
							j_log_r, j_log_alpha, j_log_c0, j_log_tinfl, j_log_K, j_logit_p, j_log_a, j_log_b);
	
	if (predict_lci_lii_lrt_flag(0) == 1)
	{
	    vector<int> wl_new(1);
	    wl_new(0) = t_new.size();
	    vector<Type> log_int_inc_new = eval_log_int_inc(log_cum_inc_new, wl_new);
	    if (se_flag == 1)
	    {
	        ADREPORT(log_int_inc_new);
	    }
	    else
	    {
	        REPORT(log_int_inc_new);
	    }
	}
	if (predict_lci_lii_lrt_flag(1) == 1)
	{
	    if (se_flag == 1)
	    {
	        ADREPORT(log_cum_inc_new);
	    }
	    else
	    {
	        REPORT(log_cum_inc_new);
	    }
	}
	if (predict_lci_lii_lrt_flag(2) == 1)
	{
	    vector<Type> log_rt_new = eval_log_rt(log_cum_inc_new, Y_new, curve_flag,
						  j_log_r, j_log_alpha, j_log_K, j_logit_p, j_log_a);
	    if (se_flag == 1)
	    {
	        ADREPORT(log_rt_new);
	    }
	    else
	    {
	        REPORT(log_rt_new);
	    }
	}
    }

    return nll;
}
