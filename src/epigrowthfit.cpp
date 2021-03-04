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
    DATA_INTEGER(excess_flag);   // excess    (1=yes,  0=no)
    DATA_INTEGER(distr_flag);    // distr     (enum distr)
    DATA_INTEGER(weekday_flag);  // weekday   (1=yes,  0=no)
    DATA_INTEGER(sparse_X_flag); // sparse X  (1=yes,  0=no)
    DATA_INTEGER(predict_flag);  // predict   (1=yes,  0=no)
    bool excess   = (excess_flag   == 1);
    bool weekday  = (weekday_flag  == 1);
    bool sparse_X = (sparse_X_flag == 1);
    bool predict  = (predict_flag  == 1);
    
    // Data
    // time series
    DATA_VECTOR(t); // length=N
    DATA_VECTOR(x); // length=N-w
    // segment lengths
    DATA_IVECTOR(slen); // length=w
    // earliest day-of-week
    DATA_IVECTOR(dow0); // length=w or 0 if weekday=false,  val={0,...,6}  (0=reference)
    // number of coefficients for each nonlinear model parameter
    DATA_IVECTOR(fncoef); // length=p
    DATA_IVECTOR(rncoef); // length=p
    // coefficients factored by nonlinear model parameter
    DATA_IVECTOR(fpar); // length=sum(fncoef),  val={0,...,p-1}
    DATA_IVECTOR(rpar); // length=sum(rncoef),  val={0,...,p-1}
    // dimensions each random effects block
    DATA_IVECTOR(rnrow); // length=M,  val={1,...,p}
    DATA_IVECTOR(rncol); // length=M
    // model matrices
    DATA_MATRIX(Xd);        // nrow=w,  ncol=sum(fncoef) or 0 if sparse_X=true
    DATA_SPARSE_MATRIX(Xs); // nrow=w,  ncol=sum(fncoef) or 0 if sparse_X=false
    DATA_SPARSE_MATRIX(Z);  // nrow=w,  ncol=sum(rncoef)
    DATA_MATRIX(Yo);        // nrow=w,  ncol=p
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
    int N = t.size();
    int M = rnrow.size();
    int p = fncoef.size();
    int w = slen.size();
    bool anyRE = (rncoef.sum() > 0);
    
    // Parameters
    // fixed effects coefficients
    PARAMETER_VECTOR(beta); // length=sum(fncoef)
    // random effects coefficients, unit variance
    PARAMETER_VECTOR(b); // length=sum(rncoef)
    // log sd random effects coefficients
    PARAMETER_VECTOR(log_sd_b); // length=sum(rnrow)
    
    
    // Compute nonlinear model parameter values ================================
    // (link scale) in each window

    // (w x p) matrix of parameter value initialized to matrix of offsets
    matrix<Type> Y = Yo;
    
    // Fixed effects component:
    // 1. construct block matrix `beta_as_matrix` from vector `beta`
    // 2. matrix multiply
    //    X * beta_as_matrix
    //    = (w x sum(fncoef)) * (sum(fncoef) x p)
    //    = (w x p)

    Eigen::SparseMatrix<Type> beta_as_matrix(beta.size(), p);
    beta_as_matrix.reserve(fncoef); // declare nnz
    
    for (int i = 0; i < beta.size(); i++) // loop over `beta` elements
    {
        beta_as_matrix.insert(i, fpar(i)) = beta(i);
    }
    
    if (sparse_X)
    {
        Y += Xs * beta_as_matrix;
    }
    else
    {
        Y += Xd * beta_as_matrix;
    }
    matrix<Type> Y_sim = Y; // preserve (offsets + fixed effects) component for SIMULATE
    
    // Random effects component:
    // 1. Construct block matrix `b_scaled_as_matrix` from vectors `b` and `log_sd_b`
    // 2. Matrix multiply
    //    Z * b_scaled_as_matrix
    //    = (w x sum(rncoef)) * (sum(rncoef) x p)
    //    = (w x p) 

    vector<Type> b_scaled(b.size());
    matrix<Type> b_scaled_as_matrix(b.size(), p);

    // A list of random effects blocks gathering related elements of `b`:
    // columns of each block are samples from a zero-mean, unit-variance,
    // multivariate normal distribution with a common covariance matrix
    vector< matrix<Type> > block_list(M);

    // A list of s.d. vectors corresponding elementwise to `block_list`
    vector< vector<Type> > sd_list(M);

    // A list of correlation matrices corresponding elementwise to `block_list`
    vector< matrix<Type> > cor_list(M);
    
    if (anyRE)
    {
        // Correlation matrix in univariate case
        matrix<Type> cor1d(1, 1);
	cor1d(0, 0) = Type(1);

	// Initialize list elements
        // NB: Assuming here that `b` and `log_sd_b` are ordered
	//     so that no permutation is needed when filling the
	//     random effects blocks (in column-major order) and
	//     s.d. vectors
        for (int m = 0, k1 = 0, k2 = 0; m < M; m++) // loop over list elements
	{
	    // Form random effects block 
	    matrix<Type> block(rnrow(m), rncol(m));
	    for (int j = 0; j < rncol(m); j++) // loop over block columns
	    {
		block.col(j) = b.segment(k1, rnrow(m));
		k1 += rnrow(m); // increment `b` index
	    }
	    block_list(m) = block;

	    // Form s.d. vector
	    vector<Type> sd = exp(log_sd_b.segment(k2, rnrow(m)));
	    sd_list(m) = sd;
	    k2 += rnrow(m); // increment `log_sd_b` index

	    // Form correlation matrix
	    // NB: UNSTRUCTURED_CORR() does not tolerate `ltri.size() == 0`
	    if (rnrow(m) == 1)
	    {
		cor_list(m) = cor1d;
	    }
	    else
	    {
	        vector<Type> ltri(rnrow(m) * (rnrow(m) - 1) / 2);
		cor_list(m) = density::UNSTRUCTURED_CORR(ltri).cov();
	    }
	}
	REPORT(block_list);
	REPORT(sd_list);
	REPORT(cor_list);

	// Construct `b_scaled` ------------------------------------------------

	for (int m = 0, k = 0; m < M; m++) // loop over blocks
	{
	    for (int j = 0; j < rncol(m); j++) // loop over block columns
	    {
	        b_scaled.segment(k, rnrow(m)) = sd_list(m) * block_list(m).col(j);
		k += rnrow(m); // increment `b_scaled` index
	    }
	}

	// Construct `b_scaled_as_matrix` --------------------------------------

	b_scaled_as_matrix.fill(Type(0));
        for (int i = 0; i < b_scaled.size(); i++) // loop over `b_scaled` elements
	{
	    b_scaled_as_matrix(i, rpar(i)) = b_scaled(i);
	}
	Y += Z * b_scaled_as_matrix;
    }
    vector<Type> Y_as_vector = Y.vec();
    ADREPORT(Y_as_vector);

    
    // Compute likelihood ======================================================

    // Log curve
    vector<Type> log_curve =
      eval_log_curve(t, Y, slen, curve_flag, excess,
		     j_log_r, j_log_alpha, j_log_c0,
		     j_log_tinfl, j_log_K,
		     j_logit_p, j_log_a, j_log_b);
    // Log cases
    vector<Type> log_cases =
      compute_log_cases(log_curve, slen, Y, dow0, weekday,
		        j_log_w1, j_log_w2, j_log_w3,
		        j_log_w4, j_log_w5, j_log_w6);
    
    // Negative log likelihood
    Type nll = Type(0);
    Type log_var_minus_mu;
    
    for (int s = 0, k = 0; s < w; s++) // loop over segments
    {
        for (int i = 0; i < slen(s); i++) // loop over within-segment index
	{
	    if (!isNA_real_(x(k+i)))
	    {
	        switch (distr_flag)
		{
		case pois:
		    nll -= dpois_robust(x(k+i), log_cases(k+i), true);
		    // usage: dpois_robust(x, log_lambda, give_log)
		    break;
		case nbinom:
		    log_var_minus_mu = Type(2) * log_cases(k+i) - Y(s, j_log_nbdisp);
		    nll -= dnbinom_robust(x(k+i), log_cases(k+i), log_var_minus_mu, true);
		    // usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)
		    break;
		}
	    }
	}
	k += slen(s); // increment `x` index
    }
    if (anyRE)
    {
    	for (int m = 0; m < M; m++) // loop over blocks
    	{
	    for (int j = 0; j < rncols(m); l++) // loop over block columns
    	    {
    	        nll += density::MVNORM(cor_list(m))(block_list(m).col(j));
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
    
    if (predict)
    {
        // Flags
        DATA_IVECTOR(what_flag); // length=3,  val={0,1}
	DATA_INTEGER(se_flag); // report  (1=ADREPORT,  0=REPORT)
	bool report_lii = (what_flag(1) == 1);
	bool report_lci = (what_flag(0) == 1);
	bool report_lrt = (what_flag(2) == 1);
	bool report_se  = (se_flag == 1);

	// Data
	// time points
        DATA_VECTOR(t_predict); // length=N'
	// segment lengths
	DATA_IVECTOR(slen_predict); // length=w'
	// earliest day-of-week
        DATA_IVECTOR(dow0_predict); // length=w' or 0 if weekday=false,  val={0,...,6}  (0=reference)
	// model matrices
	// (subsets of rows of original model matrices)
	DATA_MATRIX(Xd_predict);        // nrow=w',  ncol=sum(fncoef) or 0 if sparse_X=true
	DATA_SPARSE_MATRIX(Xs_predict); // nrow=w',  ncol=sum(fncoef) or 0 if sparse_X=false
	DATA_SPARSE_MATRIX(Z_predict);  // nrow=w',  ncol=sum(rncoef)
	DATA_MATRIX(Yo_predict);        // nrow=w',  ncol=p
    
	matrix<Type> Y_predict = Yo_predict;
	if (sparse_X)
	{
	    Y_predict += Xs_predict * beta_as_matrix;
	}
	else
	{
	    Y_predict += Xd_predict * beta_as_matrix;
	}
	if (anyRE)
    	{
	    Y_new += Z_new * b_scaled_as_matrix;
    	}

	// Log curve
        vector<Type> log_curve_predict =
	  eval_log_curve(t_predict, Y_predict, slen_predict, curve_flag, excess,
			 j_log_r, j_log_alpha, j_log_c0,
			 j_log_tinfl, j_log_K, j_logit_p,
			 j_log_a, j_log_b);

        // Log interval incidence
	vector<Type> log_int_inc =
	  compute_log_cases(log_curve_predict, slen_predict, Y_predict, dow0_predict, weekday,
			    j_log_w1, j_log_w2, j_log_w3,
			    j_log_w4, j_log_w5, j_log_w6);

        // Log cumulative incidence
        vector<Type> log_cum_inc(log_int_inc.size());

	if (weekday)
	{
	    for (int s = 0, k = 0; s < slen_predict.size(); s++) // loop over segments
	    {
	        log_cum_inc(k) = log_int_inc(k);
		for (int i = 1; i < slen_predict(s); s++) // loop over within-segment index
		{
		    log_cum_inc(k+i) = logspace_add(log_cum_inc(k+i-1), log_int_inc(k+i));
		}
		k += slen_predict(s); // increment `log_cum_inc` index
	    }
	}
	else
	{
	    for (int s = 0, k1 = 0, k2 = 0; s < slen_predict.size(); s++) // loop over segments
	    {
	        for (int i = 0; i < slen_predict(s); s++) // loop over within-segment index
		{
		    log_cum_inc(k1+i) = logspace_sub(log_curve(k2+i+1), log_curve(k2));
		}
		k1 += slen_predict(s); // increment `log_cum_inc` index
		k2 += slen_predict(s) + 1; // increment `log_curve` index
	    }
	}

	if (report_lii)
	{
	    if (report_se)
	    {
	        ADREPORT(log_int_inc);
	    }
	    else
	    {
	        REPORT(log_int_inc);
	    }
	}
	if (report_lci)
	{
	    if (report_se)
	    {
	        ADREPORT(log_cum_inc);
	    }
	    else
	    {
	        REPORT(log_cum_inc);
	    }
	}
	
	if (report_lrt)
	{
	    vector<Type> log_rt_new = eval_log_rt(log_cum_inc_new, Y_new, curve_flag,
						  j_log_r, j_log_alpha, j_log_K, j_logit_p, j_log_a);
	    if (report_se)
	    {
	        ADREPORT(log_rt);
	    }
	    else
	    {
	        REPORT(log_rt);
	    }
	}
    }

    return nll;
}
