#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "utils.h"
#include "enum.h"
#include "curve.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Flags
    DATA_INTEGER(curve_flag);    // prediction model     (enum curve)
    DATA_INTEGER(excess_flag);   // excess mortality     (1=yes,  0=no)
    DATA_INTEGER(distr_flag);    // observation model    (enum distr)
    DATA_INTEGER(weekday_flag);  // day-of-week effects  (1=yes,  0=no)
    DATA_INTEGER(sparse_X_flag); // sparse X matrix      (1=yes,  0=no)
    DATA_INTEGER(predict_flag);  // predict              (1=yes,  0=no)
    bool excess   = (excess_flag   == 1);
    bool weekday  = (weekday_flag  == 1);
    bool sparse_X = (sparse_X_flag == 1);
    bool predict  = (predict_flag  == 1);
    
    // Data
    // time series
    DATA_VECTOR(t); // length=n
    DATA_VECTOR(x); // length=n-N
    // model matrices
    DATA_MATRIX(Xd);        // nrow=N,  ncol=sum(fncoef) or 0 if sparse_X=true
    DATA_SPARSE_MATRIX(Xs); // nrow=N,  ncol=sum(fncoef) or 0 if sparse_X=false
    DATA_SPARSE_MATRIX(Z);  // nrow=N,  ncol=sum(rncoef)
    DATA_MATRIX(Yo);        // nrow=N,  ncol=p

    // Metadata
    // segment lengths
    DATA_IVECTOR(slen); // length=N
    // earliest day-of-week
    DATA_IVECTOR(dow); // length=N or 0 if weekday=false,  val={0,...,6}  (0=reference)
    // number of coefficients
    DATA_IVECTOR(fncoef); // length=p
    DATA_IVECTOR(rncoef); // length=p
    // coefficients factored by response vector
    DATA_IVECTOR(fpar); // length=sum(fncoef),  val={0,...,p-1}
    DATA_IVECTOR(rpar); // length=sum(rncoef),  val={0,...,p-1}
    // dimensions of random effects blocks
    DATA_IVECTOR(rnrow); // length=M,  val={1,...,p}
    DATA_IVECTOR(rncol); // length=M
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
    int M = rnrow.size(); // number of random effects blocks
    int N = slen.size(); // number of segments (fitting windows)
    int p = Yo.cols(); // number of response vectors (nonlinear model parameters)
    bool anyRE = (rncoef.sum() > 0);
    
    // Parameters
    // fixed effects coefficients
    PARAMETER_VECTOR(beta); // length=sum(fncoef)
    // random effects coefficients (unit variance)
    PARAMETER_VECTOR(b); // length=sum(rncoef)
    // log sd random effects coefficients
    PARAMETER_VECTOR(log_sd_b); // length=sum(rnrow)
    
    
    // Compute nonlinear model parameter values ================================
    // (link scale) in each segment

    // 1. Initialize (N x p) matrix of parameter values with offsets
    matrix<Type> Y = Yo;
    
    // 2. Add fixed effects component:

    // 2.(a) Construct block matrix `beta_as_matrix` from vector `beta`
    Eigen::SparseMatrix<Type> beta_as_matrix(beta.size(), p);
    beta_as_matrix.reserve(fncoef); // declare nnz
    for (int i = 0; i < beta.size(); i++) // loop over `beta` elements
    {
        beta_as_matrix.insert(i, fpar(i)) = beta(i);
    }

    // 2.(b) Matrix multiply  X * beta_as_matrix
    //                          = (N x sum(fncoef)) * (sum(fncoef) x p)
    //                          = (N x p)
    if (sparse_X)
    {
        Y += Xs * beta_as_matrix;
    }
    else
    {
        Y += Xd * beta_as_matrix;
    }
    
    // 3. Preserve (offsets + fixed effects) component for SIMULATE
    matrix<Type> Y_simulate = Y;

    // 4. Add random effects component:

    // 4.(a) Set up random effects infrastructure

    // A list of random effects blocks gathering related elements of `b`:
    // columns of each block are samples from a zero-mean, unit-variance,
    // multivariate normal distribution with a common covariance matrix
    vector< matrix<Type> > block_list(M);

    // A list of s.d. vectors corresponding elementwise to `block_list`
    vector< vector<Type> > sd_list(M);

    // A list of correlation matrices corresponding elementwise to `block_list`
    vector< matrix<Type> > cor_list(M);

    // `b` with elements scaled by the appropriate s.d.
    vector<Type> b_scaled(b.size());

    // `b_scaled` with elements arranged in `p` columns of a matrix
    // by relation to a nonlinear model parameter
    matrix<Type> b_scaled_as_matrix(b.size(), p);

    if (anyRE)
    {
        // Correlation matrix in univariate case
        matrix<Type> cor1d(1, 1);
	cor1d(0, 0) = Type(1);

	// Initialize list elements
        // NB: Below assumes that `b` and `log_sd_b` are ordered
	//     so that no permutation is needed when filling the
	//     random effects blocks (in column-major order) and
	//     s.d. vectors.
        for (int m = 0, i1 = 0, i2 = 0; m < M; m++) // loop over list elements
	{
	    // Form random effects block 
	    matrix<Type> block(rnrow(m), rncol(m));
	    for (int j = 0; j < rncol(m); j++) // loop over block columns
	    {
		block.col(j) = b.segment(i1, rnrow(m));
		i1 += rnrow(m); // increment `b` index
	    }
	    block_list(m) = block;

	    // Form s.d. vector
	    sd_list(m) = exp(log_sd_b.segment(i2, rnrow(m)));
	    i2 += rnrow(m); // increment `log_sd_b` index

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

	// 4.(b) Scale random effects coefficients (unit variance)
	//       by corresponding standard deviations
	for (int m = 0, i = 0; m < M; m++) // loop over blocks
	{
	    vector<Type> v(rnrow(m));
	    for (int j = 0; j < rncol(m); j++) // loop over block columns
	    {
	        v = block_list(m).col(j);
	        b_scaled.segment(i, rnrow(m)) = sd_list(m) * v;
		i += rnrow(m); // increment reference index
	    }
	}

	// 4.(c) Construct block matrix `b_scaled_as_matrix` from vector `b_scaled`
	b_scaled_as_matrix.fill(Type(0));
        for (int i = 0; i < b.size(); i++)
	{
	    b_scaled_as_matrix(i, rpar(i)) = b_scaled(i);
	}

	// 4.(d) Matrix multiply  Z * b_scaled_as_matrix
	//                          = (N x sum(rncoef)) * (sum(rncoef) x p)
	//                          = (N x p) 
	Y += Z * b_scaled_as_matrix;
    }
    vector<Type> Y_as_vector = Y.vec();
    ADREPORT(Y_as_vector);


    // Compute likelihood ======================================================

    // Log curve
    vector<Type> log_curve = eval_log_curve(t, slen, curve_flag, excess, Y,
					    j_log_r, j_log_alpha, j_log_c0,
					    j_log_tinfl, j_log_K,
					    j_logit_p, j_log_a, j_log_b);

    // Log cases
    vector<Type> log_cases = eval_log_cases(log_curve, slen, weekday, dow, Y,
					    j_log_w1, j_log_w2, j_log_w3,
					    j_log_w4, j_log_w5, j_log_w6);

    // Negative log likelihood
    Type nll = Type(0);
    Type log_var_minus_mu;

    for (int s = 0, i = 0; s < N; s++) // loop over segments
    {
        for (int k = 0; k < slen(s) - 1; k++) // loop over within-segment index
	{
	    if (!is_NA_real_(x(i+k)))
	    {
	        switch (distr_flag)
		{
		case pois:
		    nll -= dpois_robust(x(i+k), log_cases(i+k), true);
		    // usage: dpois_robust(x, log_lambda, give_log)
		    break;
		case nbinom:
		    log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
		    nll -= dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
		    // usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)
		    break;
		}
	    }
	}
	i += slen(s) - 1; // increment reference index
    }
    if (anyRE)
    {
    	for (int m = 0; m < M; m++) // loop over blocks
    	{
	    for (int j = 0; j < rncol(m); j++) // loop over block columns
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
	    vector<Type> b_simulate_scaled(b.size());
	    for (int m = 0, i = 0; m < M; m++) // loop over blocks
	    {
	        for (int j = 0; j < rncol(m); j++) // loop over block columns
		{
		    b_simulate_scaled.segment(i, rnrow(m)) = sd_list(m) * density::MVNORM(cor_list(m)).simulate();
		    i += rnrow(m); // increment reference index
		}
	    }
	    matrix<Type> b_simulate_scaled_as_matrix(b.size(), p);
	    b_simulate_scaled_as_matrix.fill(Type(0));
	    for (int i = 0; i < b.size(); i++)
	    {
		b_simulate_scaled_as_matrix(i, rpar(i)) = b_simulate_scaled(i);
	    }
	    Y_simulate += Z * b_simulate_scaled_as_matrix;
	}
        // vector<Type> Y_simulate_as_vector = Y_simulate.vec();
	// ADREPORT(Y_simulate_as_vector);

        log_curve = eval_log_curve(t, slen, curve_flag, excess, Y_simulate,
				   j_log_r, j_log_alpha, j_log_c0,
				   j_log_tinfl, j_log_K,
				   j_logit_p, j_log_a, j_log_b);
	
        log_cases = eval_log_cases(log_curve, slen, weekday, dow, Y_simulate,
				   j_log_w1, j_log_w2, j_log_w3,
				   j_log_w4, j_log_w5, j_log_w6);

	for (int s = 0, i = 0; s < N; s++) // loop over segments
	{
	    for (int k = 0; k < slen(s) - 1; k++) // loop over within-segment index
	    {
	        switch(distr_flag)
		{
		case pois:
		    x(i+k) = rpois(exp(log_cases(i+k)));
		    // usage: rpois(lambda)
		    break;
		case nbinom:
		    x(i+k) = rnbinom_robust(log_cases(i+k), Y_simulate(s, j_log_nbdisp));
		    // usage: rnbinom_robust(log_mu, log_size)
		    break;
		}
	    }
	    i += slen(s) - 1; // increment reference index
	}
	REPORT(x);
    }


    // Predict incidence =======================================================
    
    if (predict)
    {
        // Flags
        DATA_IVECTOR(what_flag); // predict  (1=yes,       0=no)
	DATA_INTEGER(se_flag);   // report   (1=ADREPORT,  0=REPORT)
	bool predict_lii = (what_flag(0) == 1);
	bool predict_lci = (what_flag(1) == 1);
	bool predict_lrt = (what_flag(2) == 1);
	bool report_se   = (se_flag      == 1);

	// Data
	// time points
        DATA_VECTOR(t_predict); // length=n'
	// model matrices
	// (subsets of rows of original model matrices)
	DATA_MATRIX(Xd_predict);        // nrow=N',  ncol=sum(fncoef) or 0 if sparse_X=true
	DATA_SPARSE_MATRIX(Xs_predict); // nrow=N',  ncol=sum(fncoef) or 0 if sparse_X=false
	DATA_SPARSE_MATRIX(Z_predict);  // nrow=N',  ncol=sum(rncoef)
	DATA_MATRIX(Yo_predict);        // nrow=N',  ncol=p

	// Metadata
	// segment lengths
	DATA_IVECTOR(slen_predict); // length=N'
	// earliest day-of-week
        DATA_IVECTOR(dow_predict); // length=N' or 0 if weekday=false,  val={0,...,6}  (0=reference)
	// misc.
	int N_predict = slen_predict.size();

	// (N' x p) matrix of nonlinear model parameter values (link scale)
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
	    Y_predict += Z_predict * b_scaled_as_matrix;
    	}

	// Log curve
        vector<Type> log_curve_predict =
	    eval_log_curve(t_predict, slen_predict, curve_flag, excess, Y_predict, 
			   j_log_r, j_log_alpha, j_log_c0,
			   j_log_tinfl, j_log_K,
			   j_logit_p, j_log_a, j_log_b);

        // Log cases
	vector<Type> log_cases_predict =
	    eval_log_cases(log_curve_predict, slen_predict, weekday, dow_predict, Y_predict,
			   j_log_w1, j_log_w2, j_log_w3,
			   j_log_w4, j_log_w5, j_log_w6);


	if (predict_lii)
	{
	    // Log interval incidence
	    vector<Type> log_int_inc = log_cases_predict;
	    if (report_se)
	    {
	        ADREPORT(log_int_inc);
	    }
	    else
	    {
	        REPORT(log_int_inc);
	    }
	}
	
	if (predict_lci)
	{
	    // Log cumulative incidence
	    // (since the earliest time point in the current segment)
	    vector<Type> log_cum_inc(log_cases.size());
	    for (int s = 0, i = 0; s < N_predict; s++) // loop over segments
	    {
	        vector<Type> log_cases_predict_segment = log_cases_predict.segment(i, slen_predict(s) - 1);
		log_cum_inc.segment(i, slen_predict(s) - 1) = logspace_cumsum_1(log_cases_predict_segment);
		i += slen_predict(s) - 1; // increment reference index
	    }
	    
	    if (report_se)
	    {
	        ADREPORT(log_cum_inc);
	    }
	    else
	    {
	        REPORT(log_cum_inc);
	    }
	}
	
	if (predict_lrt)
	{
	    // Log per capita growth rate
	    // (approximated by local linear regression if weekend=true)
	    vector<Type> log_rt =
	        eval_log_rt(t_predict, log_curve_predict, log_cases_predict, slen_predict,
			    curve_flag, excess, weekday, Y_predict,
			    j_log_r, j_log_alpha, j_log_K, j_logit_p, j_log_a, j_log_b);
	    
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
