#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "utils.h"
#include "enum.h"
#include "curve.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Parameters --------------------------------------------------------------
  
    // Concatenated fixed effects coefficient vectors
    // * length=sum(beta_seg_len)
    PARAMETER_VECTOR(beta);

    // Concatenated random effects coefficient vectors
    // * unit variance scale
    // * length=sum(b_seg_len)
    PARAMETER_VECTOR(b); 

    // Concatenated random effects covariance parameters
    // * length=sum(f(block_rows)), f(n)=n*(n+1)/2
    PARAMETER_VECTOR(theta);

    // Data --------------------------------------------------------------------

    // Concatenated time series segments
    // * length=n
    DATA_VECTOR(t);
    // * length=n-N
    DATA_VECTOR(x);

    // Segment lengths
    // * length=N
    DATA_IVECTOR(t_seg_len);

    // Segment initial days-of-week
    // * length=N if weekday=true, otherwise length=0
    // * val={0,...,6} where 0=reference
    DATA_IVECTOR(dow);

    // Number of segments
    int N = t_seg_len.size();

    // Offset component of response matrix `Y`
    // * dim=(N, p)
    DATA_MATRIX(Yo);

    // Number of nonlinear model parameters
    int p = Yo.cols();
    
    // Column indices of nonlinear model parameters in response matrix `Y`
    // * val={0,...p-1} if parameter is used, otherwise val=-1
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

    // Combined fixed effects model matrix
    // * dense format
    // * dim=(N, sum(beta_seg_len)) if sparse_X=false, otherwise dim=(N, 0)
    DATA_MATRIX(Xd);
    // * sparse format
    // * dim=(N, sum(beta_seg_len)) if sparse_X=true,  otherwise dim=(N, 0)
    DATA_SPARSE_MATRIX(Xs);

    // Combined random effects model matrix
    // * dim=(N, sum(b_seg_len))
    DATA_SPARSE_MATRIX(Z);

    // Indicator for random effects model
    bool any_RE = (Z.cols() > 0);

    // "Factors" splitting model matrix columns by relation
    // to a nonlinear model parameter
    // * length=sum(beta_seg_len)
    // * val={0,...,p-1}
    DATA_IVECTOR(beta_seg_index);
    // * length=sum(b_seg_len)
    // * val={0,...,p-1}
    DATA_IVECTOR(b_seg_index);

    // Number of model matrix columns by nonlinear model parameter
    // * in R: foo_seg_len = as.integer(table(foo_seg_index))
    // * length=p
    DATA_IVECTOR(beta_seg_len); 
    // * length=p
    DATA_IVECTOR(b_seg_len);
    
    // Dimensions of random effects blocks
    // * length=M
    // * val={1,...,p}
    DATA_IVECTOR(block_rows);
    // * length=M
    // * val={nlevels(g)}
    DATA_IVECTOR(block_cols);

    // Number of random effects blocks
    int M = block_rows.size();

    // Flags
    DATA_INTEGER(curve_flag);    // prediction model     (enum.h)
    DATA_INTEGER(distr_flag);    // observation model    (enum.h)
    DATA_INTEGER(excess_flag);   // excess mortality     (1=yes, 0=no)
    DATA_INTEGER(weekday_flag);  // day-of-week effects  (1=yes, 0=no)
    DATA_INTEGER(sparse_X_flag); // sparse format X      (1=yes, 0=no)
    DATA_INTEGER(trace_flag);    // tracing              (0=nothing, 1=degenerate nll terms, 2=1+Y, 3=all nll terms, 4=3+Y)
    DATA_INTEGER(predict_flag);  // predict              (1=yes, 0=no)
    bool excess   = (excess_flag   == 1);
    bool weekday  = (weekday_flag  == 1);
    bool sparse_X = (sparse_X_flag == 1);
    bool trace    = (trace_flag    >  0);
    bool predict  = (predict_flag  == 1);
    

    // Prepare random effects infrastructure ===================================

    // A list of matrices in which to arrange elements of `b`:
    // * blocks, block rows, and block columns should correspond to
    //   random effect terms, nonlinear model parameters, and group levels,
    //   respectively
    // * block columns should represent realizations of a multivariate
    //   normal random variable with zero mean, unit variance, and
    //   unstructured covariance
    vector< matrix<Type> > block_list(M);

    // A list of s.d. vectors corresponding elementwise to `block_list`
    vector< vector<Type> > sd_list(M);

    // A list of correlation matrices corresponding elementwise to `block_list`
    vector< matrix<Type> > cor_list(M);

    // `b` with elements scaled by corresponding standard deviations
    vector<Type> b_scaled(b.size());

    if (any_RE)
    {
        // Correlation matrix of random vector of length 1
        matrix<Type> cor1(1, 1);
	cor1(0, 0) = Type(1);
      
        // Initialize list elements
	int nr;
	int nc;
        for (int m = 0, i1 = 0, i2 = 0; m < M; m++) // loop over list elements
	{
	    nr = block_rows(m);
	    nc = block_cols(m);
	    
	    // Form s.d. vector
	    vector<Type> sd_list_el = exp(theta.segment(i1, nr));
	    sd_list(m) = sd_list_el;
	    i1 += nr; // increment `theta` index

	    // Form correlation matrix
	    // NB: UNSTRUCTURED_CORR() does not tolerate `ltri.size() == 0`
	    if (nr == 1)
	    {
		cor_list(m) = cor1;
	    }
	    else
	    {
	        vector<Type> ltri = theta.segment(i1, nr * (nr - 1) / 2);
		cor_list(m) = density::UNSTRUCTURED_CORR(ltri).cov();
		i1 += ltri.size();
	    }
	    
	    // Form random effects block
	    matrix<Type> block_list_el(nr, nc);
	    vector<Type> v(nr);
	    for (int j = 0; j < nc; j++) // loop over block columns
	    {
	        v = b.segment(i2, nr);
	        block_list_el.col(j) = v;
		b_scaled.segment(i2, nr) = sd_list_el * v;
		i2 += nr; // increment `b` index
	    }
	    block_list(m) = block_list_el;
	}
	REPORT(block_list);
	REPORT(sd_list);
	REPORT(cor_list);
    }

    
    // Compute response matrix =================================================

    // 1. Initialize (N, p) response matrix
    matrix<Type> Y = Yo;
    
    // 2. Add fixed effects component:

    // Arrange elements of `beta` in (sum(beta_seg_len), p) matrix
    Eigen::SparseMatrix<Type> beta_as_matrix(beta.size(), p);
    beta_as_matrix.reserve(beta_seg_len); // declare nnz
    for (int i = 0; i < beta.size(); i++)
    {
        beta_as_matrix.insert(i, beta_seg_index(i)) = beta(i);
    }

    // Increment response matrix with matrix product
    // X * beta_as_matrix
    // = (N, sum(beta_seg_len)) * (sum(beta_seg_len), p)
    // = (N, p)
    Y += (sparse_X ? Xs : Xd) * beta_as_matrix;
    
    // Preserve response matrix without random effects component
    // for simulation code
    matrix<Type> Y_simulate = Y;

    // 3. Add random effects component:

    Eigen::SparseMatrix<Type> b_scaled_as_matrix(b.size(), p);

    if (any_RE)
    {
        // Arrange elements of `b_scaled` in (sum(b_seg_len), p) matrix
        // Eigen::SparseMatrix<Type> b_scaled_as_matrix(b.size(), p); // must be declared outside of conditional scope for prediction code
        b_scaled_as_matrix.reserve(b_seg_len); // declare nnz
	for (int i = 0; i < b.size(); i++)
	{
	    b_scaled_as_matrix.insert(i, b_seg_index(i)) = b_scaled(i);
	}

	// Increment response matrix with matrix product
	// Z * b_scaled_as_matrix
	// = (N, sum(b_seg_len)) * (sum(b_seg_len), p)
	// = (N, p)
	Y += Z * b_scaled_as_matrix;
    }
    vector<Type> Y_as_vector = Y.vec();
    REPORT(Y_as_vector);


    // Compute predictions =====================================================

    // Log curve
    vector<Type> log_curve = eval_log_curve(t, t_seg_len, curve_flag, excess, Y,
					    j_log_r, j_log_alpha, j_log_c0,
					    j_log_tinfl, j_log_K,
					    j_logit_p, j_log_a, j_log_b);

    // Log cases
    vector<Type> log_cases = eval_log_cases(log_curve, t_seg_len, weekday, dow, Y,
					    j_log_w1, j_log_w2, j_log_w3,
					    j_log_w4, j_log_w5, j_log_w6);


    // Compute likelihood ======================================================

    // Negative log likelihood
    Type nll = Type(0);
    Type nll_term;
    Type log_var_minus_mu;
    int width_s;
    int width_k;
    bool print_Y_row;

    if (trace)
    {
        if (trace_flag % 2 == 0)
	{
	    Rcout << "Y = \n" << Y << "\n";
	}
	Rprintf("nll initialized to 0\ncommencing loop over observations\n");
        width_s = nchar(N - 1);
	width_k = nchar(t_seg_len.maxCoeff() - 1);
    }

    for (int s = 0, i = 0; s < N; s++) // loop over segments
    {
        print_Y_row = (trace_flag >= 3);
	
        for (int k = 0; k < t_seg_len(s) - 1; k++) // loop over within-segment index
	{
	    if (!is_NA_real_(x(i+k)))
	    {
	        switch (distr_flag)
		{
		case pois:
		    nll_term = -dpois_robust(x(i+k), log_cases(i+k), true);
		    // usage: dpois_robust(x, log_lambda, give_log)
		    break;
		case nbinom:
		    log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
		    nll_term = -dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
		    // usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)
		    break;
		}

		if (trace)
		{
		    if (trace_flag >= 3)
		    {
		        Rprintf("at index %*d of segment %*d: nll term is %5.6e\n",
			        width_k, k, width_s, s, asDouble(nll_term));
		    }
		    else
		    {
		        if (!is_finite(nll_term))
			{
			    Rprintf("at index %d of segment %d: nll term is non-finite\n", k, s);
			    print_Y_row = true;
			}
			else if (asDouble(nll_term) > 1.0e12)
			{
			    Rprintf("at index %d of segment %d: nll term exceeds 1.0e12\n", k, s);
			    print_Y_row = true;
			}
		    }
		}

		nll += nll_term;
	    }
	}

	if (print_Y_row)
	{
	    Rcout << "Y.row(" << s << ") = " << Y.row(s) << "\n";
	}
	
	i += t_seg_len(s) - 1; // increment reference index
    }

    if (trace)
    {
        Rprintf("loop over observations complete\nnll is %5.6e\n", asDouble(nll));
    }

    if (any_RE) {
        int width_m;
	int width_j;
	
	if (trace)
	{
	    Rprintf("commencing loop over random effects\n");
	    width_m = nchar(M - 1);
	    width_j = nchar(block_cols.maxCoeff() - 1);
	}
      
        for (int m = 0; m < M; m++) // loop over blocks
    	{
	    for (int j = 0; j < block_cols(m); j++) // loop over block columns
    	    {
    	        nll_term = density::MVNORM(cor_list(m))(block_list(m).col(j));

		if (trace)
		{
		    if (trace_flag >= 3)
		    {
		        Rprintf("at column %*d of block %*d: nll term is %5.6e\n",
			       width_j, j, width_m, m, asDouble(nll_term));
		    }
		    else
		    {
			if (!is_finite(nll_term))
			{
			    Rprintf("at column %*d of block %*d: nll term is non-finite\n",
				   width_j, j, width_m, m);
			}
			else if (asDouble(nll_term) > 1.0e12)
			{
			    Rprintf("at column %*d of block %*d: nll term exceeds 10^12\n",
				   width_j, j, width_m, m);
			}
		    }
		}
		
		nll += nll_term;
    	    }
    	}

	if (trace)
	{
	    Rprintf("loop over random effects complete\nnll is %5.6e\n", asDouble(nll));
	}
    }


    // Simulate incidence ======================================================

    SIMULATE
    {
        if (any_RE)
	{
	    vector<Type> b_simulate_scaled(b.size());
	    for (int m = 0, i = 0; m < M; m++) // loop over blocks
	    {
	        for (int j = 0; j < block_cols(m); j++) // loop over block columns
		{
		    b_simulate_scaled.segment(i, block_rows(m)) = sd_list(m) * density::MVNORM(cor_list(m)).simulate();
		    i += block_rows(m); // increment reference index
		}
	    }
	    Eigen::SparseMatrix<Type> b_simulate_scaled_as_matrix(b.size(), p);
	    b_simulate_scaled_as_matrix.reserve(b_seg_len); // declare nnz
	    for (int i = 0; i < b.size(); i++)
	    {
		b_simulate_scaled_as_matrix.insert(i, b_seg_index(i)) = b_simulate_scaled(i);
	    }
	    Y_simulate += Z * b_simulate_scaled_as_matrix;
	}

        log_curve = eval_log_curve(t, t_seg_len, curve_flag, excess, Y_simulate,
				   j_log_r, j_log_alpha, j_log_c0,
				   j_log_tinfl, j_log_K,
				   j_logit_p, j_log_a, j_log_b);
	
        log_cases = eval_log_cases(log_curve, t_seg_len, weekday, dow, Y_simulate,
				   j_log_w1, j_log_w2, j_log_w3,
				   j_log_w4, j_log_w5, j_log_w6);

	for (int s = 0, i = 0; s < N; s++) // loop over segments
	{
	    for (int k = 0; k < t_seg_len(s) - 1; k++) // loop over within-segment index
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
	    i += t_seg_len(s) - 1; // increment reference index
	}
	REPORT(x);
    }


    // Predict incidence =======================================================
    
    if (predict)
    {
        // New data ------------------------------------------------------------

        // Concatenated time series segments
        // * length=n'
        DATA_VECTOR(t_predict);

	// Segment lengths
	// * length=N'
	DATA_IVECTOR(t_predict_seg_len); 
	
	// Segment initial days-of-week
	// * length=N' if weekday=true, otherwise length=0
	// * val={0,...,6} where 0=reference
        DATA_IVECTOR(dow_predict);

	// Number of segments
	int N_predict = t_predict_seg_len.size();

	// Combined fixed effects model matrix
	// * dense format
	// * dim=(N', sum(beta_seg_len)) if sparse_X=false, otherwise dim=(N', 0)
	DATA_MATRIX(Xd_predict);
	// * sparse format
	// * dim=(N', sum(beta_seg_len)) if sparse_X=true,  otherwise dim=(N', 0)
	DATA_SPARSE_MATRIX(Xs_predict);

	// Combined random effects model matrix
	// * dim=(N', sum(b_seg_len))
	DATA_SPARSE_MATRIX(Z_predict);

	// Offset component of response matrix `Y_predict`
	// * dim=(N', p)
	DATA_MATRIX(Yo_predict);

	// Flags
	DATA_IVECTOR(what_flag); // predict  (1=yes, 0=no)
	DATA_INTEGER(se_flag);   // report   (1=ADREPORT, 0=REPORT)
	bool predict_lii = (what_flag(0) == 1);
	bool predict_lci = (what_flag(1) == 1);
	bool predict_lrt = (what_flag(2) == 1);
	bool report_se   = (se_flag      == 1);


	// Compute response matrix ---------------------------------------------
	
	matrix<Type> Y_predict = Yo_predict;
	Y_predict += (sparse_X ? Xs_predict : Xd_predict) * beta_as_matrix;
	if (any_RE)
    	{
	    Y_predict += Z_predict * b_scaled_as_matrix;
    	}

	// Compute predictions -------------------------------------------------

	// Log curve
        vector<Type> log_curve_predict =
	    eval_log_curve(t_predict, t_predict_seg_len, curve_flag, excess, Y_predict, 
			   j_log_r, j_log_alpha, j_log_c0,
			   j_log_tinfl, j_log_K,
			   j_logit_p, j_log_a, j_log_b);

        // Log cases
	vector<Type> log_cases_predict =
	    eval_log_cases(log_curve_predict, t_predict_seg_len, weekday, dow_predict, Y_predict,
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
	    vector<Type> log_cum_inc(log_cases_predict.size());
	    for (int s = 0, i = 0; s < N_predict; s++) // loop over segments
	    {
	        vector<Type> v = log_cases_predict.segment(i, t_predict_seg_len(s) - 1);
		log_cum_inc.segment(i, t_predict_seg_len(s) - 1) = logspace_cumsum_1(v);
		i += t_predict_seg_len(s) - 1; // increment reference index
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
	    // (approximated by local linear regression if weekday=true)
	    vector<Type> log_rt =
	        eval_log_rt(t_predict, log_curve_predict, log_cases_predict, t_predict_seg_len,
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
