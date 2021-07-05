#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "utils.h"
#include "enum.h"
#include "curve.h"

template<class Type>
struct list_of_vectors_t : vector< vector<Type> > {
    list_of_vectors_t(SEXP x)
    {
        (*this).resize(LENGTH(x));
        for (int i = 0; i < LENGTH(x); i++)
	{
            SEXP v = VECTOR_ELT(x, i);
            (*this)(i) = asVector<Type>(v);
        }
    }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
    /* Set up =============================================================== */

    /* Parameters ----------------------------------------------------------- */
  
    /* Concatenated fixed effects coefficient vectors
       - length=sum(beta_seg_len)
    */
    PARAMETER_VECTOR(beta);

    /* Concatenated random effects coefficient vectors
       - unit variance scale
       - length=sum(b_seg_len)
    */
    PARAMETER_VECTOR(b);

    /* Concatenated random effects covariance parameters
       - length=sum(f(block_rows)), f(n)=n*(n+1)/2
    */
    PARAMETER_VECTOR(theta);

    /* Data ----------------------------------------------------------------- */

    /* Concatenated time series segments
       - length=n
    */
    DATA_VECTOR(t);
    /* 
       - length=n-N 
    */
    DATA_VECTOR(x);

    /* Segment lengths
       - length=N
    */
    DATA_IVECTOR(t_seg_len);
    vector<int> x_seg_len = decrement(t_seg_len, 1);

    /* Segment initial days of week
       - length=N if day_of_week=true, otherwise length=0
       - val={0,...,6} where 0=reference
    */
    DATA_IVECTOR(day1);

    /* Number of segments */
    int N = t_seg_len.size();

    /* Offset component of response matrix `Y`
       - dim=(N, p)
    */
    DATA_MATRIX(Yo);

    /* Number of nonlinear model parameters */
    int p = Yo.cols();
    
    /* Column indices of nonlinear model parameters in response matrix `Y`
       - val={0,...p-1} if parameter is used, otherwise val=-1
    */
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

    /* Combined fixed effects model matrix
       - dense format
       - dim=(N, sum(beta_seg_len)) if do_sparse_X=false, otherwise dim=(N, 0)
    */
    DATA_MATRIX(Xd);
    /* 
       - sparse format
       - dim=(N, sum(beta_seg_len)) if do_sparse_X=true,  otherwise dim=(N, 0)
    */
    DATA_SPARSE_MATRIX(Xs);

    /* Combined random effects model matrix
       - sparse format
       - dim=(N, sum(b_seg_len))
    */
    DATA_SPARSE_MATRIX(Z);

    /* Indicator for random effects model */
    bool any_random_effects = (Z.cols() > 0);

    /* "Factors" splitting model matrix columns by relation 
       to a nonlinear model parameter
       - length=sum(beta_seg_len)
       - val={0,...,p-1}
    */
    DATA_IVECTOR(beta_seg_index);
    /*
       - length=sum(b_seg_len)
       - val={0,...,p-1}
    */
    DATA_IVECTOR(b_seg_index);

    /* Number of model matrix columns by nonlinear model parameter
       - in R: *_seg_len = c(table(*_seg_index))
       - length=p
    */
    DATA_IVECTOR(beta_seg_len);
    DATA_IVECTOR(b_seg_len);
    
    /* Dimensions of random effects blocks
       - length=M
       - val={1,...,p}
    */
    DATA_IVECTOR(block_rows);
    /* 
       - length=M
       - val={nlevels(g)}
    */
    DATA_IVECTOR(block_cols);

    /* Number of random effects blocks */
    int M = block_rows.size();

    /* List of vectors of hyperparameters for regularization */
    DATA_STRUCT(regularize_hyperpar, list_of_vectors_t);
    

    /* Flags ---------------------------------------------------------------- */
    
    /* Model of cumulative incidence (enum.h) */
    DATA_INTEGER(curve_flag);
    /* Excess mortality term in cumulative mortality model (1=yes, 0=no) */
    DATA_INTEGER(excess_flag);
    /* Model of observation error (enum.h) */
    DATA_INTEGER(family_flag);
    /* Day of week effects (1=yes, 0=no) */
    DATA_INTEGER(day_of_week_flag);
    /* Priors on nonlinear model parameters (enum.h)
       - length=p
       - val=-1 if no prior
    */
    DATA_IVECTOR(regularize_flag);
    /* Trace
       0=nothing
       1=degenerate nll terms
       2=1+Y
       3=all nll terms
       4=3+Y
    */
    DATA_INTEGER(trace_flag);
    /* X format (1=sparse, 0=dense) */
    DATA_INTEGER(sparse_X_flag);
    /* predict (1=yes, 0=no) */
    DATA_INTEGER(predict_flag);

    bool do_excess      = (excess_flag      == 1);
    bool do_day_of_week = (day_of_week_flag == 1);
    bool do_regularize  = any_geq_zero(regularize_flag);
    bool do_trace       = (trace_flag       >  0);
    bool do_sparse_X    = (sparse_X_flag    == 1);
    bool do_predict     = (predict_flag     == 1);
    
    
    /* Prepare random effects infrastructure ================================ */

    /* List of matrices in which to arrange elements of `b`:
       - blocks, block rows, and block columns should correspond to
         random effect terms, nonlinear model parameters, and group levels,
         respectively
       - block columns should represent realizations of a multivariate
         normal random variable with zero mean, unit variance, and
         unstructured covariance
    */
    vector< matrix<Type> > block_list(M);

    /* List of s.d. vectors corresponding elementwise to `block_list` */
    vector< vector<Type> > sd_list(M);

    /* List of correlation matrices corresponding elementwise to `block_list` */
    vector< matrix<Type> > cor_list(M);

    /* `b` with elements scaled by corresponding standard deviations */
    vector<Type> b_scaled(b.size());

    if (any_random_effects)
    {
        /* Correlation matrix of random vector of length 1 */
        matrix<Type> cor1(1, 1);
	cor1(0, 0) = Type(1);
      
        /* Initialize list elements */
	int nr;
	int nc;
        for (int m = 0, i1 = 0, i2 = 0; m < M; m++) /* loop over list elements */
	{
	    nr = block_rows(m);
	    nc = block_cols(m);
	    
	    /* Form s.d. vector */
	    vector<Type> sd_list_el = exp(theta.segment(i1, nr));
	    sd_list(m) = sd_list_el;
	    i1 += nr; /* increment `theta` index */

	    /* Form correlation matrix */
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
	    
	    /* Form random effects block */
	    matrix<Type> block_list_el(nr, nc);
	    vector<Type> v(nr);
	    for (int j = 0; j < nc; j++) /* loop over block columns */
	    {
	        v = b.segment(i2, nr);
	        block_list_el.col(j) = v;
		b_scaled.segment(i2, nr) = sd_list_el * v;
		i2 += nr; /* increment `b` index */
	    }
	    block_list(m) = block_list_el;
	}
	if (!do_predict)
	{
	    REPORT(block_list);
	    REPORT(sd_list);
	    REPORT(cor_list);
	}
    }

    
    /* Compute response matrix ============================================== */

    /* Initialize (N, p) response matrix ------------------------------------ */

    matrix<Type> Y = Yo;
    
    /* Add fixed effects component ------------------------------------------ */

    /* Arrange elements of `beta` in (sum(beta_seg_len), p) matrix */
    Eigen::SparseMatrix<Type> beta_as_matrix(beta.size(), p);
    beta_as_matrix.reserve(beta_seg_len); /* declaring nnz */
    for (int i = 0; i < beta.size(); i++)
    {
        beta_as_matrix.insert(i, beta_seg_index(i)) = beta(i);
    }

    /* Increment response matrix with matrix product
       X * beta_as_matrix
       = (N, sum(beta_seg_len)) * (sum(beta_seg_len), p)
       = (N, p)
    */
    Y += (do_sparse_X ? Xs : Xd) * beta_as_matrix;
    
    /* Preserve response matrix without random effects component
       for simulation code
    */
    matrix<Type> Y_simulate = Y;

    /* Add random effects component ----------------------------------------- */

    /* (Declare outside of conditional scope for reuse by prediction code) */
    Eigen::SparseMatrix<Type> b_scaled_as_matrix(b.size(), p);

    if (any_random_effects)
    {
        /* Arrange elements of `b_scaled` in (sum(b_seg_len), p) matrix */
        b_scaled_as_matrix.reserve(b_seg_len); /* declaring nnz */
	for (int i = 0; i < b.size(); i++)
	{
	    b_scaled_as_matrix.insert(i, b_seg_index(i)) = b_scaled(i);
	}

	/* Increment response matrix with matrix product
	   Z * b_scaled_as_matrix
	   = (N, sum(b_seg_len)) * (sum(b_seg_len), p)
	   = (N, p)
	*/
	Y += Z * b_scaled_as_matrix;
    }
    if (!do_predict)
    {
        REPORT(Y);
        ADREPORT(Y);
    }

    
    /* Compute predictions ================================================== */

    /* Log cumulative incidence */
    vector<Type> log_curve = eval_log_curve(t, t_seg_len,
					    curve_flag, do_excess,
					    Y,
					    j_log_r, j_log_alpha, j_log_c0,
					    j_log_tinfl, j_log_K,
					    j_logit_p, j_log_a, j_log_b);

    /* Log interval incidence */
    vector<Type> log_cases = logspace_diff_n(log_curve, t_seg_len);
    if (do_day_of_week)
    {
        vector<Type> offsets =
	    get_day_of_week_offsets(x_seg_len, day1,
				    Y,
				    j_log_w1, j_log_w2, j_log_w3,
				    j_log_w4, j_log_w5, j_log_w6);
        log_cases += offsets;
    }


    /* Compute negative log likelihood ====================================== */

    /* Negative log likelihood */
    parallel_accumulator<Type> nll(this);
    Type nll_term;
    Type log_var_minus_mu;
    int width_s;
    int width_k;
    bool print_Y_row;

    if (do_trace)
    {
        if (trace_flag % 2 == 0)
	{
	    std::cout << "Y = \n" << Y << "\n";
	}
	printf("nll initialized to 0\ncommencing loop over observations\n");
        width_s = nchar(N - 1);
	width_k = nchar(t_seg_len.maxCoeff() - 1);
    }

    for (int s = 0, i = 0; s < N; s++) /* loop over segments */
    {
        print_Y_row = (trace_flag >= 3);
	
        for (int k = 0; k < x_seg_len(s); k++) /* loop over within-segment index */
	{
	    if (!is_NA_real_(x(i+k)))
	    {
	        switch (family_flag)
		{
		case pois:
		    nll_term = -dpois_robust(x(i+k), log_cases(i+k), true);
		    /* usage: dpois_robust(x, log_lambda, give_log) */
		    break;
		case nbinom:
		    log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
		    nll_term = -dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
		    /* usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log) */
		    break;
		}

		if (do_trace)
		{
		    if (trace_flag >= 3)
		    {
		        printf("at index %*d of segment %*d: nll term is %5.6e\n",
			       width_k, k, width_s, s, asDouble(nll_term));
		    }
		    else
		    {
		        if (!is_finite(nll_term))
			{
			    printf("at index %d of segment %d: nll term is non-finite\n", k, s);
			    print_Y_row = true;
			}
			else if (asDouble(nll_term) > 1.0e+09)
			{
			    printf("at index %d of segment %d: nll term exceeds 1.0e+09\n", k, s);
			    print_Y_row = true;
			}
		    }
		}

		nll += nll_term;
	    }
	}

	if (print_Y_row)
	{
	    std::cout << "Y.row(" << s << ") = " << Y.row(s) << "\n";
	}
	
	i += x_seg_len(s);
    }

    if (do_trace)
    {
        printf("loop over observations complete\nnll is %5.6e\n", asDouble(nll));
    }

    if (any_random_effects) {
        int width_m;
	int width_j;
	
	if (do_trace)
	{
	    printf("commencing loop over random effects\n");
	    width_m = nchar(M - 1);
	    width_j = nchar(block_cols.maxCoeff() - 1);
	}
      
        for (int m = 0; m < M; m++) /* loop over blocks */
    	{
	    density::MVNORM_t<Type> N_0_Sigma(cor_list(m)); /* function returning negative log density */
	    for (int j = 0; j < block_cols(m); j++) /* loop over block columns */
    	    {
	        nll_term = N_0_Sigma(block_list(m).col(j));

		if (do_trace)
		{
		    if (trace_flag >= 3)
		    {
		        printf("at column %*d of block %*d: nll term is %5.6e\n",
			       width_j, j, width_m, m, asDouble(nll_term));
		    }
		    else
		    {
			if (!is_finite(nll_term))
			{
			    printf("at column %*d of block %*d: nll term is non-finite\n",
				   width_j, j, width_m, m);
			}
			else if (asDouble(nll_term) > 1.0e+09)
			{
			    printf("at column %*d of block %*d: nll term exceeds 1.0e+09\n",
				   width_j, j, width_m, m);
			}
		    }
		}
		
		nll += nll_term;
    	    }
    	}

	if (do_trace)
	{
	    printf("loop over random effects complete\nnll is %5.6e\n", asDouble(nll));
	}
    }

    if (do_regularize)
    {
        vector<Type> hp;
        int width_i;
        int width_j;
      
        if (do_trace)
	{
	    printf("commencing loop over regularized parameters\n");
	    width_i = nchar(N - 1);
	    width_j = nchar(p - 1);
	}
	
        for (int j = 0; j < p; j++)
	{
	    if (regularize_flag(j) < 0)
	    {
	        continue;
	    }
	    hp = regularize_hyperpar(j);
	    for (int i = 0; i < N; i++)
	    {
	        switch (regularize_flag(j))
		{
		case norm:
		    nll_term = -dnorm(Y(i, j), hp(0), hp(1), true);
		    /* usage: dnorm(x, mean, sd, give_log) */
		    break;
		}
		
		if (do_trace)
		{
		    if (trace_flag >= 3)
		    {
		        printf("for parameter %*d in segment %*d: nll term is %5.6e\n",
			       width_j, j, width_i, i, asDouble(nll_term));
		    }
		    else
		    {
			if (!is_finite(nll_term))
			{
			    printf("for parameter %*d in segment %*d: nll term is non-finite\n",
				   width_j, j, width_i, i);
			}
			else if (asDouble(nll_term) > 1.0e+09)
			{
			    printf("for parameter %*d in segment %*d: nll term exceeds 1.0e+09\n",
				   width_j, j, width_i, i);
			}
		    }
		}

		nll += nll_term;
	    }
	}
	
	if (do_trace)
	{
	    printf("loop over regularized parameters complete\nnll is %5.6e\n", asDouble(nll));
	}
    }
    


    /* Simulate incidence =================================================== */

    SIMULATE
    {
        if (any_random_effects)
	{
	    vector<Type> b_simulate_scaled(b.size());
	    for (int m = 0, i = 0; m < M; m++) /* loop over blocks */
	    {
	        for (int j = 0; j < block_cols(m); j++) /* loop over block columns */
		{
		    b_simulate_scaled.segment(i, block_rows(m)) = sd_list(m) * density::MVNORM(cor_list(m)).simulate();
		    i += block_rows(m);
		}
	    }
	    Eigen::SparseMatrix<Type> b_simulate_scaled_as_matrix(b.size(), p);
	    b_simulate_scaled_as_matrix.reserve(b_seg_len); /* declaring nnz */
	    for (int i = 0; i < b.size(); i++)
	    {
		b_simulate_scaled_as_matrix.insert(i, b_seg_index(i)) = b_simulate_scaled(i);
	    }
	    Y_simulate += Z * b_simulate_scaled_as_matrix;
	}

        log_curve = eval_log_curve(t, t_seg_len,
				   curve_flag, do_excess,
				   Y_simulate,
				   j_log_r, j_log_alpha, j_log_c0,
				   j_log_tinfl, j_log_K,
				   j_logit_p, j_log_a, j_log_b);
	vector<Type> log_cases = logspace_diff_n(log_curve, t_seg_len);
	if (do_day_of_week)
	{
	    vector<Type> offsets =
		get_day_of_week_offsets(x_seg_len, day1,
					Y,
					j_log_w1, j_log_w2, j_log_w3,
					j_log_w4, j_log_w5, j_log_w6);
	    log_cases += offsets;
	}
	
	for (int s = 0, i = 0; s < N; s++) /* loop over segments */
	{
	    for (int k = 0; k < x_seg_len(s); k++) /* loop over within-segment index */
	    {
	        switch(family_flag)
		{
		case pois:
		    x(i+k) = rpois(exp(log_cases(i+k)));
		    /* usage: rpois(lambda) */
		    break;
		case nbinom:
		    x(i+k) = rnbinom_robust(log_cases(i+k), Y_simulate(s, j_log_nbdisp));
		    /* usage: rnbinom_robust(log_mu, log_size) */
		    break;
		}
	    }
	    i += x_seg_len(s);
	}
	REPORT(x);
    }


    /* Predict incidence ==================================================== */
    
    if (do_predict)
    {
        /* Set up =========================================================== */

        /* New data --------------------------------------------------------- */

        /* Concatenated time series segments
           - length=n'
	*/
        DATA_VECTOR(t_predict);

	/* Segment lengths
	   - length=N'
	*/
	DATA_IVECTOR(t_predict_seg_len);
	vector<int> x_predict_seg_len = decrement(t_predict_seg_len, 1);
	
	/* Segment initial days of week
	   - length=N' if do_day_of_week=true, otherwise length=0
	   - val={0,...,6} where 0=reference
	*/
        DATA_IVECTOR(day1_predict);

	/* Number of segments */
	int N_predict = t_predict_seg_len.size();

	/* Combined fixed effects model matrix
	   - dense format
	   - dim=(N', sum(beta_seg_len)) if do_sparse_X=false, otherwise dim=(N', 0)
	*/
	DATA_MATRIX(Xd_predict);
	/* 
	   - sparse format
	   - dim=(N', sum(beta_seg_len)) if do_sparse_X=true,  otherwise dim=(N', 0)
	*/
	DATA_SPARSE_MATRIX(Xs_predict);

	/* Combined random effects model matrix
	   - dim=(N', sum(b_seg_len))
	*/
	DATA_SPARSE_MATRIX(Z_predict);

	/* Offset component of response matrix `Y_predict`
	   - dim=(N', p)
	*/
	DATA_MATRIX(Yo_predict);

	/* New flags -------------------------------------------------------- */
	
	DATA_IVECTOR(what_flag); // predict  (1=yes, 0=no)
	bool do_predict_lii = (what_flag(0) == 1);
	bool do_predict_lci = (what_flag(1) == 1);
	bool do_predict_lrt = (what_flag(2) == 1);

	
	/* Compute response matrix ========================================== */
	
	matrix<Type> Y_predict = Yo_predict;
	Y_predict += (do_sparse_X ? Xs_predict : Xd_predict) * beta_as_matrix;
	if (any_random_effects)
    	{
	    Y_predict += Z_predict * b_scaled_as_matrix;
    	}

	
	/* Compute predictions ============================================== */

	/* Log cumulative incidence */
        vector<Type> log_curve_predict =
	    eval_log_curve(t_predict, t_predict_seg_len,
			   curve_flag, do_excess,
			   Y_predict, 
			   j_log_r, j_log_alpha, j_log_c0,
			   j_log_tinfl, j_log_K,
			   j_logit_p, j_log_a, j_log_b);

        /* Log interval incidence */
	vector<Type> log_cases_predict = logspace_diff_n(log_curve_predict, t_predict_seg_len);
	if (do_day_of_week)
	{
	    vector<Type> offsets =
	        get_day_of_week_offsets(x_predict_seg_len, day1_predict,
					Y_predict,
					j_log_w1, j_log_w2, j_log_w3,
					j_log_w4, j_log_w5, j_log_w6);
	    log_cases_predict += offsets;
	}

	if (do_predict_lii)
	{
	    /* Log interval incidence */
	    vector<Type> log_int_inc = log_cases_predict;
	    REPORT(log_int_inc);
	    ADREPORT(log_int_inc);
	}
	
	if (do_predict_lci)
	{
	    /* Log cumulative incidence
	       (since earliest time point in current segment, rather than t=-Inf)
	    */
	    int n;
	    vector<Type> log_cum_inc(log_cases_predict.size());
	    for (int s = 0, i = 0; s < N_predict; s++) /* loop over segments */
	    {
	        n = x_predict_seg_len(s);
	        vector<Type> v = log_cases_predict.segment(i, n);
		log_cum_inc.segment(i, n) = logspace_cumsum_1(v);
		i += n;
	    }
	    REPORT(log_cum_inc);
	    ADREPORT(log_cum_inc);
	}
	
	if (do_predict_lrt)
	{
	    /* Log per capita growth rate */
	    if (do_day_of_week)
	    {
	        /* Approximate rate from local linear regression */ 
	        vector<Type> log_rt = eval_log_rt_approx(log_cases_predict, x_predict_seg_len);
		REPORT(log_rt);
		ADREPORT(log_rt);
	    }
	    else
	    {
	        /* Exact rate from differential equation */ 
	        vector<Type> log_rt =
		    eval_log_rt_exact(t_predict, log_curve_predict, t_predict_seg_len,
				      curve_flag, do_excess,
				      Y_predict,
				      j_log_r, j_log_alpha, j_log_K,
				      j_logit_p, j_log_a, j_log_b);
		REPORT(log_rt);
		ADREPORT(log_rt);
	    }
	}
    }
    
    return nll;
}
