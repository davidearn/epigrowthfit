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
        for (int i = 0; i < LENGTH(x); ++i)
	{
            SEXP v = VECTOR_ELT(x, i);
            (*this)(i) = asVector<Type>(v);
        }
    }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
    /* Parameters =========================================================== */
  
    /* Concatenated fixed effects coefficient vectors */
    PARAMETER_VECTOR(beta);

    /* Concatenated random effects coefficient vectors (unit variance scale) */
    PARAMETER_VECTOR(b);

    /* Concatenated random effects covariance parameters */
    PARAMETER_VECTOR(theta);
    

    /* Data ================================================================= */

    /* Concatenated time series segments 
       - length=n
    */
    DATA_VECTOR(time);
    /* 
       - length=n-N 
    */
    DATA_VECTOR(x);

    /* Segment lengths
       - length=N
    */
    DATA_IVECTOR(time_seg_len);
    vector<int> x_seg_len = decrement(time_seg_len);
    
    /* Number of segments */
    int N = time_seg_len.size();

    /* Segment initial days of week
       - length=N if day_of_week=true, length=0 otherwise
       - val={0,...,6} where 0=reference
    */
    DATA_IVECTOR(day1);

    /* Offset component of response matrix `Y`
       - dim=(N, p)
    */
    DATA_MATRIX(Yo);

    /* Number of nonlinear model parameters */
    int p = Yo.cols();
    
    /* Column indices of nonlinear model parameters in response matrix `Y`
       - val={0,...p-1} if parameter is used, val=-1 otherwise
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
       - dim=(N, length(beta)) if do_sparse_X=false, dim=(N, 0) otherwise
    */
    DATA_MATRIX(Xd);
    /* 
       - sparse format
       - dim=(N, length(beta)) if do_sparse_X=true,  dim=(N, 0) otherwise
    */
    DATA_SPARSE_MATRIX(Xs);

    /* Combined random effects model matrix
       - sparse format
       - dim=(N, length(b))
    */
    DATA_SPARSE_MATRIX(Z);

    /* Indicator for random effects model */
    bool any_random_effects = (Z.cols() > 0);

    /* "Factors" splitting model matrix columns by relation 
       to a nonlinear model parameter
       - length=length(beta)
       - val={0,...,p-1}
    */
    DATA_IVECTOR(beta_index);
    /*
       - length=length(b)
       - val={0,...,p-1}
    */
    DATA_IVECTOR(b_index);

    /* Number of model matrix columns by nonlinear model parameter
       - in R: *_index_tab = c(table(*_index))
       - length=p
    */
    DATA_IVECTOR(beta_index_tab);
    DATA_IVECTOR(b_index_tab);
    
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
    

    /* Flags ================================================================ */
    
    /* Model of cumulative incidence (enum.h) */
    DATA_INTEGER(curve_flag);
    /* Baseline term in model of cumulative incidence (1=yes, 0=no) */
    DATA_INTEGER(excess_flag);
    /* Model of observation error (enum.h) */
    DATA_INTEGER(family_flag);
    /* Day of week effects (1=yes, 0=no) */
    DATA_INTEGER(day_of_week_flag);
    /* Trace
       0=nothing
       1=degenerate nll terms
       2=all nll terms
    */
    DATA_INTEGER(trace_flag);
    /* X format (1=sparse, 0=dense) */
    DATA_INTEGER(sparse_X_flag);
    /* predict (1=yes, 0=no) */
    DATA_INTEGER(predict_flag);
    /* Priors on nonlinear model parameters (enum.h)
       - length=p
       - val=-1 if no prior
    */
    DATA_IVECTOR(regularize_flag);

    bool do_excess        = (excess_flag      == 1);
    bool do_day_of_week   = (day_of_week_flag == 1);
    bool do_trace         = (trace_flag       >= 1);
    bool do_verbose_trace = (trace_flag       >= 2);
    bool do_sparse_X      = (sparse_X_flag    == 1);
    bool do_predict       = (predict_flag     == 1);
    bool do_regularize    = false;
    for (int i = 0; !do_regularize && i < p; ++i)
    {
        do_regularize = regularize_flag(i) >= 0;
    }
    
    
    /* Random effects infrastructure ======================================== */

    /* List of matrices in which to arrange elements of `b`:
       - blocks, block rows, and block columns correspond to random effect 
         terms, nonlinear model parameters, and group levels, respectively
       - block columns represent realizations of a multivariate normal random 
         variable with zero mean, unit variance, and unstructured covariance
    */
    vector< matrix<Type> > block_list(M);

    /* List of corresponding vectors of standard deviations */
    vector< vector<Type> > sd_list(M);

    /* List of corresponding vectors of other covariance parameters */
    vector< vector<Type> > chol_list(M);

    /* List of corresponding negative log density functions */
    vector< density::MVNORM_t<Type> > nld_list(M);

    if (any_random_effects)
    {
        /* Random vectors of length 1 are handled specially */
        vector<Type> chol1(0);
        matrix<Type> cor1(1, 1);
	cor1(0, 0) = Type(1.0);
	density::MVNORM_t<Type> nld1(cor1);
      
        int nr;
	int nc;
        for (int m = 0, i1 = 0, i2 = 0; m < M; ++m)
	{ /* loop over list elements */
	    nr = block_rows(m);
	    nc = block_cols(m);

	    /* Vector of standard deviations */
	    vector<Type> sd = exp(theta.segment(i1, nr));
	    sd_list(m) = sd;
	    i1 += nr;
	    
	    if (nr == 1)
	    {
	        /* Vector of other covariance parameters */
	        chol_list(m) = chol1;

	        /* Negative log density function */
	        nld_list(m) = nld1;
	    }
	    else
	    {
	        /* Vector of other covariance parameters */
	        vector<Type> chol = theta.segment(i1, nr * (nr - 1) / 2);
		chol_list(m) = chol;
		i1 += chol.size();

		/* Negative log density function */
		density::UNSTRUCTURED_CORR_t<Type> nld(chol);
		nld_list(m) = nld;
	    }
	    
	    /* Random effects block */
	    matrix<Type> block = b.segment(i2, nr * nc).matrix();
	    block.resize(nr, nc);
	    block_list(m) = block;

	    /* Scale `b` in place */
	    for (int j = 0; j < nc; ++j)
	    {
	        b.segment(i2, nr) *= sd;
	        i2 += nr;
	    }
	    block_list(m) = block;
	} /* loop over list elements */

	if (!do_predict)
	{
	    REPORT(block_list);
	    REPORT(sd_list);
	    REPORT(chol_list);
	}
    }


    /* Coefficient matrices ================================================= */

    /* Arrange elements of `beta` in (length(beta), p) matrix */
    Eigen::SparseMatrix<Type> beta_matrix =
        as_sparse_matrix(beta, beta_index, beta_index_tab);

    /* Arrange elements of `b` in (length(b), p) matrix */
    Eigen::SparseMatrix<Type> b_matrix =
        as_sparse_matrix(b, b_index, b_index_tab);
    
    
    /* Response matrix ====================================================== */

    /* Initialize (N, p) response matrix */
    matrix<Type> Y = Yo;
    
    /* Add fixed effects component */
    Y += (do_sparse_X ? Xs : Xd) * beta_matrix;
    
    /* Preserve value without random effects component for simulation code */
    matrix<Type> simulate_Y = Y;

    /* Add random effects component */
    if (any_random_effects)
    {
        Y += Z * b_matrix;
    }
    
    if (!do_predict)
    {
        REPORT(Y);
        ADREPORT(Y);
    }

    
    /* Predicted incidence ================================================== */

    /* Log cumulative incidence (since time -Inf) */
    vector<Type> log_curve = eval_log_curve(time, time_seg_len,
					    curve_flag, do_excess,
					    Y,
					    j_log_r, j_log_alpha, j_log_c0,
					    j_log_tinfl, j_log_K,
					    j_logit_p, j_log_a, j_log_b);

    /* Log interval incidence */
    vector<Type> log_cases = logspace_diff(log_curve, time_seg_len);
    if (do_day_of_week)
    {
        vector<Type> offsets =
	    get_day_of_week_offsets(x_seg_len, day1,
				    Y,
				    j_log_w1, j_log_w2, j_log_w3,
				    j_log_w4, j_log_w5, j_log_w6);
        log_cases += offsets;
    }


    /* Negative log likelihood ============================================== */

    parallel_accumulator<Type> nll(this);
    Type nll_term;
    
    
    /* Observation likelihood ----------------------------------------------- */

    Type log_var_minus_mu;
    int x0;
    double lambda;
    double mu;
    double size;
    int width_k;
    int width_s;
    bool print_Y_row;
    
    if (do_trace)
    {
        std::cout << "Y = \n" << Y << "\n";
	printf("nll initialized to %.6e\n", asDouble(nll));
	std::cout << "commencing loop over observations\n";
	width_k = nchar(x_seg_len.maxCoeff() - 1);
	width_s = nchar(N - 1);
    }

    for (int s = 0, i = 0; s < N; i += x_seg_len(s), ++s)
    { /* loop over segments */
        print_Y_row = do_verbose_trace;
	
        for (int k = 0; k < x_seg_len(s); ++k)
	{ /* loop over within-segment index */
	    if (is_na(x(i+k)))
	    {
	        continue;
	    }
	    
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
	    nll += nll_term;

	    if (do_trace && (do_verbose_trace || !is_nll_term_ok(nll_term)))
	    {
	        print_Y_row = true;

		printf("at index %*d of segment %*d: nll term is %.6e\n",
		       width_k, k, width_s, s, asDouble(nll_term));
		x0 = asDouble(x(i+k));
		switch(family_flag)
		{
		case pois:
		    lambda = asDouble(exp(log_cases(i+k)));
		    printf("  \\__ -log(dpois(x = %d, lambda = %.6e))\n", x0, lambda);
		    break;
		case nbinom:
		    mu = asDouble(exp(log_cases(i+k)));
		    size = asDouble(exp(Y(s, j_log_nbdisp)));
		    printf("  \\__ -log(dnbinom(x = %d, mu = %.6e, size = %.6e))\n", x0, mu, size);
		    break;
		}
	    }
	} /* loop over within-segment index */

	if (print_Y_row)
	{
	    std::cout << "Y.row(" << s << ") = " << Y.row(s) << "\n";
	}
    } /* loop over segments */

    if (do_trace)
    {
        std::cout << "loop over observations complete\n";
	printf("nll is %.6e\n", asDouble(nll));
    }


    /* Random effect likelihood --------------------------------------------- */

    if (any_random_effects) {

        int width_j;
        int width_m;

        if (do_trace)
	{
	    std::cout << "commencing loop over random effects\n";
	    width_j = nchar(block_cols.maxCoeff() - 1);
	    width_m = nchar(M - 1);
	}
      
        for (int m = 0; m < M; ++m)
    	{ /* loop over blocks */
	    for (int j = 0; j < block_cols(m); ++j)
    	    { /* loop over block columns */
	        nll_term = nld_list(m)(block_list(m).col(j));
		nll += nll_term;

		if (do_trace && (do_verbose_trace || !is_nll_term_ok(nll_term)))
		{
		    printf("at column %*d of block %*d: nll term is %.6e\n",
			   width_j, j, width_m, m, asDouble(nll_term));
		}
	    } /* loop over block columns */
	} /* loop over blocks */

	if (do_trace)
	{
	    std::cout << "loop over random effects complete\n";
	    printf("nll is %.6e\n", asDouble(nll));
	}
    }


    /* Nonlinear model parameter likelihood --------------------------------- */

    if (do_regularize)
    {
        double x0;
        double mu;
	double sigma;
        int width_i;
        int width_j;
      
        if (do_trace)
	{
	    std::cout << "commencing loop over regularized parameters\n";
	    width_i = nchar(N - 1);
	    width_j = nchar(p - 1);
	}
	
        for (int j = 0; j < p; ++j)
	{ /* loop over parameters */
	    if (regularize_flag(j) < 0)
	    {
	        continue;
	    }
	    
	    vector<Type> hp = regularize_hyperpar(j);
	    for (int i = 0; i < N; ++i)
	    { /* loop over values */
	        switch (regularize_flag(j))
		{
		case norm:
		    nll_term = -dnorm(Y(i, j), hp(0), hp(1), true);
		    /* usage: dnorm(x, mean, sd, give_log) */
		    break;
		}
		nll += nll_term;

		if (do_trace && (do_verbose_trace || !is_nll_term_ok(nll_term)))
		{
		    printf("parameter %*d in segment %*d: nll term is %.6e\n",
			   width_j, j, width_i, i, asDouble(nll_term));
		    x0 = asDouble(Y(i, j));
		    switch (regularize_flag(j))
		    {
		    case norm:
		        mu = asDouble(hp(0));
			sigma = asDouble(hp(1));
			printf("  \\__ -log(dnorm(x = %.6e, mean = %.6e, sd = %.6e))\n", x0, mu, sigma);
			break;
		    }
		}
	    } /* loop over values */
	} /* loop over parameters */
	
	if (do_trace)
	{
	    std::cout << "loop over regularized parameters complete\n";
	    printf("nll is %.6e\n", asDouble(nll));
	}
    }


    /* Prediction =========================================================== */
    
    if (do_predict)
    {
        /* New data --------------------------------------------------------- */

        /* Concatenated time series segments
           - length=n'
	*/
        DATA_VECTOR(predict_time);

	/* Segment lengths
	   - length=N'
	*/
	DATA_IVECTOR(predict_time_seg_len);
	vector<int> predict_x_seg_len = decrement(predict_time_seg_len);

	/* Segment initial days of week
	   - length=N' if do_day_of_week=true, length=0 otherwise
	   - val={0,...,6} where 0=reference
	*/
        DATA_IVECTOR(predict_day1);

	/* Combined fixed effects model matrix
	   - dense format
	   - dim=(N', length(beta)) if do_sparse_X=false, dim=(N', 0) otherwise
	*/
	DATA_MATRIX(predict_Xd);
	/* 
	   - sparse format
	   - dim=(N', length(beta)) if do_sparse_X=true,  dim=(N', 0) otherwise
	*/
	DATA_SPARSE_MATRIX(predict_Xs);

	/* Combined random effects model matrix
	   - dim=(N', length(b))
	*/
	DATA_SPARSE_MATRIX(predict_Z);

	/* Offset component of response matrix `predict_Y`
	   - dim=(N', p)
	*/
	DATA_MATRIX(predict_Yo);

	
	/* New flags -------------------------------------------------------- */

	/* predict (1=yes, 0=no) */
	DATA_IVECTOR(what_flag);
	bool do_predict_lii = (what_flag(0) == 1);
	bool do_predict_lci = (what_flag(1) == 1);
	bool do_predict_lrt = (what_flag(2) == 1);

	
	/* Response matrix -------------------------------------------------- */
	
	matrix<Type> predict_Y = predict_Yo;
	predict_Y += (do_sparse_X ? predict_Xs : predict_Xd) * beta_matrix;
	if (any_random_effects)
    	{
	    predict_Y += predict_Z * b_matrix;
    	}

	
	/* Predicted incidence ---------------------------------------------- */

	/* Log cumulative incidence (since time -Inf) */
        vector<Type> predict_log_curve =
	    eval_log_curve(predict_time, predict_time_seg_len,
			   curve_flag, do_excess,
			   predict_Y, 
			   j_log_r, j_log_alpha, j_log_c0,
			   j_log_tinfl, j_log_K,
			   j_logit_p, j_log_a, j_log_b);

        /* Log interval incidence */
	vector<Type> predict_log_cases =
	    logspace_diff(predict_log_curve, predict_time_seg_len);
	if (do_day_of_week)
	{
	    vector<Type> offsets =
	        get_day_of_week_offsets(predict_x_seg_len, predict_day1,
				        predict_Y,
					j_log_w1, j_log_w2, j_log_w3,
					j_log_w4, j_log_w5, j_log_w6);
	    predict_log_cases += offsets;
	}


	/* Reported variables ----------------------------------------------- */

	if (do_predict_lii)
	{
	    /* Log interval incidence */
	    vector<Type> log_int_inc = predict_log_cases;
	    REPORT(log_int_inc);
	    ADREPORT(log_int_inc);
	}
	
	if (do_predict_lci)
	{
	    /* Log cumulative incidence (since time 0) */
	    vector<Type> log_cum_inc =
	        logspace_cumsum(predict_log_cases, predict_x_seg_len);
	    REPORT(log_cum_inc);
	    ADREPORT(log_cum_inc);
	}
	
	if (do_predict_lrt)
	{
	    /* Log per capita growth rate */
	    if (do_day_of_week)
	    {
	        /* Approximate rate from local linear regression */ 
	        vector<Type> log_rt =
		    eval_log_rt_approx(predict_log_cases, predict_x_seg_len);
		REPORT(log_rt);
		ADREPORT(log_rt);
	    }
	    else
	    {
	        /* Exact rate from differential equation */ 
	        vector<Type> log_rt =
		    eval_log_rt_exact(predict_time, predict_log_curve,
				      predict_time_seg_len,
				      curve_flag, do_excess,
				      predict_Y,
				      j_log_r, j_log_alpha, j_log_K,
				      j_logit_p, j_log_a, j_log_b);
		REPORT(log_rt);
		ADREPORT(log_rt);
	    }
	}
    }


    /* Simulation =========================================================== */

    SIMULATE
    {
        /* Response matrix -------------------------------------------------- */

        if (any_random_effects)
	{
	    /* Update `b` in place with simulated random effects */
	    int nr;
	    int nc;
	    for (int m = 0, i = 0; m < M; ++m)
	    {
	        nr = block_rows(m);
	        nc = block_cols(m);
	        for (int j = 0; j < nc; ++j)
		{
		    b.segment(i, nr) = sd_list(m) * nld_list(m).simulate();
		    i += nr;
		}
	    }

	    /* Update `b_matrix` in place */
	    for (int i = 0; i < b.size(); ++i)
	    {
		b_matrix.coeffRef(i, b_index(i)) = b(i);
	    }

	    /* Increment incomplete response matrix `simulate_Y`
	       with simulated random effects component
	    */
	    simulate_Y += Z * b_matrix;
	}


	/* Predicted incidence ---------------------------------------------- */
	
	/* Log cumulative incidence (since time -Inf) */
        log_curve = eval_log_curve(time, time_seg_len,
				   curve_flag, do_excess,
				   simulate_Y,
				   j_log_r, j_log_alpha, j_log_c0,
				   j_log_tinfl, j_log_K,
				   j_logit_p, j_log_a, j_log_b);

	/* Log interval incidence */
	log_cases = logspace_diff(log_curve, time_seg_len);
	if (do_day_of_week)
	{
	    vector<Type> offsets =
		get_day_of_week_offsets(x_seg_len, day1,
					Y,
					j_log_w1, j_log_w2, j_log_w3,
					j_log_w4, j_log_w5, j_log_w6);
	    log_cases += offsets;
	}


	/* Simulated time series -------------------------------------------- */

	/* Update `x` in place with simulated observations */
	int n;
	for (int s = 0, i = 0; s < N; ++s)
	{
	    n = x_seg_len(s);
	    for (int k = 0; k < n; ++k)
	    {
	        switch(family_flag)
		{
		case pois:
		    x(i+k) = rpois(exp(log_cases(i+k)));
		    /* usage: rpois(lambda) */
		    break;
		case nbinom:
		    x(i+k) = rnbinom_robust(log_cases(i+k), simulate_Y(s, j_log_nbdisp));
		    /* usage: rnbinom_robust(log_mu, log_size) */
		    break;
		}
	    }
	    i += n;
	}
	REPORT(x);
    }

    return nll;
}
