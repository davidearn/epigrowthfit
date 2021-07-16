#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "utils.h"
#include "enum.h"
#include "curve.h"
#include "nll.h"


template<class Type>
Type objective_function<Type>::operator() ()
{
#ifdef _OPENMP
    this->max_parallel_regions = omp_get_max_threads();
#endif
  
    /* Parameters =========================================================== */
  
    /* Concatenated vectors of fixed effects coefficients */
    PARAMETER_VECTOR(beta);

    /* Concatenated vectors of random effects coefficients (unit variance) */
    PARAMETER_VECTOR(b);

    /* Concatenated vectors of random effects covariance parameters */
    PARAMETER_VECTOR(theta);
    

    /* Data ================================================================= */

    /* Concatenated vectors of time points 
       - length=n
    */
    DATA_VECTOR(time);
    
    /* Concatenated vectors of observed counts 
       - length=n-N 
    */
    DATA_VECTOR(x);

    /* Segment lengths
       - length=N
    */
    DATA_IVECTOR(time_seg_len);

    /* Number of segments */
    int N = time_seg_len.size();
    
    /* Segment initial days of week
       - length=N if do_day_of_week=true, length=0 otherwise
       - val={0,...,6} where 0=reference
    */
    DATA_IVECTOR(day1);

    /* Offset component of response matrix
       - dim=(N, p)
    */
    DATA_MATRIX(Y);

    /* Number of nonlinear model parameters */
    int p = Yo.cols();
    
    /* Column indices of nonlinear model parameters in response matrix
       - val={0,...p-1} if parameter is used, val=-1 otherwise
    */
    DATA_INTEGER(index_log_r);
    DATA_INTEGER(index_log_alpha);
    DATA_INTEGER(index_log_c0);
    DATA_INTEGER(index_log_tinfl);
    DATA_INTEGER(index_log_K);
    DATA_INTEGER(index_logit_p);
    DATA_INTEGER(index_log_a);
    DATA_INTEGER(index_log_b);
    DATA_INTEGER(index_log_nbdisp);
    DATA_INTEGER(index_log_w1);
    DATA_INTEGER(index_log_w2);
    DATA_INTEGER(index_log_w3);
    DATA_INTEGER(index_log_w4);
    DATA_INTEGER(index_log_w5);
    DATA_INTEGER(index_log_w6);

    /* Combined fixed effects model matrix
       - dim=(N, length(beta)) if do_sparse_X=false, dim=(N, 0) otherwise
    */
    DATA_MATRIX(Xd);
    /* 
       - dim=(N, length(beta)) if do_sparse_X=true,  dim=(N, 0) otherwise
    */
    DATA_SPARSE_MATRIX(Xs);

    /* Combined random effects model matrix
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
    DATA_STRUCT(regularize_hyperpar, egf::list_of_vectors_t);
    

    /* Flags ================================================================ */
    
    /* Model of cumulative incidence (enum.h) */
    DATA_INTEGER(flag_curve);
    /* Baseline term in model of cumulative incidence (1=yes, 0=no) */
    DATA_INTEGER(flag_excess);
    /* Model of observation error (enum.h) */
    DATA_INTEGER(flag_family);
    /* Day of week effects (1=yes, 0=no) */
    DATA_INTEGER(flag_day_of_week);
    /* Priors on nonlinear model parameters (enum.h)
       - length=p
       - val=-1 if no prior
    */
    DATA_IVECTOR(flag_regularize);
    /* Trace
       0=nothing
       1=degenerate nll terms
       2=all nll terms
    */
    DATA_INTEGER(flag_trace);
    /* X format (1=sparse, 0=dense) */
    DATA_INTEGER(flag_sparse_X);
    /* predict (1=yes, 0=no) */
    DATA_INTEGER(flag_predict);

    bool do_excess        = (flag_excess      == 1);
    bool do_day_of_week   = (flag_day_of_week == 1);
    bool do_trace         = (flag_trace       >= 1);
    bool do_trace_verbose = (flag_trace       >= 2);
    bool do_sparse_X      = (flag_sparse_X    == 1);
    bool do_predict       = (flag_predict     == 1);
    bool do_simulate      = this->do_simulate;
    bool do_regularize    = false;
    for (int i = 0; !do_regularize && i < p; ++i)
    {
        do_regularize = flag_regularize(i) >= 0;
    }
    

    /* Coefficient matrices ================================================= */

    /* Fixed effects, to be multiplied on the left by `X` */
    Eigen::SparseMatrix<Type> beta_matrix(beta.size(), p);
    beta_matrix.reserve(beta_index_tab);
    for (int i = 0; i < beta.size(); ++i)
    {
        beta_matrix.insert(i, beta_index(i)) = beta(i);
    }

    /* Random effects, to be multiplied on the left by `Z` */
    Eigen::SparseMatrix<Type> b_matrix;
    
    
    /* Random effects infrastructure ======================================== */

    /* List of matrices in which to arrange elements of `b`:
       - blocks, block rows, and block columns correspond to random effect
         terms, nonlinear model parameters, and group levels, respectively
       - block columns represent realizations of a multivariate normal random 
         variable with zero mean, unit variance, and unstructured covariance
    */
    vector< matrix<Type> > list_of_blocks(M);

    /* List of corresponding vectors of standard deviations */
    vector< vector<Type> > list_of_sd(M);

    /* List of corresponding vectors of other covariance parameters */
    vector< vector<Type> > list_of_chol(M);

    /* List of corresponding negative log density functions */
    vector< density::MVNORM_t<Type> > list_of_nld(M);

    if (any_random_effects)
    {
        b_matrix.resize(b.size(), p);
        b_matrix.reserve(b_index_tab);
	
        matrix<Type> cor1(1, 1);
	cor1(0, 0) = Type(1.0);
	density::MVNORM_t<Type> nld1(cor1);

	int nr;
	int nc;
	int nt;
	vector<Type> u;
	for (int m = 0, i1 = 0, i2 = 0; m < M; ++m)
	{ /* loop over random effects terms */
	    nr = block_rows(m);
	    nc = block_cols(m);
	    nt = nr * (nr - 1) / 2;
	    
	    /* Vector of standard deviations */
	    list_of_sd(m) = exp(theta.segment(i1, nr));
	    i1 += nr;
	    
	    if (nr == 1)
	    {
	        /* Negative log density function */
	        list_of_nld(m) = nld1;
	    }
	    else
	    {
	        /* Vector of other covariance parameters */
	        list_of_chol(m) = theta.segment(i1, nt);
		i1 += nt;

		/* Negative log density function */
		list_of_nld(m) = density::UNSTRUCTURED_CORR(list_of_chol(m));
	    }

	    list_of_blocks(m).resize(nr, nc);
	    for (int j = 0; j < nc; ++j)
	    { /* loop over group levels */
	        /* Random effects block column */
	        u = (do_simulate) ? list_of_nld(m).simulate() : b.segment(i2, nr);
		list_of_blocks(m).col(j) = u;

		/* Insertion of scaled elements in coefficient matrix */
		for (int k = 0; k < nr; ++k)
		{
		    b_matrix.insert(i2, b_index(i2)) = list_of_sd(m)(k) * u(k);
		    ++i2;
		}
	    } /* loop over group levels */
	} /* loop over random effects terms */

	REPORT(list_of_blocks);
	REPORT(list_of_sd);
	REPORT(list_of_chol);
    }


    /* Response matrix ====================================================== */

    /* Add fixed effects component */
    Y += (do_sparse_X ? Xs : Xd) * beta_matrix;

    /* Add random effects component */
    if (any_random_effects)
    {
        Y += Z * b_matrix;
    }
    
    REPORT(Y);
    if (!do_predict)
    {
        ADREPORT(Y);
    }


    /* Negative log likelihood ============================================== */

    Type nll = 0.0;    
    if (do_trace)
    {
        std::cout << "Y = \n" << Y << "\n";
	printf("nll initialized to %.6e\n", asDouble(nll));
    }
    
    /* Observation likelihood ----------------------------------------------- */

    nll += egf::get_nll_observations(time,
				     x,
				     len,
				     Y,
				     flag_curve,
				     index_log_r,
				     index_log_alpha,
				     index_log_c0,
				     index_log_tinfl,
				     index_log_K,
				     index_logit_p,
				     index_log_a,
				     do_excess,
				     index_log_b,
				     do_day_of_week,
				     index_log_w1,
				     index_log_w2,
				     index_log_w3,
				     index_log_w4,
				     index_log_w5,
				     index_log_w6,
				     day1,
				     flag_family,
				     index_log_nbdisp,
				     do_trace,
				     do_trace_verbose);
    if (do_trace)
    {
        printf("nll is %.6e\n", asDouble(nll));
    }

    /* Random effect likelihood --------------------------------------------- */

    if (any_random_effects) {

        nll += egf::get_nll_random_effects(list_of_blocks,
					   list_of_nld,
					   do_trace,
					   do_trace_verbose);
	if (do_trace)
	{
	    printf("nll is %.6e\n", asDouble(nll));
	}
    }

    /* Nonlinear model parameter likelihood --------------------------------- */

    if (do_regularize)
    {
        nll += egf::get_nll_parameter_values(Y,
					     regularize_hyperpar,
					     do_trace,
					     do_trace_verbose);
        if (do_trace)
	{
	    printf("nll is %.6e\n", asDouble(nll));
	}
    }


    /* Prediction =========================================================== */
    
    if (do_predict)
    {
        /* New data --------------------------------------------------------- */

        /* Concatenated vectors of time points
           - length=n'
	*/
        DATA_VECTOR(predict_time);

	/* Segment lengths
	   - length=N'
	*/
	DATA_IVECTOR(predict_time_seg_len);
	int predict_N = predict_time_seg_len.size();

	/* Segment initial days of week
	   - length=N' if do_day_of_week=true, length=0 otherwise
	   - val={0,...,6} where 0=reference
	*/
        DATA_IVECTOR(predict_day1);

	/* Combined fixed effects model matrix
	   - dim=(N', length(beta)) if do_sparse_X=false, dim=(N', 0) otherwise
	*/
	DATA_MATRIX(predict_Xd);
	/* 
	   - dim=(N', length(beta)) if do_sparse_X=true,  dim=(N', 0) otherwise
	*/
	DATA_SPARSE_MATRIX(predict_Xs);

	/* Combined random effects model matrix
	   - dim=(N', length(b))
	*/
	DATA_SPARSE_MATRIX(predict_Z);

	/* Offset component of response matrix
	   - dim=(N', p)
	*/
	DATA_MATRIX(predict_Y);

	
	/* New flags -------------------------------------------------------- */

	/* predict (1=yes, 0=no) */
	DATA_IVECTOR(what_flag);
	bool do_predict_lii = (what_flag(0) == 1);
	bool do_predict_lci = (what_flag(1) == 1);
	bool do_predict_lrt = (what_flag(2) == 1);

	
	/* Response matrix -------------------------------------------------- */
	
	predict_Y += (do_sparse_X ? predict_Xs : predict_Xd) * beta_matrix;
	if (any_random_effects)
    	{
	    predict_Y += predict_Z * b_matrix;
    	}

	
	/* Predicted incidence ---------------------------------------------- */

	/* Log cumulative incidence since time -Inf */
	vector<Type> foo(predict_time.size());
	vector<Type> bar;
	int n;
	for (int s = 0, i = 0; s < predict_N; ++s)
	{
	    n = predict_time_seg_len(s);
	    bar = predict_time.segment(i, n);

	    egf::eval_log_curve(bar,
				Y_row,
				flag_curve,
				index_log_r,
				index_log_alpha,
				index_log_c0,
				index_log_tinfl,
				index_log_K,
				index_logit_p,
				index_log_a);
	    
	}
	

	
	vector<Type> predict_log_diff_curve(predict_time.size() - predict_N);
	vector<Type> predict_log_diff_curve(
	int n;
	for (int s = 0, i = 0; s < N_predict; ++s)
	{
	    n = predict_time_seg_len(s);
	    
	}

	
	
	
        vector<Type> predict_log_curve =
	    eval_log_curve(predict_time, predict_time_seg_len,
			   curve_flag, do_excess,
			   predict_Y, 
			   index_log_r, index_log_alpha, index_log_c0,
			   index_log_tinfl, index_log_K,
			   index_logit_p, index_log_a, index_log_b);

        /* Log interval incidence */
	vector<Type> predict_log_cases =
	    logspace_diff(predict_log_curve, predict_time_seg_len);
	if (do_day_of_week)
	{
	    vector<Type> offsets =
	        get_day_of_week_offsets(predict_x_seg_len, predict_day1,
				        predict_Y,
					index_log_w1, index_log_w2, index_log_w3,
					index_log_w4, index_log_w5, index_log_w6);
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
				      index_log_r, index_log_alpha, index_log_K,
				      index_logit_p, index_log_a, index_log_b);
		REPORT(log_rt);
		ADREPORT(log_rt);
	    }
	}
    }


    /* Simulation =========================================================== */

    SIMULATE
    {
        


	/* Predicted incidence ---------------------------------------------- */
	
	


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
		    x(i+k) = rnbinom_robust(log_cases(i+k), simulate_Y(s, index_log_nbdisp));
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
