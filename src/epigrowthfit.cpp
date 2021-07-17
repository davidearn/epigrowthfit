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
    
    /* Structure storing integer column indices of nonlinear model parameters 
       in response matrix
       - val={0,...p-1} if parameter is used, val=-1 otherwise
    */
    DATA_STRUCT(indices, egf::indices_t);

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
    DATA_STRUCT(hyperparameters, egf::list_of_vectors_t);

    /* Structure storing integer and Boolean flags specifying model
       and template behaviour
    */
    DATA_STRUCT(flags, egf::flags_t);
    flags.do_simulate = this->do_simulate;
    flags.do_random_effects = (Z.cols() > 0);
    
    
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

    if (flags.do_random_effects)
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
	        u = (flags.do_simulate) ? list_of_nld(m).simulate() : b.segment(i2, nr);
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
    Y += (flags.do_sparse_X ? Xs : Xd) * beta_matrix;

    /* Add random effects component */
    if (flags.do_random_effects)
    {
        Y += Z * b_matrix;
    }
    
    REPORT(Y);
    if (!flags.do_predict)
    {
        ADREPORT(Y);
    }


    /* Negative log likelihood ============================================== */

    Type nll = Type(0.0);
    
    /* Observation likelihood ----------------------------------------------- */

    if (flags.do_trace && !flags.do_simulate)
    {
        std::cout << "Y = \n" << Y << "\n";
	printf("nll initialized to %.6e\n", asDouble(nll));
	std::cout << "commencing loop over observations\n";
    }
    add_nll_ob(nll, this, time, time_seg_len, x, Y, indices, flags, day1);
    if (flags.do_simulate)
    {
        REPORT(x);
	return nll;
    }
    if (flags.do_trace)
    {
        std::cout << "loop over observations complete\n";
        printf("nll is %.6e\n", asDouble(nll));
    }
    
    /* Random effect likelihood --------------------------------------------- */

    if (flags.do_random_effects) {
        if (flags.do_trace)
	{
	    std::cout << "commencing loop over random effects\n";
	}
	add_nll_re(nll, this, list_of_blocks, list_of_nld, flags);
	if (flags.do_trace)
	{
	    std::cout << "loop over random effects complete\n";
	    printf("nll is %.6e\n", asDouble(nll));
	}
    }

    /* Parameter value likelihood ------------------------------------------- */

    if (flags.do_regularize)
    {
        if (flags.do_trace)
	{
	    std::cout << "commencing loop over regularized parameters\n";
	}
        add_nll_pv(nll, this, Y, hyperparameters, flags);
        if (flags.do_trace)
	{
	    std::cout << "loop over regularized parameters complete\n";
	    printf("nll is %.6e\n", asDouble(nll));
	}
    }


    /* Prediction =========================================================== */
    
    if (flags.do_predict)
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

	/* Flags */
	DATA_IVECTOR(flag_what);
	flags.do_predict_lii = (flag_what(0) == 1);
	flags.do_predict_lci = (flag_what(1) == 1);
	flags.do_predict_lrt = (flag_what(2) == 1);

	
	/* Response matrix -------------------------------------------------- */
	
	predict_Y += (flags.do_sparse_X ? predict_Xs : predict_Xd) * beta_matrix;
	if (flags.do_random_effects)
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



    return nll;
}
