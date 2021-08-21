#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>
#include "enums.h"
#include "structs.h"
#include "utils.h"
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

    /* Concatenated vectors of random effects covariance parameters */
    PARAMETER_VECTOR(theta);

    /* Concatenated vectors of random effects coefficients (unit variance) */
    PARAMETER_VECTOR(b);
    

    /* Flags ================================================================ */

    /* Structure storing integer and Boolean flags specifying model
       being estimated and template behaviour
    */
    DATA_STRUCT(flags, egf::flags_t);
    flags.do_simulate = this->do_simulate;
    
    
    /* Fixed effects infrastructure ========================================= */

    /* Combined fixed effects model matrix
       - dim=(N, length(beta)) if do_sparse_X=false, dim=(N, 0) otherwise
    */
    DATA_MATRIX(Xd);
    /* 
       - dim=(N, length(beta)) if do_sparse_X=true,  dim=(N, 0) otherwise
    */
    DATA_SPARSE_MATRIX(Xs);

    /* "Factor" splitting fixed effects model matrix columns by relation 
       to a nonlinear model parameter
       - length=length(beta)
       - val={0,...,p-1}
    */
    DATA_IVECTOR(beta_index);

    /* Number of fixed effects model matrix columns by nonlinear model
       parameter
       - in R: beta_index_tab = c(table(beta_index))
       - length=p
    */
    DATA_IVECTOR(beta_index_tab);

    /* Number of nonlinear model parameters */
    int p = beta_index_tab.size();

    /* Fixed effects coefficient matrix, to be multiplied on the left by `X` */
    Eigen::SparseMatrix<Type> beta_matrix(beta.size(), p);
    beta_matrix.reserve(beta_index_tab);
    for (int i = 0; i < beta.size(); ++i)
    {
        beta_matrix.insert(i, beta_index(i)) = beta(i);
    }

    
    /* Random effects infrastructure ======================================== */

    /* Combined random effects model matrix
       - dim=(N, length(b))
    */
    DATA_SPARSE_MATRIX(Z);
    flags.do_random_effects = (Z.cols() > 0);    

    /* Random effects coefficient matrix, to be multiplied on the left by `Z` */
    Eigen::SparseMatrix<Type> b_matrix;

    /* Number of random effects blocks */
    int M;
    
    /* List of random effects blocks
       - each block is a matrix whose columns vectors are contiguous segments
         of parameter vector `b`
       - blocks, row vectors, and column vectors correspond to random effect
         terms, nonlinear model parameters, and group levels, respectively
       - column vectors represent realizations of a multivariate normal random 
         variable with zero mean, unit variance, and unstructured covariance
    */
    vector< matrix<Type> > list_of_blocks;

    /* List of corresponding vectors of standard deviations */
    vector< vector<Type> > list_of_sd;

    /* List of corresponding vectors of other covariance parameters */
    vector< vector<Type> > list_of_chol;

    /* List of corresponding negative log density functions */
    vector< density::MVNORM_t<Type> > list_of_nld;
    
    if (flags.do_random_effects)
    {
        /* Data ------------------------------------------------------------- */

        /* "Factor" splitting random effects model matrix
	   columns by relation to a nonlinear model parameter
	   - length=length(b)
	   - val={0,...,p-1}
	*/
	DATA_IVECTOR(b_index);

	/* Number of random effects model matrix columns 
	   by nonlinear model parameter
	   - in R: b_index_tab = c(table(b_index))
	   - length=p
	*/
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
	M = block_cols.size();

	
	/* Initialization --------------------------------------------------- */

	b_matrix.resize(b.size(), p);
        b_matrix.reserve(b_index_tab);
	list_of_blocks.resize(M);
	list_of_sd.resize(M);
	list_of_chol.resize(M);
	list_of_nld.resize(M);
	
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

    /* Offset component of response matrix
       - dim=(N, p)
    */
    DATA_MATRIX(Y);

    /* Structure storing integer column indices of nonlinear model parameters 
       in response matrix
       - val={0,...p-1} if parameter is used, val=-1 otherwise
    */
    DATA_STRUCT(indices, egf::indices_t);
    
    /* Add fixed effects component */
    Y += (flags.do_sparse_X ? Xs : Xd) * beta_matrix;

    /* Add random effects component */
    if (flags.do_random_effects)
    {
        Y += Z * b_matrix;
    }
    
    REPORT(Y);
    if (!(flags.do_simulate || flags.do_predict))
    {
        ADREPORT(Y);
    }

    
    /* Time series infrastructure =========================================== */

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
       - length=N
       - val={0,...,6} if do_day_of_week=true, val=-1 otherwise
    */
    DATA_IVECTOR(day1);
    

    /* Negative log likelihood ============================================== */

    Type res = Type(0.0);

    
    /* Observation likelihood ----------------------------------------------- */

    if (flags.do_trace && !flags.do_simulate)
    {
        std::cout << "Y = \n" << Y << "\n";
	printf("nll initialized to %.6e\n", asDouble(res));
	std::cout << "commencing loop over observations\n";
    }
    res += nll_obs(this, time, time_seg_len, x, Y, indices, flags, day1);
    if (flags.do_simulate)
    {
        REPORT(x);
	return Type(0.0);
    }
    if (flags.do_trace)
    {
        std::cout << "loop over observations complete\n";
        printf("nll is %.6e\n", asDouble(res));
    }

    
    /* Random effect likelihood --------------------------------------------- */

    if (flags.do_random_effects) {
        if (flags.do_trace)
	{
	    std::cout << "commencing loop over random effects\n";
	}
	res += nll_ran(this, list_of_blocks, list_of_nld, flags);
	if (flags.do_trace)
	{
	    std::cout << "loop over random effects complete\n";
	    printf("nll is %.6e\n", asDouble(res));
	}
    }

    
    /* Top level parameter value likelihood --------------------------------- */

    if (flags.do_regularize_top)
    {
        /* List of vectors of hyperparameters for regularization
	   - length=p
	*/
        DATA_STRUCT(hyperparameters_top, egf::list_of_vectors_t);

	if (flags.do_trace)
	{
	    std::cout << "commencing loop over regularized top level parameters\n";
	}
        res += nll_top(this, Y, hyperparameters_top, flags);
        if (flags.do_trace)
	{
	    std::cout << "loop over regularized top level parameters complete\n";
	    printf("nll is %.6e\n", asDouble(res));
	}
    }


    /* Bottom level parameter value likelihood ------------------------------ */

    if (flags.do_regularize_bottom)
    {
        /* List of vectors of hyperparameters for regularization
	   - length=length(beta)+length(theta)
	*/
        DATA_STRUCT(hyperparameters_bottom, egf::list_of_vectors_t);

	if (flags.do_trace)
	{
	    std::cout << "commencing loop over regularized bottom level parameters\n";
	}
        res += nll_bot(this, beta, theta, hyperparameters_bottom, flags);
        if (flags.do_trace)
	{
	    std::cout << "loop over regularized bottom level parameters complete\n";
	    printf("nll is %.6e\n", asDouble(res));
	}
    }

    
    /* Prediction =========================================================== */

    if (flags.do_predict)
    {
        /* Flags ------------------------------------------------------------ */
        
        /* What should be predicted? */
	DATA_IVECTOR(flag_what);

	/* Let c(t) be cumulative incidence since time -Inf. Then:
	   0. Log interval incidence
	      
	      = log(diff(c(t)))
	*/
	bool do_predict_log_interval = (flag_what(0) == 1);
	/* 1. Log cumulative incidence since time -Inf
	      
	      = log(c(t))
	*/
        bool do_predict_log_cumulative = (flag_what(1) == 1);
	/* 2. Log per capita growth rate
	      
	      = log(c'(t) / c(t))
	*/
	bool do_predict_log_rt = (flag_what(2) == 1);

      
        /* Data ------------------------------------------------------------- */
      
        /* Time series segments for which predictions should be computed 
	   - length=N'
	   - val={0,...,N-1}
	*/
        DATA_IVECTOR(subset);

	/* Concatenated vectors of time points 
	   - length=n'
	*/
	DATA_VECTOR(new_time);

	/* Segment lengths
	   - length=N'
	*/
	DATA_IVECTOR(new_time_seg_len);

	/* Number of segments */
	N = new_time_seg_len.size();

	/* Segment initial days of week
	   - length=N'
	   - val={0,...,6} if do_day_of_week=true, val=-1 otherwise
	*/
	DATA_IVECTOR(new_day1);


	/* Log cumulative incidence since time -Inf, no day of week effects - */

	int n;
	vector<Type> Y_row;
	vector< vector<Type> > list_of_predict(N);
	for (int s = 0, i = 0; s < N; ++s)
	{
	    n = new_time_seg_len(s);
	    /* t */
	    list_of_predict(s) = new_time.segment(i, n);
	    /* t <- log(c(t)) */
	    egf::eval_log_curve(list_of_predict(s),
				(vector<Type>) Y.row(subset(s)),
				indices,
				flags.flag_curve);
	    i += n;
	}


	/* Prediction variables --------------------------------------------- */

	if (do_predict_log_rt)
	{
	    vector<Type> log_rt(new_time.size());
	    vector<Type> tmp;
	    for (int s = 0, i = 0; s < N; ++s)
	    {
		n = new_time_seg_len(s);
		Y_row = Y.row(subset(s));
		/* log(c(t)) <- log(c'(t) / c(t)) */
		if (do_predict_log_interval || do_predict_log_cumulative)
		{
		    /* Modify a copy of log(c(t)) in place, 
		       because we still need log(c(t)) 
		    */ 
		    tmp = list_of_predict(s);
		    egf::eval_log_rt_exact(tmp,
					   Y_row,
					   indices,
					   flags.flag_curve);
		    log_rt.segment(i, n) = tmp;
		}
		else
		{
		    /* Modify log(c(t)) in place,
		       because we no longer need log(c(t)) 
		    */
		    egf::eval_log_rt_exact(list_of_predict(s),
					   Y_row,
					   indices,
					   flags.flag_curve);
		    log_rt.segment(i, n) = list_of_predict(s);
		}
		i += n;
	    }
	    REPORT(log_rt);
	    ADREPORT(log_rt);
	}

	if (do_predict_log_interval || do_predict_log_cumulative)
	{
	    vector<Type> log_interval;
	    vector<Type> log_cumulative;
	    if (do_predict_log_interval)
	    {
	        log_interval.resize(new_time.size() - N);
	    }
	    if (do_predict_log_cumulative)
	    {
	        log_cumulative.resize(new_time.size());
	    }
	    for (int s = 0, i = 0; s < N; ++s)
	    {
	        n = new_time_seg_len(s);
		Y_row = Y.row(subset(s));
		if (flags.do_excess)
		{
		    /* log(c(t)) <- log(b * t + c(t)) */
		    egf::logspace_add_baseline(list_of_predict(s),
					       (vector<Type>) new_time.segment(i, n),
					       Y_row(indices.index_log_b));
		}
		if (do_predict_log_cumulative)
		{
		    if (flags.do_day_of_week)
		    {
		        log_cumulative(i) = list_of_predict(s)(0);
		    }
		    else
		    {
			log_cumulative.segment(i, n) = list_of_predict(s);
		    }
		}
		if (n > 1 && (do_predict_log_interval || (do_predict_log_cumulative && flags.do_day_of_week)))
		{  
		    /* log(c(t)) <- log(diff(c(t))) */
		    egf::logspace_diff(list_of_predict(s));
		    if (flags.do_day_of_week)
		    {
			/* log(diff(c(t))) <- log(diff(c(t)) * w(t[-n], t[-1])) */
			egf::logspace_add_offsets(list_of_predict(s),
						  Y_row(indices.index_log_w1),
						  Y_row(indices.index_log_w2),
						  Y_row(indices.index_log_w3),
						  Y_row(indices.index_log_w4),
						  Y_row(indices.index_log_w5),
						  Y_row(indices.index_log_w6),
						  new_day1(s));
		    }
		    if (do_predict_log_interval)
		    {
			log_interval.segment(i - s, n - 1) = list_of_predict(s);
		    }
		    if (do_predict_log_cumulative)
		    {
			for (int k = 0; k < n - 1; ++k)
			{
			    log_cumulative(i + k + 1) = logspace_add(log_cumulative(i + k), list_of_predict(s)(k));
			}
		    }
		}
		i += n;
	    }
	    if (do_predict_log_interval)
	    {
	        REPORT(log_interval);
		ADREPORT(log_interval);
	    }
	    if (do_predict_log_cumulative)
	    {
	        REPORT(log_cumulative);
		ADREPORT(log_cumulative);
	    }
	}
    }

    return res;
}
