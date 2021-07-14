#include <TMB.hpp>

template<class Type>
vector<Type> eval_log_logistic(vector<Type> time,
			       Type log_r,
			       Type log_tinfl,
			       Type log_K)
{
    int n = time.size();
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);
    vector<Type> res(n);
    for (int i = 0; i < n; i++)
    {
        res(i) = log_K - logspace_add(Type(0), -r * (time(i) - tinfl));
    }
    return res;
}

template<class Type>
vector<Type> eval_log_logistic(vector<Type> time,
			       vector<int> len,
			       vector<Type> log_r,
			       vector<Type> log_tinfl,
			       vector<Type> log_K)
{
    int n;
    vector<Type> res(time.size());
    for (int s = 0, i = 0; s < len.size(); s++)
    {
        n = len(s);
        vector<Type> time0 = time.segment(i, n);
        res.segment(i, n) = 
	  eval_log_logistic(time0, log_r(s), log_tinfl(s), log_K(s));
	i += n;
    }
    return res;
}

template<class Type>
vector<Type> logspace_diff(vector<Type> log_x)
{
    int n = log_x.size();
    vector<Type> log_diff_x(n - 1);
    for (int i = 0; i < n - 1; i++)
    {
        log_diff_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    return log_diff_x;
}

template<class Type>
vector<Type> logspace_diff(vector<Type> log_x, vector<int> len)
{
    int n;
    vector<Type> log_diff_x(log_x.size() - len.size());
    for (int s = 0, i = 0; s < len.size(); s++)
    {
        n = len(s);
        vector<Type> log_x_segment = log_x.segment(i + s, n);
        log_diff_x.segment(i, n - 1) = logspace_diff(log_x_segment);
	i += n - 1;
    }
    return log_diff_x;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    /* Parameters =========================================================== */
  
    /* Concatenated fixed effects coefficient vectors */
    PARAMETER_VECTOR(beta); // length=p

    /* Concatenated random effects coefficient vectors (unit variance scale) */
    PARAMETER_VECTOR(b); // length=p*N

    /* Concatenated random effects covariance parameters */
    PARAMETER_VECTOR(theta); // length=p*(p+1)/2

    
    /* Data ================================================================= */

    /* Concatenated time series segments */
    DATA_VECTOR(time); // length=n
    DATA_VECTOR(x); // length=n-N

    /* Segment lengths */
    DATA_IVECTOR(time_seg_len); // length=N
    
    /* Number of segments */
    int N = time_seg_len.size();

    /* Random effects model matrix */
    DATA_SPARSE_MATRIX(Z); // dim=(N, p*N)

    /* Number of nonlinear and dispersion model parameters */
    int p = 4;

    /* Column indices in response matrix `Y` */
    DATA_INTEGER(j_log_r);
    DATA_INTEGER(j_log_tinfl);
    DATA_INTEGER(j_log_K);
    DATA_INTEGER(j_log_nbdisp);

    /* Trace (1=yes, 0=no) */
    DATA_INTEGER(trace);
    bool do_trace = (trace > 0);

    
    /* Random effects infrastructure ======================================== */
    
    /* Vector of standard deviations */
    vector<Type> sd = exp(theta.segment(0, p));
    
    /* Vector of other covariance parameters from Cholesky factor */
    vector<Type> chol = theta.segment(p, p * (p - 1) / 2);

    /* Function evaluating negative log density */
    density::UNSTRUCTURED_CORR_t<Type> nld(chol);

    /* Random effects block */
    matrix<Type> block(p, N);
    for (int i = 0, j = 0; j < N; i += p, j++)
    {
        block.col(j) = (vector<Type>) b.segment(i, p);
    }

    /* Random effects coefficient matrix */
    Eigen::SparseMatrix<Type> b_scaled_matrix(p * N, p);
    vector<int> nnz(p);
    nnz.fill(N);
    b_scaled_matrix.reserve(nnz);
    for (int i = 0, j = 0; i < p * N; i++, j = i % 4)
    {
        b_scaled_matrix.insert(i, j) = sd(j) * b(i);
    }
    
    
    /* Response matrix ====================================================== */

    matrix<Type> Y(N, p);
    for (int j = 0; j < p; j++)
    {
        Y.col(j).fill(beta(j));
    }
    Y += Z * b_scaled_matrix;

    
    /* Predictions ========================================================== */

    /* Log cumulative incidence */
    vector<Type> log_logistic =
        eval_log_logistic(time,
			  time_seg_len,
			  (vector<Type>) Y.col(j_log_r),
			  (vector<Type>) Y.col(j_log_tinfl),
			  (vector<Type>) Y.col(j_log_K));
      
    /* Log interval incidence */
    vector<Type> log_cases = logspace_diff(log_logistic, time_seg_len);


    /* Negative log likelihood ============================================== */

    parallel_accumulator<Type> nll(this);
    Type nll_term;

    if (do_trace)
    {
        printf("nll initialized to %.6e\n", asDouble(nll));
	std::cout << "commencing loop over observations\n";
    }

    Type log_var_minus_mu;
    int n;
    for (int s = 0, i = 0; s < N; s++)
    {
        n = time_seg_len(s) - 1;
	for (int k = 0; k < n; k++)
	{
	    log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
	    nll_term = -dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
	    nll += nll_term;

	    if (do_trace)
	    {
	        printf("at index %d of segment %d: nll term is %.6e\n",
		       k, s, asDouble(nll_term));
		printf("  \\__ -log(dnbinom(x = %d, mu = %.6e, size = %.6e))\n",
		       (int) asDouble(x(i+k)),
		       asDouble(exp(log_cases(i+k))),
		       asDouble(exp(Y(s, j_log_nbdisp))));
	    }
	}
	i += n;
    }

    if (do_trace)
    {
        std::cout << "loop over observations complete\n";
        printf("nll is %.6e\n", asDouble(nll));
        std::cout << "commencing loop over random effects\n";
    }

    for (int j = 0; j < N; j++)
    {
	nll_term = nld(block.col(j));
	nll += nll_term;

	if (do_trace)
	{
	    printf("at column %d: nll term is %.6e\n", j, asDouble(nll_term));
	}
    }

    if (do_trace)
    {
        std::cout << "loop over random effects complete\n";
        printf("nll is %.6e\n", asDouble(nll));
    }

    return nll;
}
