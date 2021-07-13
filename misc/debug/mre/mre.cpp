#include <TMB.hpp>

template<class Type>
vector<Type> eval_log_logistic(vector<Type> t,
			       Type log_r,
			       Type log_tinfl,
			       Type log_K)
{
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);
    vector<Type> res(t.size());
    for (int i = 0; i < t.size(); i++)
    {
        res(i) = log_K - logspace_add(Type(0), -r * (t(i) - tinfl));
    }
    return res;
}

template<class Type>
vector<Type> eval_log_logistic(vector<Type> t,
			       vector<int> len,
			       vector<Type> log_r,
			       vector<Type> log_tinfl,
			       vector<Type> log_K)
{
    vector<Type> res(t.size());
    for (int i = 0, s = 0; s < len.size(); i += len(s), s++)
    {
        vector<Type> t0 = t.segment(i, len(s));
        res.segment(i, len(s)) =
	  eval_log_logistic(t0, log_r(s), log_tinfl(s), log_K(s));
    }
    return res;
}

template<class Type>
vector<Type> logspace_diff(vector<Type> log_x)
{
    vector<Type> log_diff_x(log_x.size() - 1);
    for (int i = 0; i < log_x.size() - 1; i++)
    {
        log_diff_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    return log_diff_x;
}

template<class Type>
vector<Type> logspace_diff(vector<Type> log_x, vector<int> len)
{
    vector<Type> log_diff_x(log_x.size() - len.size());
    for (int s = 0, i = 0; s < len.size(); i += len(s) - 1, s++)
    {
        vector<Type> log_x_segment = log_x.segment(i + s, len(s));
        log_diff_x.segment(i, len(s) - 1) = logspace_diff(log_x_segment);
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
    DATA_VECTOR(t); // length=n
    DATA_VECTOR(x); // length=n-N

    /* Segment lengths */
    DATA_IVECTOR(t_seg_len); // length=N
    
    /* Random effects model matrix */
    DATA_SPARSE_MATRIX(Z); // dim=(N, p*N)

    /* Indices */
    DATA_INTEGER(j_log_r);
    DATA_INTEGER(j_log_tinfl);
    DATA_INTEGER(j_log_K);
    DATA_INTEGER(j_log_nbdisp);


    /* Set up =============================================================== */
    
    /* Number of segments */
    int N = t_seg_len.size();

    /* Number of nonlinear and dispersion model parameters */
    int p = 4;

    /* Standard deviations */
    vector<Type> sd = exp(theta.segment(0, p));
    
    /* Other covariance parameters from Cholesky factor */
    vector<Type> chol = theta.segment(p, p * (p - 1) / 2);

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
    
    
    /* Compute (N, p) response matrix ======================================= */

    matrix<Type> Y(N, p);
    for (int j = 0; j < p; j++)
    {
        Y.col(j).fill(beta(j));
    }
    Y += Z * b_scaled_matrix;    
    
    /* Compute predictions ================================================== */

    /* Log cumulative incidence */
    vector<Type> log_logistic = eval_log_logistic(t,
						  t_seg_len,
						  (vector<Type>) Y.col(j_log_r),
						  (vector<Type>) Y.col(j_log_tinfl),
						  (vector<Type>) Y.col(j_log_K));
      
    /* Log interval incidence */
    vector<Type> log_cases = logspace_diff(log_logistic, t_seg_len);


    /* Compute negative log likelihood ====================================== */

    parallel_accumulator<Type> nll(this);
    Type nll_term;
    Type log_var_minus_mu;

    // printf("nll initialized to 0\ncommencing loop over observations\n");

    for (int i = 0, s = 0; s < N; i += t_seg_len(s) - 1, s++)
    {
        for (int k = 0; k < t_seg_len(s) - 1; k++)
	{
	    log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
	    nll_term = -dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
	    nll += nll_term;
	    // printf("at index %d of segment %d: nll %.6e x %d mu %.6e size %.6e\n", k, s,
	    //        asDouble(nll_term),
	    //        (int) asDouble(x(i+k)),
	    //        asDouble(exp(log_cases(i+k))),
	    //        asDouble(exp(Y(s, j_log_nbdisp))));
	}
    }

    // printf("loop over observations complete\nnll is %.6e\n", asDouble(nll));
    // printf("commencing loop over random effects\n");

    density::UNSTRUCTURED_CORR_t<Type> N_0_Sigma(chol); /* function returning negative log density */
    for (int j = 0; j < N; j++)
    {
	nll_term = N_0_Sigma(block.col(j));
	nll += nll_term;
	// printf("at column %d: nll %.6e\n", j, asDouble(nll_term));
    }

    // printf("loop over random effects complete\nnll is %.6e\n", asDouble(nll));

    return nll;
}
