#include <TMB.hpp>

template<class Type>
bool is_NA_real_(Type x)
{
    return R_IsNA(asDouble(x));
}

template<class Type>
bool is_finite(Type x)
{
    return R_finite(asDouble(x));
}

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
    
    /* Combined random effects model matrix */
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
    REPORT(sd);

    /* Other covariance parameters from Cholesky factor */
    vector<Type> chol = theta.segment(p, p * (p - 1) / 2);
    REPORT(chol);

    /* Correlation matrix */
    matrix<Type> cor = density::UNSTRUCTURED_CORR(chol).cov();
    REPORT(cor);

    /* Random effects block */
    matrix<Type> block(p, N);
    for (int i = 0, j = 0; j < N; i += p, j++)
    {
        block.col(j) = (vector<Type>) b.segment(i, p);
    }
    REPORT(block);

    /* Random effects coefficient matrix */
    Eigen::SparseMatrix<Type> b_scaled_matrix(p * N, p);
    vector<int> nnz(p);
    nnz.fill(N);
    b_scaled_matrix.reserve(nnz);
    for (int i = 0, j = 0; i < p * N; i++, j = i % 4)
    {
        b_scaled_matrix.insert(i, j) = sd(j) * b(i);
    }
    REPORT(b_scaled_matrix);
    
    
    /* Compute (N, p) response matrix ======================================= */

    matrix<Type> Y(N, p);
    for (int j = 0; j < p; j++)
    {
        Y.col(j).fill(beta(j));
    }
    Y += Z * b_scaled_matrix;
    REPORT(Y);
    
    
    /* Compute predictions ================================================== */

    /* Log cumulative incidence */
    vector<Type> log_logistic = eval_log_logistic_n(t,
						    t_seg_len,
						    (vector<Type>) Y.col(j_log_r),
						    (vector<Type>) Y.col(j_log_tinfl),
						    (vector<Type>) Y.col(j_log_K));
    REPORT(log_logistic);
      
    /* Log interval incidence */
    vector<Type> log_cases = logspace_diff_n(log_logistic, t_seg_len);
    REPORT(log_cases);


    /* Compute negative log likelihood ====================================== */

    parallel_accumulator<Type> nll(this);
    Type nll_term;
    Type log_var_minus_mu;

    // std::cout << "Y = \n" << Y << "\n";
    // printf("nll initialized to 0\ncommencing loop over observations\n");

    for (int i = 0, s = 0; s < N; i += t_seg_len(s) - 1, s++) /* loop over segments */
    {
        for (int k = 0; k < t_seg_len(s) - 1; k++) /* loop over within-segment index */
	{
	    if (!is_NA_real_(x(i+k)))
	    {
	        log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
		nll_term = -dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
		nll += nll_term;

		// if (!is_finite(nll_term))
		// {
		//     printf("at index %d of segment %d: nll term is non-finite\n", k, s);
		// }
		// else if (asDouble(nll_term) > 1.0e+09)
		// {
		//     printf("at index %d of segment %d: nll term exceeds 1.0e+09\n", k, s);
		// }
		
		// printf("at index %d of segment %d: nll %.6e x %d mu %.6e size %.6e\n", k, s,
		//        asDouble(nll_term),
		//        (int) asDouble(x(i+k)),
		//        asDouble(exp(log_cases(i+k))),
		//        asDouble(exp(Y(s, j_log_nbdisp))));
	    }
	}

	// std::cout << "Y.row(" << s << ") = " << Y.row(s) << "\n";
    }

    printf("loop over observations complete\nnll is %5.6e\n", asDouble(nll));
    printf("commencing loop over random effects\n");

    density::MVNORM_t<Type> N_0_Sigma(cor); /* function returning negative log density */
    for (int j = 0; j < N; j++)
    {
	nll_term = N_0_Sigma(block.col(j));
	nll += nll_term;

	// if (!is_finite(nll_term))
	// {
	//     printf("at column %d: nll term is non-finite\n", j);
	// }
	// else if (asDouble(nll_term) > 1.0e+09)
	// {
	//     printf("at column %d: nll term exceeds 1.0e+09\n", j);
	// }

	// printf("at column %d: nll %.6e\n", j, asDouble(nll_term));
    }

    // printf("loop over random effects complete\nnll is %5.6e\n", asDouble(nll));

    return nll;
}
