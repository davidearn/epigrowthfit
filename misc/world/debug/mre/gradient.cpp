#include <TMB.hpp>

namespace mre
{

/* Similar to R function `is.na` */
template<class Type>
bool is_na(Type x)
{
    return R_IsNA(asDouble(x));
}

/* Evaluate log logistic function at `time` (in place) */
template<class Type>
void eval_log_logistic(vector<Type> &time,
		       Type log_r,
		       Type log_tinfl,
		       Type log_K)
{
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);
    int n = time.size();
    for (int i = 0; i < n; ++i)
    {
        time(i) = log_K - logspace_add(Type(0.0), -r * (time(i) - tinfl));
    }
}

/* Compute `log(diff(x))` given increasing `log(x)` (in place) */
template<class Type>
void logspace_diff(vector<Type> &log_x)
{
    int n = log_x.size() - 1;
    for (int i = 0; i < n; ++i)
    {
        log_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    log_x.conservativeResize(n);
}

} // namespace egf

template<class Type>
Type objective_function<Type>::operator() ()
{
#ifdef _OPENMP
    this->max_parallel_regions = omp_get_max_threads();
#endif
  
    /* Data ================================================================= */

    /* Concatenated vectors of time points */
    DATA_VECTOR(time); // length=n

    /* Segment lengths */
    DATA_IVECTOR(time_seg_len); // length=N

    /* Number of segments */
    int N = time_seg_len.size();

    /* Concatenated vectors of observations */
    DATA_VECTOR(x); // length=n-N

    /* Number of top level (logistic and dispersion) parameters */
    int p = 4;

    /* Random effects model matrix */
    DATA_SPARSE_MATRIX(Z); // dim=(N, p*N)

    
    /* Parameters =========================================================== */

    /* Concatenated fixed effects coefficient vectors */
    PARAMETER_VECTOR(beta); // length=p

    /* Concatenated random effects coefficient vectors (unit variance scale) */
    PARAMETER_VECTOR(b); // length=p*N

    /* Concatenated random effects covariance parameters */
    PARAMETER_VECTOR(theta); // length=p*(p+1)/2

    
    /* Random effects infrastructure ======================================== */

    /* Random effects block */
    matrix<Type> block = b.matrix();
    block.resize(p, N);
    
    /* Vector of standard deviations */
    vector<Type> sd = exp(theta.head(p));
    
    /* Vector of other covariance parameters from Cholesky factor */
    vector<Type> chol = theta.tail(p * (p - 1) / 2);
    
    /* Random effects coefficient matrix */
    Eigen::SparseMatrix<Type> b_matrix(p * N, p);
    vector<int> nnz(p);
    nnz.fill(N);
    b_matrix.reserve(nnz);
    for (int i = 0, j = 0; i < p * N; ++i, j = i % p)
    {
        b_matrix.insert(i, j) = sd(j) * b(i);
    }
    
    
    /* Response matrix ====================================================== */

    matrix<Type> Y(N, p);
    for (int j = 0; j < p; ++j)
    {
        Y.col(j).fill(beta(j));
    }
    Y += Z * b_matrix;


    /* Negative log likelihood ============================================== */

    Type nll = Type(0.0);
    
    /* Observation likelihood ----------------------------------------------- */

    vector<Type> tmp;
    Type log_mu;
    Type log_size;
    Type log_var_minus_mu;
    int n;
    for (int s = 0, i = 0; s < N; ++s) // loop over segments
    {
        n = time_seg_len(s);

	/* t */
	tmp = time.segment(i, n);
	/* t <- log(c(t)) ... where `c` is the logistic function */
	mre::eval_log_logistic(tmp, Y(s, 0), Y(s, 1), Y(s, 2));
	/* log(c(t)) <- log(diff(c(t))) */
	mre::logspace_diff(tmp);
	
	log_size = Y(s, 3);
	for (int k = 0; k < n - 1; ++k) // loop over observations
	{
	    if (this->parallel_region() && !mre::is_na(x(i-s+k)))
	    {
	        log_mu = tmp(k);
		log_var_minus_mu = Type(2.0) * log_mu - log_size;
		nll -= dnbinom_robust(x(i-s+k), log_mu, log_var_minus_mu, true);
	    }
	}
	i += n;
    }
    
    /* Random effect likelihood --------------------------------------------- */

    density::UNSTRUCTURED_CORR_t<Type> nld(chol);
    for (int j = 0; j < N; ++j)
    {
        if (this->parallel_region())
	{
	    nll += nld(block.col(j));
	}
    }

    return nll;
}
