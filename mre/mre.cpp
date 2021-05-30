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
Type eval_log_logistic(Type t, Type log_r, Type log_tinfl, Type log_K)
{
    // log(c(t))
    // = log(K / (1 + exp(-r * (t - tinfl))))
    // = log(K) - log(1 + exp(-r * (t - tinfl)))
    return log_K - logspace_add(Type(0), -exp(log_r) * (t - exp(log_tinfl)));
}

template<class Type>
vector<Type> logspace_diff_1(vector<Type> log_x)
{
    vector<Type> log_diff_x(log_x.size() - 1);
    for (int i = 0; i < log_x.size() - 1; i++)
    {
        log_diff_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    return log_diff_x;
}

template<class Type>
vector<Type> logspace_diff_n(vector<Type> log_x, vector<int> len)
{
    vector<Type> log_diff_x(log_x.size() - len.size());
    for (int s = 0, i = 0; s < len.size(); s++) // loop over segments
    {
        vector<Type> log_x_segment = log_x.segment(i + s, len(s));
        log_diff_x.segment(i, len(s) - 1) = logspace_diff_1(log_x_segment); 
        i += len(s) - 1; // increment reference index
    }
    return log_diff_x;
}

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
    DATA_INTEGER(j_log_tinfl);
    DATA_INTEGER(j_log_K);
    DATA_INTEGER(j_log_nbdisp);
    
    // Combined fixed effects model matrix
    // * dim=(N, sum(beta_seg_len))
    DATA_MATRIX(X);

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

    // Flag for trace
    DATA_INTEGER(trace_flag);
    bool trace = (trace_flag == 1);
    

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
    Y += X * beta_as_matrix;

    // 3. Add random effects component:

    if (any_RE)
    {
        // Arrange elements of `b_scaled` in (sum(b_seg_len), p) matrix
        Eigen::SparseMatrix<Type> b_scaled_as_matrix(b.size(), p);
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

    vector<Type> log_curve(t.size());
    for (int s = 0, i = 0; s < t_seg_len.size(); s++) // loop over segments
    {
        for (int k = 0; k < t_seg_len(s); k++) // loop over within-segment index
	{
	    log_curve(i+k) = eval_log_logistic(t(i+k), Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K));
	}
	i += t_seg_len(s);
    }
    vector<Type> log_cases = logspace_diff_n(log_curve, t_seg_len);


    // Compute likelihood ======================================================

    // Negative log likelihood
    Type nll = Type(0);
    Type nll_term;
    Type log_var_minus_mu;
    bool print_Y_row;

    if (trace)
    {
        Rprintf("nll initialized to 0\ncommencing loop over observations\n");
    }

    for (int s = 0, i = 0; s < N; s++) // loop over segments
    {
        print_Y_row = false;
	
        for (int k = 0; k < t_seg_len(s) - 1; k++) // loop over within-segment index
	{
	    if (!is_NA_real_(x(i+k)))
	    {
	        log_var_minus_mu = Type(2) * log_cases(i+k) - Y(s, j_log_nbdisp);
		nll_term = -dnbinom_robust(x(i+k), log_cases(i+k), log_var_minus_mu, true);
		// usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)

		if (trace)
		{
		    if (!is_finite(nll_term))
		    {
			Rprintf("at index %d of segment %d: nll term is non-finite\n", k, s);
			print_Y_row = true;
		    }
		    else if (asDouble(nll_term) > 1.0e09)
		    {
			Rprintf("at index %d of segment %d: nll term exceeds 1.0e09\n", k, s);
			print_Y_row = true;
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
        if (trace)
	{
	    Rprintf("commencing loop over random effects\n");
	}
      
        for (int m = 0; m < M; m++) // loop over blocks
    	{
	    for (int j = 0; j < block_cols(m); j++) // loop over block columns
    	    {
    	        nll_term = density::MVNORM(cor_list(m))(block_list(m).col(j));

		if (trace)
		{
		    if (!is_finite(nll_term))
		    {
			Rprintf("at column %d of block %d: nll term is non-finite\n", j, m);
		    }
		    else if (asDouble(nll_term) > 1.0e12)
		    {
			Rprintf("at column %d of block %d: nll term exceeds 10^12\n", j, m);
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
    
    return nll;
}
