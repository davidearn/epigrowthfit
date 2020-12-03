#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>

template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA(Type x)
{
    return R_IsNA(asDouble(x));
}



// https://kaskr.github.io/adcomp/Introduction.html
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Flags
    DATA_INTEGER(predict_flag); // predict   (1=predict, 0=fit)
    DATA_INTEGER(spX_flag);     // sparse X  (1=sparse,  0=dense)
    DATA_INTEGER(spZ_flag);     // sparse Z  (1=sparse,  0=dense)

    // Data
    // time series
    DATA_VECTOR(t);   // length=N
    DATA_VECTOR(x);   // length=N
    // window
    DATA_IVECTOR(w);   // length=N,  {0,...,nlevels(w)-1}
    // nlevels
    DATA_IVECTOR(fnl);  // length=nFE=np
    DATA_IVECTOR(rnl);  // length=nRE
    // RE indicators
    DATA_IMATRIX(rid);  // nrow=nRE,  ncol=np,  {0,1}
    // parameter indices
    DATA_INTEGER(i_log_r);
    DATA_INTEGER(i_log_c0);     // !NA  iff  curve="exponential"
    DATA_INTEGER(i_log_thalf);  //  NA  iff  curve="exponential"
    DATA_INTEGER(i_log_K);      //  NA  iff  curve="exponential"
    DATA_INTEGER(i_log_p);      // !NA  iff  curve="richards"
    DATA_INTEGER(i_log_nbdisp); // !NA  iff  distr="nbinom"
    DATA_INTEGER(i_log_b);      //  NA  iff  include_baseline=FALSE
    
    // Parameters
    // fixed effects
    PARAMETER_VECTOR(beta); // length=sum(nlevels(FE i))
    // random effects
    PARAMETER_VECTOR(b);    // length=sum(nlevels(RE i)*np(RE i)) 
    // sd random effects
    PARAMETER_VECTOR(sd_b); // length=sum(np(RE i))

    // Number of parameters
    int np = rid.cols();

    // Any random effects?
    bool anyRE = (rid.rows() > 0);

    
    // Compute parameter values at every time point ============================
    // NOTE: parameters are constant within but not across fitting windows

    // (N x np) matrix of parameter values
    matrix<Type> Q(t.size(), np);
    
    // Fixed effects component:
    // matrix multiply
    // X * beta2
    // = (N x sum(nlevels(FE i))) * (sum(nlevels(FE i)) x np)
    // = (N x np)
    // need to construct block matrix `beta2` from vector `beta`

    Eigen::SparseMatrix<Type> beta2_sparse(fnl.sum(), np);
    matrix<Type> beta2_dense(fnl.sum(), np);
    
    if (spX_flag == 1)
    {
        beta2_sparse.reserve(fnl);
	for (int j = 0, i = 0; j < np; j++) // loop over parameters
	{
	    for (int l = 0; l < fnl(j); l++) // loop over levels
	    {
		beta2_sparse.insert(i + l, j) = beta(i + l);
	    }
	    i += fnl(j);
	}
        beta2_sparse.makeCompressed();
	DATA_SPARSE_MATRIX(X);
	Q = X * beta2_sparse;
    }
    else
    {
        beta2_dense.fill(Type(0));
	for (int j = 0, i = 0; j < np; j++) // loop over parameters
	{
	    beta2_dense.block(i, j, fnl(j), 1) = beta.segment(i, fnl(j));
	    i += fnl(j);
	}
	DATA_MATRIX(X);
	Q = X * beta2_dense;
    }
    
    // Random effects component:
    // matrix multiply
    // Z * b2
    // = (N x sum(nlevels(RE i))) * (sum(nlevels(RE i)) x np)
    // = (N x np)
    // need to construct block matrix `b2` from vector `b`
    // will require some setting up...

    Eigen::SparseMatrix<Type> b2_sparse(rnl.sum(), np);
    matrix<Type> b2_dense(rnl.sum(), np);

    // A list of (nlevels(RE i) x np(RE i)) matrices gathering related elements
    // of vector `b`: rows of each matrix follow a zero-mean multivariate normal
    // distribution with a common covariance matrix
    vector< matrix<Type> > re_list(rid.rows());

    // A list of covariance matrices corresponding elementwise to `re_list`
    vector< matrix<Type> > cov_list(rid.rows());
	
    // A list of s.d. vectors corresponding elementwise to `re_list`
    vector< vector<Type> > sd_list(rid.rows());
    
    if (anyRE)
    {

        // Form random effect matrices -----------------------------------------

	// Number of parameters per random effect
	vector<int> rnp = rid.rowwise().sum();

	for (int r = 0, k = 0; r < re_list.size(); r++) // loop over RE
	{
	    matrix<Type> re_list_el(rnl(r), rnp(r));
	    for (int j = 0; j < rnp(r); j++) // loop over RE parameters
	    {
		re_list_el.col(j) = b.segment(k, rnl(r));
		k += rnl(r);
	    }
	    re_list(r) = re_list_el;
	}
	REPORT(re_list);


	// Form covariance matrices --------------------------------------------

	// Covariance matrix in univariate case
	matrix<Type> cov_list_el(1, 1);
	cov_list_el(0, 0) = Type(1);

	for (int r = 0; r < cov_list.size(); r++) // loop over RE
	{
	    // UNSTRUCTURED_CORR() does not tolerate `cor_ltri.size() == 0`
	    if (rnp(r) == 1)
	    {
		cov_list(r) = cov_list_el;
		continue;
	    }
	    else
	    {
		vector<Type> cor_ltri(rnp(r) * (rnp(r) - 1) / 2);
		cov_list(r) = density::UNSTRUCTURED_CORR(cor_ltri).cov();
	    }
	}
	REPORT(cov_list);


	// Form scale vectors --------------------------------------------------

	for (int r = 0, k = 0; r < sd_list.size(); r++)
	{
	    vector<Type> sd_list_el = sd_b.segment(k, rnp(r));
	    sd_list(r) = sd_list_el;
	    k += rnp(r);
	}
	REPORT(sd_list);


	// Construct `b2` ------------------------------------------------------

	if (spZ_flag == 1)
	{
	    b2_sparse.reserve(rnl);
	    for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	    {
		for (int j = 0, k = 0; j < np; j++) // loop over parameters
		{
		    if (rid(r, j) == 1)
		    {
			for (int l = 0; l < rnl(r); l++) // loop over levels
			{
			    b2_sparse.insert(i + l, j) = sd_list(r)(k) * re_list(r)(l, k);
			}
			k += 1;
		    }
		}
		i += rnl(r);
	    }
	    b2_sparse.makeCompressed();
	    DATA_SPARSE_MATRIX(Z);
	    Q += Z * b2_sparse;
	}
	else
	{
	    b2_dense.fill(Type(0));
	    for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	    {
		for (int j = 0, k = 0; j < np; j++) // loop over parameters
		{
		    if (rid(r, j) == 1)
		    {
		        b2_dense.block(i, j, rnl(r), 1) = sd_list(r)(k) * re_list(r).col(k);
			k += 1;
		    }
		}
		i += rnl(r);
	    }
	    DATA_MATRIX(Z);
	    Q += Z * b2_dense;
	}
    }


    // Compute likelihood ======================================================

    // Parameter values
    vector<Type> r(t.size());
    vector<Type> log_c0(t.size());
    vector<Type> thalf(t.size());
    vector<Type> log_K(t.size());
    vector<Type> p(t.size());
    vector<Type> log_nbdisp(t.size());
    vector<Type> log_b(t.size());
    Type log_nbexcessvar; // excess variance = mean^2/dispersion

    r = Q.col(i_log_r);
    if (!isNA(i_log_c0))
    {
        log_c0 = Q.col(i_log_c0);
    }
    else
    {
        thalf = Q.col(i_log_thalf);
        log_K = Q.col(i_log_K);
	if (!isNA(i_log_p))
	{
	    p = Q.col(i_log_p);
	}
    }
    if (!isNA(i_log_nbdisp))
    {
        log_nbdisp = Q.col(i_log_nbdisp);
    }
    if (!isNA(i_log_b))
    {
        log_b = Q.col(i_log_b);
    }
    
    // Incidence curves
    vector<Type> log_cum_inc(t.size());                // log cumulative incidence
    vector<Type> log_int_inc(t.size()-w.maxCoeff()-1); // log interval incidence

    // Objective function
    Type nll = Type(0);


    // Evaluate log cumulative incidence at observation times
    for (int i = 0; i < t.size(); i++)
    {
        // curve="exponential"
        if (!isNA(i_log_c0))
	{
	    log_cum_inc(i) = log_c0(i) + r(i) * t(i);
	}
	else
	{
	    // curve="logistic"
	    if (isNA(i_log_p))
	    {
	        log_cum_inc(i) = log_K(i) - logspace_add(Type(0), -r(i) * (t(i) - thalf(i)));
	    }
	    // curve="richards"
	    else
	    {
		log_cum_inc(i) = log_K(i) - logspace_add(Type(0), logspace_sub(p(i) * log(2), Type(0)) - r(i) * p(i) * (t(i) - thalf(i))) / p(i);
	    }
        }
        // include_baseline=TRUE
	if (!isNA(i_log_b))
	{
	    log_cum_inc(i) = logspace_add(log_b(i) + log(t(i)), log_cum_inc(i));
	}
    }

    // Compute log interval incidence and increment negative log likelihood
    for (int i = 0, d = 0; i < t.size()-1; i++)
    {
        // First element of `x` in each level of `w` is ignored
        if (w(i) != w(i+1))
	{
	    d++;
	    continue;
	}
	
        log_int_inc(i-d) = logspace_sub(log_cum_inc(i+1), log_cum_inc(i));
        if (!isNA(x(i)))
        {
	    // distr="pois"
	    if (isNA(i_log_nbdisp))
	    {
	        nll -= dpois_robust(x(i+1), log_int_inc(i-d), true);
	    }
	    // distr="nbinom"
	    else
	    {
	        log_nbexcessvar = Type(2) * log_int_inc(i-d) - log_nbdisp(i+1);
		nll -= dnbinom_robust(x(i+1), log_int_inc(i-d), log_nbexcessvar, true);
	    }
        }
    }

    // Compute parameter joint density and increment negative log likelihood
    if (anyRE)
    {
	for (int r = 0; r < re_list.size(); r++) // loop over RE
	{
	    for (int l = 0; l < rnl(r); l++) // loop over levels
	    {
	        nll += density::MVNORM(cov_list(r))(re_list(r).row(l));
	    }
	}
    }

    REPORT(Q);
    REPORT(log_cum_inc);
    REPORT(log_int_inc);


    // Predict incidence for new data ==========================================

    if (predict_flag == 1)
    {
        // Data
        // NOTE: X/Z matrices are row vectors because we are enforcing on R side
        // that new data belong to one group
        DATA_VECTOR(t_new);
        DATA_MATRIX(X_new);
	DATA_MATRIX(Z_new);

	// np-vector of parameter values
	matrix<Type> Q_new(1, np);
	if (spX_flag == 1)
	{
	    Q_new = X_new * beta2_sparse;
	}
	else
	{
	    Q_new = X_new * beta2_dense;
	}
	
	if (anyRE)
	{
	    if (spZ_flag == 1)
	    {
	        Q_new += Z_new * b2_sparse;
	    }
	    else
	    {
	        Q_new += Z_new * b2_dense;
	    }
	}
	
	// Parameter values
	Type r_new;
	Type log_c0_new;
	Type thalf_new;
	Type log_K_new;
	Type p_new;
	Type log_b_new;

	r_new = exp(Q_new(1, i_log_r));
	if (!isNA(i_log_c0))
	{
	    log_c0_new = Q_new(1, i_log_c0);
	}
	else
	{
	    thalf_new = exp(Q_new(1, i_log_thalf));
	    log_K_new = Q_new(1, i_log_K);
	    if (!isNA(i_log_p))
	    {
	        p_new = exp(Q_new(1, i_log_p));
	    }
	}
	if (!isNA(i_log_b))
	{
	    log_b_new = Q_new(1, i_log_b);
	}
	
	// Incidence curves
	vector<Type> log_cum_inc_new(t_new.size());
	vector<Type> log_int_inc_new(t_new.size()-1);

	// Evaluate log cumulative incidence at new times
	for (int i = 0; i < t_new.size(); i++)
	{
	    if (!isNA(i_log_c0)) // curve="exponential"
	    {
		log_cum_inc_new(i) = log_c0_new + r_new * t_new(i);
	    }
	    else if (isNA(i_log_p)) // curve="logistic"
	    {
		log_cum_inc_new(i) = log_K_new - logspace_add(Type(0), -r_new * (t_new(i) - thalf_new));
	    }
	    else // curve="richards"
	    {
	        log_cum_inc_new(i) = log_K_new - logspace_add(Type(0), logspace_sub(p_new * log(2), Type(0)) - r_new * p_new * (t_new(i) - thalf_new)) / p_new;
	    }
	    if (!isNA(i_log_b)) // include_baseline=TRUE
	    {
	        log_cum_inc_new(i) = logspace_add(log_b_new + log(t_new(i)), log_cum_inc_new(i));
	    }
	}

        // Compute log interval incidence
	for (int i = 0; i < t_new.size()-1; i++)
	{
	    log_int_inc_new(i) = logspace_sub(log_cum_inc_new(i+1), log_cum_inc_new(i));
        }

	REPORT(Q_new);
        ADREPORT(log_cum_inc_new);
	ADREPORT(log_int_inc_new);
    }
    

    return nll;
}
