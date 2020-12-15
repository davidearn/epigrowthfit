// #define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>

template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA_real_(Type x)
{
    return R_IsNA(asDouble(x));
}


// https://kaskr.github.io/adcomp/Introduction.html
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Flags
    DATA_INTEGER(spX_flag);     // sparse X  (1=sparse,   0=dense)
    DATA_INTEGER(spZ_flag);     // sparse Z  (1=sparse,   0=dense)
    DATA_INTEGER(predict_flag); // predict   (1=predict,  0=fit)

    // Data
    // time series
    DATA_VECTOR(t); // length=N
    DATA_VECTOR(x); // length=N
    // window
    DATA_IVECTOR(w); // length=N,  {0,...,nlevels(w)-1}
    // nlevels
    DATA_IVECTOR(fnl); // length=nFE=np
    DATA_IVECTOR(rnl); // length=nRE
    // RE indicators
    DATA_IMATRIX(rid); // nrow=nRE,  ncol=np,  {0,1}
    // parameter indices
    DATA_INTEGER(i_log_r);
    DATA_INTEGER(i_log_c0);     // !NA  iff  curve="exponential"
    DATA_INTEGER(i_log_thalf);  //  NA  iff  curve="exponential"
    DATA_INTEGER(i_log_K);      //  NA  iff  curve="exponential"
    DATA_INTEGER(i_log_p);      // !NA  iff  curve="richards"
    DATA_INTEGER(i_log_nbdisp); // !NA  iff  distr="nbinom"
    DATA_INTEGER(i_log_b);      //  NA  iff  excess=FALSE
    
    // Parameters
    // fixed effects
    PARAMETER_VECTOR(beta); // length=sum(nlevels(FE i))
    // sd random effects
    PARAMETER_VECTOR(sd_b); // length=sum(np(RE i))
    // random effects
    PARAMETER_VECTOR(b);    // length=sum(nlevels(RE i)*np(RE i)) 

    // Number of parameters
    int np = rid.cols();

    // Any random effects?
    bool anyRE = (rid.rows() > 0);

    
    // Compute parameter values at every time point ============================
    // NOTE: parameters are constant within but not across fitting windows

    // (N x np) matrix of parameter values
    matrix<Type> Y(t.size(), np);
    
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
	Y = X * beta2_sparse;
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
	Y = X * beta2_dense;
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
    // of vector `b`: rows of each matrix follow a zero-mean, unit-variance
    // multivariate normal distribution with a common covariance matrix
    vector< matrix<Type> > re_list(rid.rows());

    // A list of s.d. vectors corresponding elementwise to `re_list`
    vector< vector<Type> > sd_list(rid.rows());

    // A list of correlation matrices corresponding elementwise to `re_list`
    vector< matrix<Type> > cor_list(rid.rows());
	
    
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


	// Form scale vectors --------------------------------------------------

	for (int r = 0, k = 0; r < sd_list.size(); r++)
	{
	    vector<Type> sd_list_el = sd_b.segment(k, rnp(r));
	    sd_list(r) = sd_list_el;
	    k += rnp(r);
	}
	REPORT(sd_list);

	
	// Form correlation matrices -------------------------------------------

	// Correlation matrix in univariate case
	matrix<Type> cor_list_el(1, 1);
	cor_list_el(0, 0) = Type(1);

	for (int r = 0; r < cor_list.size(); r++) // loop over RE
	{
	    // UNSTRUCTURED_CORR() does not tolerate `ltri.size() == 0`
	    if (rnp(r) == 1)
	    {
		cor_list(r) = cor_list_el;
		continue;
	    }
	    else
	    {
		vector<Type> ltri(rnp(r) * (rnp(r) - 1) / 2);
		cor_list(r) = density::UNSTRUCTURED_CORR(ltri).cov();
	    }
	}
	REPORT(cor_list);


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
	    Y += Z * b2_sparse;
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
	    Y += Z * b2_dense;
	}
    }


    // Report parameter values in each fitting window ==========================
    // NOTE: parameters are constant within but not across fitting windows

    matrix<Type> Y_short(w.maxCoeff()+1, Y.cols());
    Y_short.row(0) = Y.row(0);

    for (int i = 1, j = 1; (i < Y.rows()) && (j < Y_short.rows()); i++)
    {
        if (w(i) == w(i-1))
	{
	    continue;
	}
	Y_short.row(j) = Y.row(i);
	j++;
    }

    vector<Type> y = Y_short.vec();
    ADREPORT(y);
    

    // Compute likelihood ======================================================

    // Parameter values
    vector<Type> log_r(t.size());
    vector<Type> log_c0(t.size());
    vector<Type> log_thalf(t.size());
    vector<Type> log_K(t.size());
    vector<Type> log_p(t.size());
    vector<Type> log_nbdisp(t.size());
    vector<Type> log_b(t.size());
    Type log_nbexcessvar; // excess variance = mean^2/dispersion

    log_r = Y.col(i_log_r);
    if (i_log_c0 >= 0) // FIXME: isNA_integer_()?
    {
        log_c0 = Y.col(i_log_c0);
    }
    else
    {
        log_thalf = Y.col(i_log_thalf);
        log_K = Y.col(i_log_K);
    	if (i_log_p >= 0)
    	{
    	    log_p = Y.col(i_log_p);
    	}
    }
    if (i_log_nbdisp >= 0)
    {
        log_nbdisp = Y.col(i_log_nbdisp);
    }
    if (i_log_b >= 0)
    {
        log_b = Y.col(i_log_b);
    }

    // Inverse-link transformed parameter values
    // FIXME: vector<Type> r = exp(Y.col(i_log_r));?
    vector<Type> r = exp(log_r);
    vector<Type> thalf = exp(log_thalf);
    vector<Type> p = exp(log_p);
    
    // Incidence curves
    vector<Type> log_cum_inc(t.size()); // log cumulative incidence
    vector<Type> log_int_inc(t.size()-w.maxCoeff()-1); // log interval incidence

    // Objective function
    Type nll = 0;


    // Evaluate log cumulative incidence at observation times
    for (int i = 0; i < t.size(); i++)
    {
        // curve="exponential"
        if (i_log_c0 >= 0)
    	{
    	    log_cum_inc(i) = log_c0(i) + r(i) * t(i);
    	}
    	else
    	{
    	    // curve="logistic"
    	    if (i_log_p < 0)
    	    {
    	        log_cum_inc(i) = log_K(i) - logspace_add(Type(0), -r(i) * (t(i) - thalf(i)));
    	    }
    	    // curve="richards"
    	    else
    	    {
    		log_cum_inc(i) = log_K(i) - logspace_add(Type(0), logspace_sub(p(i) * log(2), Type(0)) - r(i) * p(i) * (t(i) - thalf(i))) / p(i);
    	    }
        }
        // excess=TRUE
    	if (i_log_b >= 0)
    	{
    	    log_cum_inc(i) = logspace_add(log_b(i) + log(t(i)), log_cum_inc(i));
    	}
    }
    REPORT(log_cum_inc);

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
        if (!isNA_real_(x(i)))
        {
    	    // distr="pois"
    	    if (i_log_nbdisp < 0)
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
    REPORT(log_int_inc);

    // Compute parameter joint density and increment negative log likelihood
    if (anyRE)
    {
    	for (int r = 0; r < re_list.size(); r++) // loop over RE
    	{
    	    for (int l = 0; l < rnl(r); l++) // loop over levels
    	    {
    	        nll += density::MVNORM(cor_list(r))(re_list(r).row(l));
    	    }
    	}
    }
    

    // Predict incidence for new data ==========================================

    if (predict_flag == 1)
    {
        // Flags
        DATA_INTEGER(se_flag); // report  (1=ADREPORT, 0=REPORT)
      
        // Data
        // NOTE: X/Z matrices are row vectors because we are enforcing on R side
        // that new data belong to one group
        DATA_VECTOR(t_new);
        DATA_MATRIX(X_new);
    	DATA_MATRIX(Z_new);

    	// np-vector of parameter values
    	matrix<Type> Y_new(1, np);
    	if (spX_flag == 1)
    	{
    	    Y_new = X_new * beta2_sparse;
    	}
    	else
    	{
    	    Y_new = X_new * beta2_dense;
    	}
	
    	if (anyRE)
    	{
    	    if (spZ_flag == 1)
    	    {
    	        Y_new += Z_new * b2_sparse;
    	    }
    	    else
    	    {
    	        Y_new += Z_new * b2_dense;
    	    }
    	}
	
    	// Parameter values
    	Type r_new;
    	Type log_c0_new;
    	Type thalf_new;
    	Type log_K_new;
    	Type p_new;
    	Type log_b_new;

        r_new = exp(Y_new(0, i_log_r));
    	if (i_log_c0 >= 0)
    	{
    	    log_c0_new = Y_new(0, i_log_c0);
    	}
    	else
    	{
    	    thalf_new = exp(Y_new(0, i_log_thalf));
    	    log_K_new = Y_new(0, i_log_K);
    	    if (i_log_p >= 0)
    	    {
    	        p_new = exp(Y_new(0, i_log_p));
    	    }
    	}
    	if (i_log_b >= 0)
    	{
    	    log_b_new = Y_new(0, i_log_b);
    	}
	
    	// Incidence curves
    	vector<Type> log_cum_inc_new(t_new.size());
    	vector<Type> log_int_inc_new(t_new.size()-1);

    	// Evaluate log cumulative incidence at new times
    	for (int i = 0; i < t_new.size(); i++)
    	{
	    // curve="exponential"
    	    if (i_log_c0 >= 0)
    	    {
    		log_cum_inc_new(i) = log_c0_new + r_new * t_new(i);
    	    }
    	    else
    	    {
	        // curve="logistic"
	        if (i_log_p < 0)
		{
		    log_cum_inc_new(i) = log_K_new - logspace_add(Type(0), -r_new * (t_new(i) - thalf_new));
		}
		// curve="richards"
		else
		{
		    log_cum_inc_new(i) = log_K_new - logspace_add(Type(0), logspace_sub(p_new * log(2), Type(0)) - r_new * p_new * (t_new(i) - thalf_new)) / p_new;
		}
    		
    	    }
	    // excess=TRUE
    	    if (i_log_b >= 0)
    	    {
    	        log_cum_inc_new(i) = logspace_add(log_b_new + log(t_new(i)), log_cum_inc_new(i));
    	    }
    	}

        // Compute log interval incidence
    	for (int i = 0; i < t_new.size()-1; i++)
    	{
    	    log_int_inc_new(i) = logspace_sub(log_cum_inc_new(i+1), log_cum_inc_new(i));
        }

	// Report to R session
	if (se_flag == 1)
	{
	    ADREPORT(log_cum_inc_new);
	    ADREPORT(log_int_inc_new);
	}
	else
	{
	    REPORT(log_cum_inc_new);
	    REPORT(log_int_inc_new);
	}
    }
    
    return nll;
}
