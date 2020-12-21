// #define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>

enum curve
{
    exponential,
    logistic,
    richards
};

template<class Type>
Type eval_exponential(Type t, Type log_r, Type log_c0)
{
    // log(c(t))
    // = log(c0 * exp(r * t))
    // = log(c0) + r * t
    return log_c0 + exp(log_r) * t;
}

template<class Type>
Type eval_logistic(Type t, Type log_r, Type log_thalf, Type log_K)
{
    // log(c(t))
    // = log(K / (1 + exp(-r * (t - thalf))))
    // = log(K) - log(1 + exp(-r * (t - thalf)))
    return log_K - logspace_add(Type(0), -exp(log_r) * (t - exp(log_thalf)));
}

template<class Type>
Type eval_richards(Type t, Type log_r, Type log_thalf, Type log_K, Type log_p)
{
    // log(c(t))
    // = log(K / (1 + (2^p - 1) * exp(-r * p * (t - thalf)))^(1 / p))
    // = log(K) - log(1 + (2^p - 1) * exp(-r * p * (t - thalf))) / p
    Type p = exp(log_p);
    return log_K - logspace_add(Type(0), logspace_sub(p * log(Type(2)), Type(0)) - exp(log_r) * p * (t - exp(log_thalf))) / p;
}

template<class Type>
vector<Type> eval_log_cum_inc(vector<Type> t, int curve_flag, bool excess,
			      matrix<Type> Y,
			      int j_log_r,
			      int j_log_c0,
			      int j_log_thalf,
			      int j_log_K,
			      int j_log_p,
			      int j_log_b)
{
    vector<Type> log_cum_inc(t.size());
    for (int i = 0; i < t.size(); i++)
    {
        switch (curve_flag)
	{
	case exponential:
	    log_cum_inc(i) = eval_exponential(t(i), Y(i, j_log_r), Y(i, j_log_c0));
	    break;
	case logistic:
	    log_cum_inc(i) = eval_logistic(t(i), Y(i, j_log_r), Y(i, j_log_thalf), Y(i, j_log_K));
	    break;
	case richards:
	    log_cum_inc(i) = eval_richards(t(i), Y(i, j_log_r), Y(i, j_log_thalf), Y(i, j_log_K), Y(i, j_log_p));
	    break;
	}
	if (excess)
	{
	    log_cum_inc(i) = logspace_add(Y(i, j_log_b) + log(t(i)), log_cum_inc(i));
	}
    }
    return log_cum_inc;
}

template<class Type>
vector<Type> eval_log_int_inc(vector<Type> log_cum_inc, vector<int> wl)
{
    vector<Type> log_int_inc(log_cum_inc.size() - wl.size());
    for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++)
    {
        for (int j = 0; j < wl(k) - 1; j++)
	{
	    log_int_inc(i1+j) = logspace_sub(log_cum_inc(i2+j), log_cum_inc(i2+j-1));
	}
	i1 += wl(k) - 1;
	i2 += wl(k);
    }
    return log_int_inc;
}

template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

template<class Type>
Type rnbinom_robust(Type log_mu, Type log_disp)
{
    Type log_p = log_disp - logspace_add(log_mu, log_disp);
    return rnbinom(exp(log_disp), exp(log_p));
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA_real_(Type x)
{
    return R_IsNA(asDouble(x));
}

template<class Type>
matrix<Type> prune_dupl_rows(matrix<Type> Y, vector<int> wl)
{
    matrix<Type> Y_short(wl.size(), Y.cols());
    for (int i = 0, k = 0; k < wl.size(); k++)
    {
        Y_short.row(k) = Y.row(i);
	i += wl(k);
    }
    return Y_short;
}





// https://kaskr.github.io/adcomp/Introduction.html
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Flags
    DATA_INTEGER(sparseX_flag); // sparse X  (1=sparse,   0=dense)
    DATA_INTEGER(sparseZ_flag); // sparse Z  (1=sparse,   0=dense)
    DATA_INTEGER(predict_flag); // predict   (1=predict,  0=fit)
    bool sparseX = (sparseX_flag == 1);
    bool sparseZ = (sparseZ_flag == 1);
    bool predict = (predict_flag == 1);
    
    // Data
    // time series
    DATA_VECTOR(t);  // length=N
    DATA_VECTOR(x);  // length=N
    // window lengths
    DATA_IVECTOR(wl); // length=nts
    int nw = wl.size();
    // factor nlevels
    DATA_IVECTOR(fnl); // length=nFE=np
    DATA_IVECTOR(rnl); // length=nRE
    // design matrices
    Eigen::SparseMatrix<Type> Xs(t.size(), fnl.sum());
    matrix<Type> Xd(t.size(), fnl.sum());
    if (sparseX)
    {
        DATA_SPARSE_MATRIX(X); // nrow=N,  ncol=sum(nlevels(FE i))
	Xs = X;
    }
    else
    {
	DATA_MATRIX(X);
	Xd = X;
    }
    Eigen::SparseMatrix<Type> Zs(t.size(), rnl.sum());
    matrix<Type> Zd(t.size(), rnl.sum());
    if (sparseZ)
    {
        DATA_SPARSE_MATRIX(Z); // nrow=N, ncol=sum(nlevels(RE i))
	Zs = Z;
    }
    else
    {
	DATA_MATRIX(Z);
	Zd = Z;
    }
    // random effect indicators
    DATA_IMATRIX(rid); // nrow=nRE,  ncol=np,  {0,1}
    bool anyRE = (rid.rows() > 0);
    int np = rid.cols();
    vector<int> rnp = rid.rowwise().sum();
    // parameter indices
    DATA_INTEGER(j_log_r);
    DATA_INTEGER(j_log_c0);
    DATA_INTEGER(j_log_thalf);
    DATA_INTEGER(j_log_K);
    DATA_INTEGER(j_log_p);
    DATA_INTEGER(j_log_nbdisp);
    DATA_INTEGER(j_log_b);
    int curve_flag = ((j_log_c0 >= 0) ? 0 : ((j_log_p < 0) ? 1 : 2));
    bool nbinom = (j_log_nbdisp >= 0);
    bool excess = (j_log_b >= 0);
    
    // Parameters
    // fixed effects
    PARAMETER_VECTOR(beta); // length=sum(nlevels(FE i))
    // sd random effects
    PARAMETER_VECTOR(sd); // length=sum(np(RE i))
    // random effects
    PARAMETER_VECTOR(b);    // length=sum(nlevels(RE i)*np(RE i)) 

    
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

    Eigen::SparseMatrix<Type> beta2(fnl.sum(), np);
    beta2.reserve(fnl);
    
    for (int j = 0, i = 0; j < np; j++) // loop over parameters
    {
	for (int l = 0; l < fnl(j); l++) // loop over levels
	{
	    beta2.insert(i + l, j) = beta(i + l);
	}
	i += fnl(j);
    }
    
    if (sparseX)
    {
        Y = Xs * beta2;
    }
    else
    {
        Y = Xd * beta2;
    }
    matrix<Type> Y_sim = Y; // preserve fixed effects component for SIMULATE
    
    // Random effects component:
    // matrix multiply
    // Z * b2
    // = (N x sum(nlevels(RE i))) * (sum(nlevels(RE i)) x np)
    // = (N x np)
    // need to construct block matrix `b2` from vector `b`
    // will require some setting up...
    
    matrix<Type> b2(rnl.sum(), np);
    
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

	for (int r = 0, k = 0; r < re_list.size(); r++) // loop over RE
	{
	    matrix<Type> re_list_el(rnl(r), rnp(r));
	    for (int j = 0; j < rnp(r); j++) // loop over parameters
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
	    vector<Type> sd_list_el = sd.segment(k, rnp(r));
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

        b2.fill(Type(0));
	for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	{
	    for (int j = 0, k = 0; j < np; j++) // loop over parameters
	    {
		if (rid(r, j) == 1)
		{
		    b2.block(i, j, rnl(r), 1) = sd_list(r)(k) * re_list(r).col(k);
		    k++;
		}
	    }
	    i += rnl(r);
	}
	
	if (sparseZ)
	{
	    Y += Zs * b2;
	}
	else
	{
	    Y += Zd * b2;
	}
    }


    // Report parameter values in each fitting window ==========================
    // NOTE: parameters are constant within but not across fitting windows

    matrix<Type> Y_short = prune_dupl_rows(Y, wl);
    vector<Type> y = Y_short.vec();
    ADREPORT(y);
    

    // Compute likelihood ======================================================

    // Log cumulative incidence
    vector<Type> log_cum_inc = eval_log_cum_inc(t, curve_flag, excess, Y, j_log_r, j_log_c0, j_log_thalf, j_log_K, j_log_p, j_log_b);
    
    // Log interval incidence
    vector<Type> log_int_inc = eval_log_int_inc(log_cum_inc, wl);

    // Negative log likelihood
    Type nll = 0;
    Type log_var_minus_mu;
    
    for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++) // loop over time series
    {
        for (int j = 0; j < wl(k) - 1; j++) // loop over time points
	{
	    if (!isNA_real_(x(i2+j)))
	    {
	        if (nbinom)
		{
		    log_var_minus_mu = Type(2) * log_int_inc(i1+j) - Y(i2+j, j_log_nbdisp);
		    nll -= dnbinom_robust(x(i2+j), log_int_inc(i1+j), log_var_minus_mu, true);
		    // usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)
		}
		else
		{
		    nll -= dpois_robust(x(i2+j), log_int_inc(i1+j), true);
		    // usage: dpois_robust(x, log_lambda, give_log)
		}
	    }
	}
	i1 += wl(k) - 1;
	i2 += wl(k);
    }
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


    // Simulate incidence ======================================================

    SIMULATE
    {
        if (anyRE)
	{
	    matrix<Type> b2_sim(rnl.sum(), np);
	    b2_sim.fill(Type(0));
	    
	    for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	    {
		vector<Type> z(sd_list(r).size());
		for (int l = 0; l < rnl(r); l++) // loop over levels
		{
		    z = sd_list(r) * density::MVNORM(cor_list(r)).simulate();
		    for (int j = 0, k = 0; j < np; j++) // loop over parameters
		    {
			if (rid(r, j) == 1)
			{
			    b2_sim(i + l, j) = z(k);
			    k++;
			}
		    }
		}
		i += rnl(r);
	    }
	    
	    if (sparseZ)
	    {
		Y_sim += Zs * b2_sim;
	    }
	    else
	    {
		Y_sim += Zd * b2_sim;
	    }
	}
	matrix<Type> Y_sim_short = prune_dupl_rows(Y_sim, wl);
	REPORT(Y_sim_short);
	
        log_cum_inc = eval_log_cum_inc(t, curve_flag, excess, Y_sim, j_log_r, j_log_c0, j_log_thalf, j_log_K, j_log_p, j_log_b);
	log_int_inc = eval_log_int_inc(log_cum_inc, wl);

	for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++) // loop over time series
	{
	    for (int j = 0; j < wl(k) - 1; j++) // loop over time points
	    {
		if (nbinom)
		{
		    x(i2+j) = rnbinom_robust(log_int_inc(i1+j), Y_sim(i2+j, j_log_nbdisp));
		    // usage: rnbinom_robust(log_mu, log_disp)
		}
		else
		{
		    x(i2+j) = rpois(exp(log_int_inc(i1+j)));
		    // usage: rpois(mu)
		}

	    }
	    i1 += wl(k) - 1;
	    i2 += wl(k);
	}
	REPORT(x);
    }
    

    // Predict incidence =======================================================
    // NOTE: enforcing on R side that there is only one time series here
    
    if (predict)
    {
        DATA_VECTOR(t_new); // length=N_new
	vector<int> wl_new(1);
	wl_new(0) = t_new.size();
	DATA_MATRIX(X_new); // nrow=N_new,  ncol=sum(nlevels(FE i))
    	DATA_MATRIX(Z_new); // nrow=N_new,  ncol=sum(nlevels(RE i))
	DATA_INTEGER(se_flag); // report  (1=ADREPORT, 0=REPORT)
	bool se = (se_flag == 1);
      
        matrix<Type> Y_new = X_new * beta2;
    	if (anyRE)
    	{
	    Y_new += Z_new * b2;
    	}

        vector<Type> log_cum_inc_new = eval_log_cum_inc(t_new, curve_flag, excess, Y_new, j_log_r, j_log_c0, j_log_thalf, j_log_K, j_log_p, j_log_b);
    	vector<Type> log_int_inc_new = eval_log_int_inc(log_cum_inc_new, wl_new);

    	if (se)
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
